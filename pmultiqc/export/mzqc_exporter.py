

import json
import logging
import os
from collections.abc import Callable
from datetime import datetime
from typing import Any

import numpy as np
from mzqc.MZQCFile import QualityMetric

from .mapping import MZQC_METRIC_MAPPING

log = logging.getLogger("pmultiqc.mzqc_generator")

class MzQcExporter:
    def __init__(self, pipeline_name: str, raw_data: dict[str, Any], output_dir: str):
        self.pipeline_name = pipeline_name.lower().replace("-", "")
        self.raw_data = raw_data
        self.output_dir = output_dir

    def sanitize_value(self,v):
            class_name = v.__class__.__name__
        
            if class_name == "DataFrame":
                # Replace all NaN, inf, and -inf values with None so they become null in JSON
                clean_df = v.replace({np.nan: None, np.inf: None, -np.inf: None})
                return clean_df.to_dict(orient="records")
                
            if class_name == "Series":
                clean_series = v.replace({np.nan: None, np.inf: None, -np.inf: None})
                return clean_series.to_list()
                
            if isinstance(v, (np.ndarray, np.generic)):
                # If it's a raw numpy array, convert it to a list and handle NaNs safely
                return [None if (isinstance(x, float) and np.isnan(x)) else x for x in v.tolist()]
            return v

    ###
    def export_to_file(self, metrics: list[Any], filename: str = "qc_metrics.mzQC") -> str:
        """
        Serializes parsed QualityMetric objects into a valid mzQC JSON format 
        and writes it to the designated output directory.
        """
        if not metrics:
            self.log_warning("No metrics provided to export.")
            return ""
    
            
        # Construct the mzQC JSON schema
        mzqc_data = {
            "mzQC": {
                "version": "1.0.0",
                "creationDate": datetime.now().isoformat(),
                "controlledVocabularies": [
                    {
                        "name": "PSI-MS Quality Control Ontology",
                        "uri": "https://raw.githubusercontent.com/HUPO-PSI/psi-ms-CV/master/psi-ms.obo",
                        "version": "4.1.130"
                    }
                ],
                "runQuality": [
                    {
                        "metadata": {
                            "inputFiles": [
                                {
                                    "name": "MaxQuant Output Directory",
                                    "fileFormat": {
                                        "accession": "MS:1000584",
                                        "name": "mzML file"  # Or your specific raw format context
                                    },
                                    "location": "local://run_data"
                                }
                            ],
                            "analysisSoftware": [
                                {
                                    "accession": "MS:1002298",
                                    "name": "MaxQuant",
                                    "version": "1.6.x (or later)",
                                    "uri": "https://www.maxquant.org"
                                }
                            ]
                        },
                        "qualityMetrics": [
                            {
                                "accession": metric.accession,
                                "name": metric.name,
                                "value": self.sanitize_value(metric.value)  
                            }
                            for metric in metrics
                        ]
                    }
                ]
            }
        }

        # Determine the target export path
        target_path = os.path.join(self.output_dir, filename)
        
        try:
            # Ensure output directory exists
            os.makedirs(self.output_dir, exist_ok=True)
            
            with open(target_path, "w", encoding="utf-8") as f:
                json.dump(mzqc_data, f, indent=2)
                
            print(f"[mzQC Export] Successfully wrote mzQC file to: {target_path}")
            return target_path
            
        except Exception as e:
            raise OSError(f"Failed writing mzQC payload to disk: {e}")
    ###
    ### SHARED CALCULATION across all 7 pipelines
    ###
    
    def _extract_missed_cleavages(self, data_source: Any, key_name: str = "Missed cleavages") -> list[dict[str, Any]] | None:
        """Calculates missed cleavages and formats as a PSI-compliant table."""
        if not data_source:
            return None
            
        missed_counts = {"0": 0, "1": 0, "2": 0, "3+": 0}
        total = 0
        
        # Handle your data structure
        items = data_source.values() if isinstance(data_source, dict) else data_source
        for item in items:
            mc = item.get(key_name)
            if mc is not None:
                total += 1
                mc_str = str(mc)
                key = mc_str if mc_str in ["0", "1", "2"] else "3+"
                missed_counts[key] += 1
        
        if total == 0:
            return None

        # Transform to table format [ {col1: val, col2: val}, ... ]
        table = [
            {"number of missed cleavages": k, "fraction": round(v / total, 4)}
            for k, v in missed_counts.items()
            if v > 0  # Optional: skip rows with 0 count
        ]
        
        return table

    def _extract_charge_state_distribution(self, data_source: Any, key_name: str = "Charge") -> dict[str, float] | None:
        """Generic calculator for charge state distributions fractions."""
        if not data_source:
            return None
        charge_counts = {}
        total = 0
        items = data_source.values() if isinstance(data_source, dict) else data_source
        for item in items:
            charge = item.get(key_name)
            if charge is not None:
                total += 1
                c_str = f"z={charge}"
                charge_counts[c_str] = charge_counts.get(c_str, 0) + 1
        return {k: round(v / total, 4) for k, v in charge_counts.items()} if total > 0 else None

    ###
    ### dispatch routing
    ###

    def _execute_extraction_loop(self, registry: dict[str, Any]) -> list[QualityMetric]:
        metrics = []
        for metric_key, value in registry.items():
            
            # Skip None 
            if value is None:
                continue 
            try:
                # Sanitize the value (Ensure sanitize_value is @staticmethod)
                # This converts DataFrames/Series to list/dict so they aren't null
                clean_value = self.sanitize_value(value)
                
                # 3. Retrieve config OR use a fallback if not in mapping
                config = MZQC_METRIC_MAPPING.get(metric_key, {
                    "accession": "MS:1000000", 
                    "name": metric_key.replace("_", " ")
                })
                
                # 4. Append to metrics (NO () here!)
                metrics.append(QualityMetric(
                    accession=config["accession"],
                    name=config["name"],
                    value=clean_value
                ))
                
            except Exception as error:
                log.warning(f"Metric '{metric_key}' failed: {error}")
                
        return metrics
    ##
    ### workflow parsers
    ###

    def _parse_maxquant(self) -> list[QualityMetric]:
        """
        Parses MaxQuant results from raw data hooks into an mzqc compatible metrics list.
        """
        protein_groups = self.raw_data.get("get_protegroups_dicts") or {}
        evidence = self.raw_data.get("get_evidence_dicts") or {}
        msms_scans = self.raw_data.get("get_msms_scans_dicts") or {}
        parameters = self.raw_data.get("get_parameter_dicts") or {}
        summary_data = self.raw_data.get("ms_ms_identified") or {}
        heatmap_data = self.raw_data.get("maxquant_heatmap") or {}
        pg_summary = protein_groups.get("protein_summary") or {}
        evidence_data = evidence.get("evidence_dicts") or {}
    
        
        maxquant_registry: dict[str, Callable[[], Any]] = {
            # --- PROTEIN GROUPS FILE METRICS ---
            "protein_group_count": len(pg_summary) if pg_summary else None,
            "potential_contaminants_per_group": protein_groups.get("pg_contaminant"),
            "intensity_distribution_box": protein_groups.get("pg_intensity_distri"),
            "lfq_intensity_distribution_box": protein_groups.get("pg_lfq_intensity_distri"),
            "pca_of_raw_intensity": protein_groups.get("raw_intensity_pca"),
            "pca_of_lfq_intensity": protein_groups.get("lfq_intensity_pca"),
            "number_of_peptides_per_proteins": protein_groups.get("num_pep_per_protein_dict"),
            "protein_quant_result": pg_summary if pg_summary else None,

            # --- EVIDENCE FILE METRICS ---
            "missed_cleavages": self._extract_missed_cleavages(evidence_data, "Missed cleavages"),
            "charge_state": evidence.get("charge_counts"),
            "delta_mass_da": evidence.get("maxquant_delta_mass_da"),
            "delta_mass_ppm": evidence.get("calibrated_mass_error"),
            "uncalibrated_mass_error_box": evidence.get("uncalibrated_mass_error"),
            "peptide_intensity_distribution_box": evidence.get("peptide_intensity"),
            "peptide_length_distribution": evidence.get("peptide_length"),
            "peptides_quantification_table": evidence.get("peptides_quant_table"),
            "top_contaminants_per_raw_file": evidence.get("top_contaminants"),
            "peptide_id_count": evidence.get("peptide_id_count"),

            # --- SUMMARY FILE METRICS ---
            "msms_identified": list(summary_data.keys()) if summary_data else None,
            "identification_summary_table": summary_data if summary_data else None,

            # --- MSMS SCANS METRICS ---
            "best_andromeda_score": msms_scans.get("best_andromeda_score"),
            "summary_of_andromeda_scores": msms_scans.get("summary_of_andromeda_scores"),

            # --- CHROMATOGRAPHY & HARDWARE METRICS ---
            "IDs_over_RT": evidence.get("rt_counts"),
            "peak_width_over_RT": evidence.get("peak_rt"),
            "oversampling_distribution": evidence.get("oversampling"),

            # --- PARAMETERS & VISUALIZATIONS ---
            "parameters": parameters.get("parameters_tb_dict"),
            "heatmap": heatmap_data if heatmap_data else None,
        }
        return self._execute_extraction_loop(maxquant_registry)

    def _parse_fragpipe(self) -> list[QualityMetric]:
        # FragPipe populates a 'pipeline_stats' dict and 'charge_states' directly 
        pipe_stats = self.raw_data.get("pipeline_stats", {}) or {}
        charge_data = self.raw_data.get("charge_states", {}) or {}
        
        fragpipe_registry = {
            "peptide_id_count": lambda: pipe_stats.get("peptide_count"),
            "charge_state": lambda: charge_data if charge_data else None,
            "missed_cleavages": lambda: self.raw_data.get("missed_cleavages"),
        }
        return self._execute_extraction_loop(fragpipe_registry)

    def _parse_diann(self, diann_payload: dict[str, Any]) -> list[QualityMetric]:
        """
        Parses DIA-NN results from the DiannModule payload into an 
        mzqc compatible metrics list.
        """
     
        diann_registry = {
            # --- Identification ---
            "peptide_id_count": diann_payload.get("peptide_id_count"),
            "protein_group_count": diann_payload.get("protein_group_count"),
            "psm_count": diann_payload.get("psm_count"),
            
            # --- Distributions ---
            "peptide_length_distribution": diann_payload.get("peptide_length_data"),
            "number_of_peptides_per_proteins": diann_payload.get("peptide_per_protein_data"),
            
            # --- Quality & Trends ---
            "modifications": diann_payload.get("modifications"),
            "IDs_over_RT": diann_payload.get("rt_trends"),
            
            # --- MS1/MS2 Stats 
            "ms1_scan_count": diann_payload.get("ms1_scan_count"),
            "ms2_scan_count": diann_payload.get("ms2_scan_count"),
            
            # --- Summary Tables ---
            "identification_summary_table": diann_payload.get("identification_summary")
        }

       
        return self._execute_extraction_loop(diann_registry)
    

    def _parse_quantms(self, quantms_payload: dict[str, Any]) -> list[QualityMetric]:
        quantms_registry = {
            "missed_cleavages": quantms_payload.get("missed_cleavages"),
            "modifications": quantms_payload.get("modifications"),
            "identified_msms_spectra": quantms_payload.get("identified_msms_spectra"),
            "charge_states": quantms_payload.get("charge_states"),
            "ids_over_rt": quantms_payload.get("ids_over_rt"),
            
            "peptide_intensity": quantms_payload.get("peptide_intensity"),
            "contaminant_percent": quantms_payload.get("contaminant_percent"),
            "top_contaminants": quantms_payload.get("top_contaminants"),
            "mass_error": quantms_payload.get("mass_error"),
            "peptide_lengths": quantms_payload.get("peptide_lengths"),
            
            "peaks_per_ms2": quantms_payload.get("peaks_per_ms2"),
            "peak_distribution": quantms_payload.get("peak_distribution"),
            "charge_distribution": quantms_payload.get("charge_distribution"),
            
            "ms1_tic": quantms_payload.get("ms1_tic"),
            "ms1_bpc": quantms_payload.get("ms1_bpc")
        }

        # Route through the extraction loop
        return self._execute_extraction_loop(quantms_registry)
    
    def _parse_mzidentml(self, mzid_payload: dict[str, Any]) -> list[QualityMetric]:
        """
        Parses mzIdentML results from the MzIdentMLModule payload 
        into an mzQC compatible metrics list.
        """
       
        mzid_registry = {
            # --- Identification Stats ---
            "protein_group_count": mzid_payload.get("protein_group_count") or self.raw_data.get("total_protein_identified"),
            "total_peptide_ids": mzid_payload.get("total_peptide_ids") or self.raw_data.get("total_peptide_count"),
            "total_psm_count": mzid_payload.get("total_psm_count") or self.raw_data.get("total_ms2_spectra_identified"),
            
            # --- Quality Metrics ---
            "missed_cleavages": mzid_payload.get("missed_cleavages") or self.raw_data.get("quantms_missed_cleavages"),
            "fdr_threshold_applied": mzid_payload.get("fdr_threshold_applied"),
            "target_decoy_ratio": mzid_payload.get("target_decoy_ratio"),
            
            # --- Distributions ---
            "peptide_length_distribution": mzid_payload.get("peptide_length_distribution") or self.raw_data.get("peptide_length_distribution"),
            "peptides_per_protein_distribution": mzid_payload.get("peptides_per_protein_distribution") or self.raw_data.get("peptides_per_protein_distribution"),
            "charge_state": mzid_payload.get("charge_state") or self.raw_data.get("heatmap_charge_score"),
            
            # --- Chromatography & Instrument Info ---
            "identified_peptides_over_rt": mzid_payload.get("identified_peptides_over_rt") or self.raw_data.get("id_rt_score"),
            "peptide_modification_frequency": mzid_payload.get("peptide_modification_frequency") or self.raw_data.get("quantms_modified"),
        }

        return self._execute_extraction_loop(mzid_registry)