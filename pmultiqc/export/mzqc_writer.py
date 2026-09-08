import datetime
import json
import logging
import os
from typing import Any

from pmultiqc.export.mapping import MZQC_METRIC_MAPPING

logger = logging.getLogger(__name__)

class MzQCExporter:
    def __init__(self, version: str = "1.0.0"):
        self.version = version
        self.cv_ref = {
            "name": "PSI-MS",
            "uri": "https://raw.githubusercontent.com/HUPO-PSI/psi-ms-CV/master/psi-ms-CV.obo",
            "version": "4.1.140"
        }

    def _create_input_file_meta(self, file_name: str) -> dict[str, Any]:
        return {
            "name": file_name,
            "location": f"file:///{file_name}",
            "fileFormat": {"accession": "MS:1000584", "name": "mzML format"}
        }

    def _format_value(self, value: Any) -> Any:
        if isinstance(value, (int, float, str, list, dict)):
            return value
        return str(value)

    def create_metric_entry(self, internal_key: str, raw_data: Any) -> dict[str, Any]:
        mapping = MZQC_METRIC_MAPPING.get(internal_key)
        if not mapping:
            raise ValueError(f"Metric key '{internal_key}' is not defined in mappings.py")
        return {
            "accession": mapping["accession"],
            "name": mapping["name"],
            "value": self._format_value(raw_data)
        }

    def generate_mzqc(self, run_name: str, run_metrics: dict[str, Any], set_metrics: dict[str, Any] = None) -> dict[str, Any]:
        input_file_meta = self._create_input_file_meta(run_name)
        run_quality_list = []
        
        if run_metrics:
            quality_metrics = []
            for key, data in run_metrics.items():
                if key in MZQC_METRIC_MAPPING and MZQC_METRIC_MAPPING[key]["scope"] == "run":
                    quality_metrics.append(self.create_metric_entry(key, data))
            
            run_quality_list.append({
                "metadata": {
                    "inputFiles": [input_file_meta],
                    "analysisSoftware": [{"accession": "MS:1003160", "name": "MultiQC", "version": "1.0"}]
                },
                "qualityMetrics": quality_metrics
            })

        set_quality_list = []
        if set_metrics:
            quality_metrics = []
            for key, data in set_metrics.items():
                if key in MZQC_METRIC_MAPPING and MZQC_METRIC_MAPPING[key]["scope"] == "set":
                    quality_metrics.append(self.create_metric_entry(key, data))
            if quality_metrics:
                set_quality_list.append({
                    "metadata": {
                        "inputFiles": [input_file_meta],
                        "analysisSoftware": [{"accession": "MS:1003160", "name": "MultiQC", "version": "1.0"}]
                    },
                    "qualityMetrics": quality_metrics
                })

        return {
            "mzQC": {
                "version": self.version,
                "creationDate": datetime.datetime.now().isoformat() + "Z",
                "controlledVocabularies": [self.cv_ref],
                "runQuality": run_quality_list,
                "setQuality": set_quality_list
            }
        }
logger = logging.getLogger(__name__)
def extract_pca_coordinate_data(plot_id, plot_obj):
            datasets = plot_obj.get("datasets", [])
            if not datasets or not isinstance(datasets, list):
                return {}
                
            first_ds = datasets[0]
            points_list = first_ds.get("points", [])
            
            if not isinstance(points_list, list) or not points_list:
                return {}
                
            extracted_samples = {}
            for pt in points_list:
                if isinstance(pt, dict):
                    sample_name = pt.get("name") or pt.get("label") or "Unknown_Sample"
                    
                    # 1. Dynamic Dictionary Parsing
                    # If MultiQC stores them as explicit keys ('x', 'y', 'z' or 'PC1', 'PC2')
                    if "PC1" in pt:
                        # Grab all present PC keys dynamically
                        pc_coords = {k: v for k, v in pt.items() if k.startswith("PC")}
                    else:
                        # Default fallback for standard standard X/Y scatter objects
                        pc_coords = {"PC1": pt.get("x"), "PC2": pt.get("y")}
                        if "z" in pt:
                            pc_coords["PC3"] = pt.get("z")

                elif isinstance(pt, list) and len(pt) >= 2:
                    # 2. Dynamic List Parsing (e.g., ['Sample_A', 32.85, -23.46, 5.12])
                    #  first element always sample name string
                    sample_name = pt[0]
                    coordinate_array = pt[1:]
                    
                    # Generate "PC1", "PC2", "PC3" dynamically based on array length
                    pc_coords = {}
                    for idx, val in enumerate(coordinate_array):
                        pc_coords[f"PC{idx + 1}"] = val
                else:
                    continue

                if sample_name and pc_coords:
                    extracted_samples[sample_name] = pc_coords
                    
            return extracted_samples


def write_mzqc(output_dir: str, input_dir: str, job_id: str, workflow: str, json_path: str = None, ) -> None:
    """
    Parses MultiQC's report_plot_data structure to extract  metrics
    and generates mzQC output file.
    """
    logger.info("Extracting metrics from MultiQCs cache.")
    
    if json_path is not None:
        multiqc_json_path = json_path
    else:
        multiqc_json_path = os.path.join(output_dir, "multiqc_report_data", "multiqc_data.json")
    
    if not os.path.exists(multiqc_json_path):
        logger.warning(f"MultiQC data not found at {multiqc_json_path}. Skipping.")
        return

    run_metrics = {}
    set_metrics = {}
    #run_name = workflow
    if not os.path.exists(multiqc_json_path):
        logger.warning(f"MultiQC data not found at {multiqc_json_path}. Skipping.")
        return

    if not os.path.exists(multiqc_json_path):
            logger.warning(f"MultiQC data not found at {multiqc_json_path}. Skipping.")
            return

        # --- UPDATED TO MATCH EXPECTED TEST PATTERN ---
        # This guarantees the file path matches what test.py is trying to assert
    target_mzqc_file = f"{multiqc_json_path}_quality.mzQC"
    try:
        with open(multiqc_json_path, "r") as f:
            mq_data = json.load(f)
            
        plot_data = mq_data.get("report_plot_data", {})
        
        def extract_categories_data(plot_id):
            plot_obj = plot_data.get(plot_id, {})
            datasets = plot_obj.get("datasets", [])
            
            if not datasets or not isinstance(datasets, list):
                return None
                
            first_dataset = datasets[0]
            samples = first_dataset.get("samples", [])
            cats_list = first_dataset.get("cats", [])
            
            # Rebuild a structured dictionary: { Sample_Name: { Category_Name: Value } }
            extracted_map = {}
            for cat_entry in cats_list:
                cat_name = cat_entry.get("name", "Value")
                values_array = cat_entry.get("data", [])
                
                for idx, sample_name in enumerate(samples):
                    if idx < len(values_array):
                        if sample_name not in extracted_map:
                            extracted_map[sample_name] = {}
                        extracted_map[sample_name][cat_name] = values_array[idx]
                        
            return extracted_map if extracted_map else None

        # --- Extract targeted metrics using the verified path ---
        from pmultiqc.export.mapping import MZQC_METRIC_MAPPING

        # --- Targeted Metric Execution - test -- TO DO - > external file
        target_metrics = [
                    # --- working ---
            "peptide_id_count",          
            "protein_group_count",      
            "msms_identified",          
            "charge_state",             
            "missed_cleavages",          
            "pca_of_raw_intensity",      
            "pca_of_lfq_intensity",      
                ##testing
            "ms1_ion_injection_time_mean",    
            "ms2_ion_injection_time_mean",      
            "ms1_ion_injection_time_quantiles", 
            "ms2_ion_injection_time_quantiles",
            "xic_fwhm_quantiles",               
            "xic_50_fraction",                  
            "ms2_count",                        
            "msms_identification_rate",         
            "best_andromeda_score",             
            "top_contaminants_per_raw_file",    
            "potential_contaminants_per_group", 
            "mbr_transferred_exclusive"         
        ]

        run_metrics = {}
        set_metrics = {}

        for metric_id in target_metrics:
            plot_obj = plot_data.get(metric_id, {})
            if not plot_obj:
                continue

            # Parse the metrics out with handlers
            if "pca_of_" in metric_id:
                parsed_data = extract_pca_coordinate_data(metric_id, plot_obj)
                metric_type_label = "PCA coordinates"
            else:
                parsed_data = extract_categories_data(metric_id)
                metric_type_label = "bar plot data"

            #  route to run_metrics or set_metrics based on metadata scope
            if parsed_data:
                metric_mapping = MZQC_METRIC_MAPPING.get(metric_id, {})
                scope = metric_mapping.get("scope", "run")  # Default back to run if undefined
                
                if scope == "set":
                    set_metrics[metric_id] = parsed_data
                    logger.info(f" routed to SET metrics: {metric_id} ({metric_type_label})")
                else:
                    run_metrics[metric_id] = parsed_data
                    logger.info(f" routed to RUN metrics: {metric_id} ({metric_type_label})")

    
        exporter = MzQCExporter()
    
        mzqc_document = exporter.generate_mzqc(
            run_name=job_id,
            run_metrics=run_metrics,
            set_metrics=set_metrics
        )

        # force output path 
        target_mzqc_file = os.path.join(output_dir, f"{job_id}.mzQC")

        with open(target_mzqc_file, "w") as f:
            json.dump(mzqc_document, f, indent=2)
            
        logger.info(f"Successfully generated mzQC report : {target_mzqc_file}")

    except Exception as e:
        logger.error(f"Critical error during mzQC export : {e}", exc_info=True)