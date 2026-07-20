###Mapping metrics to https://github.com/HUPO-PSI/psi-ms-CV/blob/master/psi-ms.obo

MZQC_METRIC_MAPPING = {
   ####SCOPE RUN
    "missed_cleavages": {
        "accession": "MS:4000215", 
        "name": "missed cleavages fractions",
        "scope": "run",
        
    },
    "protein_group_count": {
        "accession": "MS:1003327",
        "name": "number of identified protein groups",
        "scope": "run",
       
    },
    "charge_state": {
        "accession": "MS:4000208",
        "name": "identified quantification data points charges fractions",
        "scope": "run",
        
    },
    "top_contaminants_per_raw_file": {
        "accession": "MS:4000177",
        "name": "contaminant protein abundance fraction",
        "scope": "run",
        "value_type": "float"  # Expects a value between 0.0 and 1.0
    },
    "ms2_count": {
        "accession": "MS:4000060",
        "name": "number of MS2 spectra",
        "scope": "run",
        "value_type": "int"  
    },
        "delta_mass_da": {
        "accession": "MS:1001976",
        "name": "delta M",
        "scope": "run",
        "value_type": "float",
        "unit": "UO:0000221"
    },
    "delta_mass_ppm": {
        "accession": "MS:1001976",
        "name": "delta M",
        "scope": "run",
        "value_type": "float",
        "unit": "UO:0000169"
    },
    
    ####MBR gain (transferred exclusiv)/((genuine(exclusive)+genuine+transferred))
    "mbr_transferred_exclusive": {
        "accession": "MS:4000096",
        "name": "number of identified peptides via matching",
        "scope": "run",
        "value_type": "int"
    },
    "total_peptide_ids": {
        "accession": "MS:4000095",
        "name": "number of identified peptides",
        "scope": "run",
        "value_type": "int"
    },
    ##
    "unique_peptide_counts": {
        "accession": "MS:1001897",
        "name": "MaxQuant:peptide counts (unique)",
        "scope": "run",
        "value_type": "int"
    },
    "razor_unique_peptide_counts": {
        "accession": "MS:1001899",
        "name": "MaxQuant:peptide counts (razor+unique)",
        "scope": "run",
        "value_type": "int"
    },
    "best_andromeda_score": {
        "accession": "MS:1002338",
        "name": "Andromeda:score",
        "scope": "run",
        "value_type": "float"
    },

   
    ####---
    #### ion_injection_time_over_rt
    "ms1_ion_injection_time_mean": {
        "accession": "MS:4000132",
        "name": "MS1 ion collection time mean",
        "scope": "run",
        "value_type": "float"
        
    },
    "ms2_ion_injection_time_mean": {
        "accession": "MS:4000137",
        "name": "MS2 ion collection time mean",
        "scope": "run",
        "value_type": "float"
    },
    "ms1_ion_injection_time_quantiles": {
        "accession": "MS:4000131",
        "name": "MS1 ion collection time quantiles",
        "scope": "run",
        "value_type": "list"
    },
    "xic_fwhm_quantiles": {
        "accession": "MS:4000051",
        "name": "XIC-FWHM quantiles",
        "scope": "run",
        "value_type": "list"  
    },
    "xic_50_fraction": {
        "accession": "MS:4000050",
        "name": "XIC50 fraction",
        "scope": "run",
        "value_type": "float" 
    },
     "identified_peptides_over_rt": {
        "accession": "MS:4000216",
        "name": "identified peptides over retention time",
        "scope": "run",
        "value_type": "list"
    },

    "peptide_modification_frequency": {
        "accession": "MS:4000217",
        "name": "identified peptide modifications frequency",
        "scope": "run",
        "value_type": "table"
    },

    "peptides_per_protein_distribution": {
        "accession": "MS:4000218",
        "name": "number of peptides per protein distribution",
        "scope": "run",
        "value_type": "list"
    },

    "ms2_trigger_oversampling_distribution": {
        "accession": "MS:4000219",
        "name": "MS2 trigger oversampling distribution",
        "scope": "run",
        "value_type": "list"
    },

    "chromatographic_peak_width_over_rt": {
        "accession": "MS:4000220",
        "name": "chromatographic peak width over retention time",
        "scope": "run",
        "value_type": "list",
        "unit": "UO:0000010"
    },

    "total_identified_peptides": {
        "accession": "MS:4000221",
        "name": "total number of identified peptides",
        "scope": "run",
        "value_type": "int"
    },
    "peptide_length_distribution": {
        "accession": "MS:4000222",
        "name": "peptide amino acid length distribution",
        "scope": "run",
        "value_type": "list"
    },
    "top_n_acquisition_setting": {
        "accession": "MS:4000224",
        "name": "data-dependent top N acquisition setting",
        "scope": "run",
        "value_type": "int"
    },
    "ms2_scans_per_duty_cycle_over_rt": {
        "accession": "MS:4000225",
        "name": "MS2 scans triggered per duty cycle over retention time",
        "scope": "run",
        "value_type": "list"
    },
    "uncalibrated_precursor_mass_error_distribution": {
        "accession": "MS:4000226",
        "name": "uncalibrated precursor mass error distribution",
        "scope": "run",
        "value_type": "list"
    },
    "fdr_threshold_applied": {
        "accession": "MS:4000227",
        "name": "PSM-level FDR threshold",
        "scope": "run",
        "value_type": "float"
    },
    "target_decoy_ratio": {
        "accession": "MS:4000228",
        "name": "target-decoy ratio of identifications",
        "scope": "run",
        "value_type": "float"
    },
    # --- New ProteoBench Specific Metrics ---
    "expected_vs_observed_ratios": {
        "accession": "MS:4000229",  # Note: Verify if custom or sub-terms are needed for ratios
        "name": "systematic ratio bias metric",
        "scope": "set",
        "value_type": "table"  # Stores expected vs. observed ratios per species mix (e.g. E. coli, Yeast, Human)
    },
    "coefficient_of_variation_distribution": {
        "accession": "MS:4000230",
        "name": "coefficient of variation distribution of protein abundances",
        "scope": "set",
        "value_type": "list"   # Excellent for plotting quantitative reproducibility across replicates
    },
    "species_specific_identification_counts": {
        "accession": "MS:4000231", 
        "name": "species-specific protein identification counts",
        "scope": "set",
        "value_type": "table"  # Tracks how many proteins were correctly assigned to human, yeast, or E. coli
    },

    ####SCOPE SET
    "pca_of_lfq_intensity": {
        "accession": "MS:4000091",
        "name": "principal component analysis of MaxQuant's protein group lfq intensities",
        "scope": "set",
    },
    "pca_of_raw_intensity": {
        "accession": "MS:4000092",
        "name": "identified MS1 feature area principal component analysis result",
        "scope": "set",
    },
    "potential_contaminants_per_group": {
        "accession": "MS:4000177",
        "name": "contaminant protein abundance fraction",
        "scope": "set",
        "value_type": "float",
    },
     "protein_quantification_matrix_summary": {
        "accession": "MS:4000223",
        "name": "protein level quantification matrix summary",
        "scope": "set",
        "value_type": "table"
    }
}