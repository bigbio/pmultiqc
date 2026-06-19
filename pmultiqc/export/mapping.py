###Mapping metrics to https://github.com/HUPO-PSI/psi-ms-CV/blob/master/psi-ms.obo

MZQC_METRIC_MAPPING = {
   ####SCOPE RUN
    "missed_cleavages": {
        "accession": "MS:4000215", 
        "name": "missed cleavages fractions",
        "scope": "run",
        "fallback": False
    },
    "protein_group_count": {
        "accession": "MS:1003327",
        "name": "number of identified protein groups",
        "scope": "run",
        "fallback": False
    },
    "charge_state": {
        "accession": "MS:4000208",
        "name": "identified quantification data points charges fractions",
        "scope": "run",
        "fallback": False
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
        "value_type": "float"  #  enforces xsd:double/float
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
        "value_type": "list"  # Expects a list of 3 floats: [Q1, Q2, Q3]
    },
    ####---
    ####---
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
        "value_type": "float" # Single ratio value
    },####---

    ####SCOPE SET
    "pca_of_lfq_intensity": {
        "accession": "MS:4000091",
        "name": "principal component analysis of MaxQuant's protein group lfq intensities",
        "scope": "set",
        "fallback": False
    },
    "pca_of_raw_intensity": {
        "accession": "MS:4000092",
        "name": "identified MS1 feature area principal component analysis result",
        "scope": "set",
        "fallback": False
    },
    "potential_contaminants_per_group": {
        "accession": "MS:4000177",
        "name": "contaminant protein abundance fraction",
        "scope": "set",
        "value_type": "float",
        "fallback": False
    },


}