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




    ####Fallbacks
    "peptide_id_count": {
        "accession": "MS:4000001",
        "name": "number of peptide identifications",
        "scope": "run",
        "fallback": True
    },
    "msms_identified": {
        "accession": "MS:4000001",
        "name": "number of identified MS2 spectra",
        "scope": "run",
        "fallback": True
    },

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


}