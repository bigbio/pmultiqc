import os
import pandas as pd

from pmultiqc.modules.common.logging import get_logger


# Initialise the module logger via centralized logger
log = get_logger("pmultiqc.modules.fragpipe.fragpipe_io")


REQUIRED_COLS = {
    "psm": [
        "Spectrum", "Peptide", "Modified Peptide", "Charge", "Retention", "Intensity", 
        "Delta Mass", "Number of Missed Cleavages", "Is Unique", "Protein"
    ],
    "peptide": []
}


# FragPipe File Paths
def get_fragpipe_files(find_log_files):

    # FragPipe results (https://fragpipe.nesvilab.org/docs/tutorial_fragpipe_outputs.html)
    # Main report files
    # psm.tsv (from Philosopher, updated by PTM-Shepherd and IonQuant)
    # ion.tsv (from Philosopher, overwritten by IonQuant)
    # peptide.tsv (from Philosopher, overwritten by IonQuant)
    # protein.tsv (from Philosopher, overwritten by IonQuant)
    # combined_ion.tsv (from Philosopher, overwritten by IonQuant)
    # combined_modified_peptide.tsv (from IonQuant)
    # combined_peptide.tsv (from Philosopher, overwritten by IonQuant)
    # combined_protein.tsv (from Philosopher, overwritten by IonQuant)
    # diann-output files (see DIA-NN documentation)

    required_files = ["psm"]
    req_set = set(required_files)

    fragpipe_files = {req: [] for req in required_files}

    # FragPipe *tsv Data
    for file_info in find_log_files("pmultiqc/tsv", filecontents=False):
        filename = file_info["fn"]
        full_path = os.path.join(file_info["root"], filename)

        for req in req_set:
            if req in filename:

                if req == "protein" and "combined_protein" in filename:
                    continue

                if req == "peptide" and (
                    "combined_peptide" in filename or "combined_modified_peptide" in filename
                ):
                    continue

                fragpipe_files[req].append(full_path)

    if fragpipe_files:

        for k, v in fragpipe_files.items():
            log.info(f"FragPipe data loaded: {k} ({len(v)} files).")
            log.debug(f"FragPipe data loaded: {k}: {v}")

        return fragpipe_files

    else:
        log.warning("FragPipe data loading failed: Please verify the input folder path.")

        return None


def psm_reader(file_path: str):

    psm_df = pd.read_csv(file_path, sep="\t")
    psm_df["Run"] = psm_df["Spectrum"].astype(str).apply(
        lambda x: x.rsplit(".", 3)[0]
    )

    if check_columns(file_path=file_path, data_type="psm"):

        required_cols = REQUIRED_COLS["psm"].copy()
        required_cols.append("Run")

        return psm_df[required_cols].copy()

    return psm_df


def check_columns(file_path: str, data_type: str):

    try:
        df_header = pd.read_csv(file_path, sep="\t", nrows=0)
        actual_columns = set(df_header.columns.str.strip())

        missing_columns = [
            col
            for col in REQUIRED_COLS[data_type]
            if col not in actual_columns
        ]

        if not missing_columns:
            log.info("Check passed: All required columns are present.")
            return True
        else:
            log.info(
                f"Check failed: {file_path}. Missing the following columns: {missing_columns}"
            )
            return False

    except Exception as e:
        log.warning(f"Check failed: Unable to read file: {e}")
        return False

