import itertools
import numpy as np
import pandas as pd
import os
import re
from collections import OrderedDict
from sdrf_pipelines.openms.openms import UnimodDatabase
from multiqc.plots import table

from pmultiqc.modules.common.histogram import Histogram
from pmultiqc.modules.common.stats import qual_uniform
from pmultiqc.modules.common.plots import dia as dia_plots
from pmultiqc.modules.common.file_utils import file_prefix
from pmultiqc.modules.common.common_utils import evidence_rt_count, mod_group_percentage
from pmultiqc.modules.core.section_groups import add_sub_section
from pmultiqc.modules.common.plots.id import draw_ids_rt_count


from statsmodels.nonparametric.smoothers_lowess import lowess
from pmultiqc.modules.common.file_utils import drop_empty_row


DEFAULT_BINS = 500

def parse_diann_report(
        sub_sections,
        diann_report_path,
        heatmap_color_list,
        sample_df,
        file_df,
        ms_with_psm,
        cal_num_table_data,
        quantms_modified,
        ms_paths,
        msstats_input_valid=False
    ):

    log.info("Parsing {}...".format(diann_report_path))

    # parse DIA-NN report data
    if os.path.splitext(diann_report_path)[1] == ".tsv":
        report_data = pd.read_csv(
            diann_report_path, header=0, sep="\t", on_bad_lines="warn"
        )
    else:
        report_data = pd.read_parquet(diann_report_path)

    # "Decoy" appears only in DIA-NN 2.0 and later.
    # 0 or 1 based on whether the precursor is decoy, relevant when using --report-decoys
    if "Decoy" in report_data.columns:
        report_data = report_data[report_data["Decoy"] == 0].copy()

    # Normalisation.Factor: can be calculated as Precursor.Normalised/Precursor.Quantity
    required_cols = ["Precursor.Normalised", "Precursor.Quantity"]
    if "Normalisation.Factor" not in report_data.columns and all(
        col in report_data.columns for col in required_cols
    ):
        report_data["Normalisation.Factor"] = (
            report_data[required_cols[0]] / report_data[required_cols[1]]
        )

    # Draw: Standard Deviation of Intensity
    if "Precursor.Quantity" in report_data.columns:
        draw_dia_intensitys(sub_sections["quantification"], report_data)
        draw_dia_heatmap(sub_sections["summary"], report_data, heatmap_color_list)

    log.info("Draw the DIA MS1 subsection.")
    draw_dia_ms1(sub_sections["ms1"], report_data)

    log.info("Draw the DIA MS2 subsection.")
    draw_dia_ms2s(sub_sections["ms2"], report_data)

    log.info("Draw the DIA mass_error subsection.")
    draw_dia_mass_error(sub_sections["mass_error"], report_data)

    log.info("Draw the DIA rt_qc subsection.")
    draw_dia_rt_qc(sub_sections["rt_qc"], report_data)

    # Draw: Quantification Table (DIA-NN, without msstats data)
    if not msstats_input_valid:
        log.info("Draw the DIA quant table subsection.")
        draw_diann_quant_table(
            sub_sections["quantification"], report_data, sample_df, file_df
        )

    pattern = re.compile(r"\(.*?\)")
    report_data["sequence"] = [pattern.sub("", s) for s in report_data["Modified.Sequence"]]

    total_protein_quantified = len(set(report_data["Protein.Group"]))
    total_peptide_count = len(set(report_data["sequence"]))

    log.info("Processing DIA pep_plot.")
    protein_pep_map = report_data.groupby("Protein.Group")["sequence"].agg(list).to_dict()
    pep_plot = Histogram("number of peptides per proteins", plot_category="frequency")
    for _, peps in protein_pep_map.items():
        number = len(set(peps))
        pep_plot.add_value(number)

    categorys = OrderedDict()
    categorys["Frequency"] = {
        "name": "Frequency",
        "description": "number of peptides per proteins",
    }

    pep_plot.to_dict(percentage=True, cats=categorys)

    log.info("Processing DIA peptide_search_score.")
    peptide_search_score = dict()
    pattern = re.compile(r"\((.*?)\)")
    unimod_data = UnimodDatabase()
    for peptide, group in report_data.groupby("Modified.Sequence"):
        origianl_mods = re.findall(pattern, peptide)
        for mod in set(origianl_mods):
            name = unimod_data.get_by_accession(mod.upper()).get_name()
            peptide = peptide.replace(mod, name)
        if peptide.startswith("("):
            peptide = peptide + "."

        peptide_search_score[peptide] = np.min(group["Q.Value"])

    # Modifications Name
    log.info("Processing DIA Modifications.")
    mod_pattern = re.compile(r"\((.*?)\)")

    def find_diann_modified(peptide):
        if isinstance(peptide, str):
            mods = mod_pattern.findall(peptide)
            if mods:
                mod_type = [
                    unimod_data.get_by_accession(mod.upper()).get_name() for mod in set(mods)
                ]
                return ",".join(mod_type)
            else:
                return "Unmodified"
        return None

    report_data["Modifications"] = report_data["Modified.Sequence"].apply(
        lambda x: find_diann_modified(x)
    )

    log.info("Processing DIA mod_plot_dict.")
    mod_plot_dict = dict()
    modified_cats = list()
    for run_file, group in report_data.groupby("Run"):
        run_file = str(run_file)
        ms_with_psm.append(run_file)

        # Modifications
        mod_group_processed = mod_group_percentage(group.drop_duplicates())
        mod_plot_dict[run_file] = dict(
            zip(mod_group_processed["modifications"], mod_group_processed["percentage"])
        )
        modified_cats.extend(mod_group_processed["modifications"])

        cal_num_table_data[run_file] = {"protein_num": len(set(group["Protein.Group"]))}
        cal_num_table_data[run_file]["peptide_num"] = len(set(group["sequence"]))
        peptides = set(group["Modified.Sequence"])
        modified_pep = list(
            filter(lambda x: re.match(r".*?\(.*?\).*?", x) is not None, peptides)
        )
        group_peptides = group.groupby("sequence")["Protein.Group"].apply(list).to_dict()
        unique_peptides = [
            pep for pep, prots in group_peptides.items() if len(set(prots)) == 1
        ]
        cal_num_table_data[run_file]["unique_peptide_num"] = len(unique_peptides)
        cal_num_table_data[run_file]["modified_peptide_num"] = len(modified_pep)

    quantms_modified["plot_data"] = mod_plot_dict
    quantms_modified["cats"] = list(
        sorted(modified_cats, key=lambda x: (x == "Modified (Total)", x))
    )

    ms_without_psm = set([file_prefix(i) for i in ms_paths]) - set(ms_with_psm)
    for i in ms_without_psm:
        log.warning("No PSM found in '{}'!".format(i))

    for i in ms_without_psm:
        cal_num_table_data[i] = {
            "protein_num": 0,
            "peptide_num": 0,
            "unique_peptide_num": 0,
            "modified_peptide_num": 0,
        }

    return (
        total_protein_quantified,
        total_peptide_count,
        pep_plot,
        peptide_search_score,
        ms_with_psm,
        cal_num_table_data,
        quantms_modified,
        ms_without_psm
    )

def draw_dia_heatmap(sub_section, report_df, heatmap_color):

    log.info("Compute the Heatmap.")
    heatmap_data = cal_dia_heatmap(report_df)
    dia_plots.draw_heatmap(sub_section, heatmap_color, heatmap_data)
    log.info("Heatmap calculation is done.")


def draw_dia_intensitys(sub_section, report_df):

    df_sub = report_df[report_df["Precursor.Quantity"] > 0].copy()
    df_sub["log_intensity"] = np.log2(df_sub["Precursor.Quantity"])

    dia_plots.draw_dia_intensity_dis(sub_section, df_sub)

    if dia_plots.can_groupby_for_std(report_df, "Run"):
        dia_plots.draw_dia_intensity_std(sub_section, df_sub)


def draw_dia_ms1(sub_section, df):

    # Ms1.Area: non-normalised MS1 peak area
    if "Ms1.Area" in df.columns:
        df_sub = df[df["Ms1.Area"] > 0][["Ms1.Area", "Run"]].copy()
        if len(df_sub) > 0:
            df_sub["log_ms1_area"] = np.log2(df_sub["Ms1.Area"])
            dia_plots.draw_dia_ms1_area(sub_section, df_sub)


def draw_dia_ms2s(sub_section, df):

    # Distribution of Precursor Charges
    if "Precursor.Charge" in df.columns:
        dia_plots.draw_dia_whole_exp_charge(sub_section, df)

    # Charge-state of Per File
    dia_plots.draw_dia_ms2_charge(sub_section, df)


def draw_dia_mass_error(sub_section, df):

    # Ms1.Apex.Mz.Delta: difference between observed precursor m/z and the theoretical value
    if "Ms1.Apex.Mz.Delta" in df.columns:
        dia_plots.draw_dia_delta_mass(sub_section, df)


# DIA-NN: RT Quality Control
def draw_dia_rt_qc(sub_section, report_df):

    df = report_df.copy()

    # IDs over RT
    draw_dia_ids_rt(sub_section, df)

    # 1. Normalisation Factor over RT
    #   Normalisation.Factor: normalisation factor applied to the precursor in the specific run,
    #   i.e. normalised quantity = normalisation factor X non-normalised quantity
    if "Normalisation.Factor" in df.columns:
        log.info("Draw[rt_qc]: norm_factor_rt")
        norm_factor_rt = cal_feature_avg_rt(df, "Normalisation.Factor")
        dia_plots.draw_norm_factor_rt(sub_section, norm_factor_rt)

    # 2. FWHM over RT
    #   FWHM: estimated peak width at half-maximum
    if "FWHM" in df.columns:
        log.info("Draw[rt_qc]: draw_fwhm_rt")
        fwhm_rt = cal_feature_avg_rt(df, "FWHM")
        dia_plots.draw_fwhm_rt(sub_section, fwhm_rt)

    # 3. Peak Width over RT
    #   RT.Start and RT.Stop peak boundaries
    if all(col in df.columns for col in ["RT.Start", "RT.Stop"]):
        log.info("Draw[rt_qc]: draw_peak_width_rt")
        df["peak_width"] = df["RT.Stop"] - df["RT.Start"]
        peak_width_rt = cal_feature_avg_rt(df, "peak_width")
        dia_plots.draw_peak_width_rt(sub_section, peak_width_rt)

    # 4. Absolute RT Error over RT
    if all(col in df.columns for col in ["RT", "Predicted.RT"]):
        log.info("Draw[rt_qc]: draw_rt_error_rt")
        df["rt_error"] = abs(df["RT"] - df["Predicted.RT"])
        rt_error_rt = cal_feature_avg_rt(df, "rt_error")
        dia_plots.draw_rt_error_rt(sub_section, rt_error_rt)

    # 5. loess(RT ~ iRT)
    if all(col in df.columns for col in ["RT", "iRT"]):
        log.info("Draw[rt_qc]: draw_loess_rt_irt")
        rt_irt_loess = cal_rt_irt_loess(df)
        if rt_irt_loess is not None:
            dia_plots.draw_loess_rt_irt(sub_section, rt_irt_loess)


# DIA-NN: IDs over RT
def draw_dia_ids_rt(sub_section, report_df):

    rt_df = report_df[["Run", "RT"]].copy()
    rt_df.rename(columns={"Run": "raw file", "RT": "retention time"}, inplace=True)
    ids_over_rt = evidence_rt_count(rt_df)
    draw_ids_rt_count(sub_section, ids_over_rt, "dia")


# DIA-NN: Quantification Table
def draw_diann_quant_table(sub_section, diann_report, sample_df, file_df):

    # Peptides Quantification Table
    peptides_table, peptides_headers = create_peptides_table(
        diann_report, sample_df, file_df
    )
    draw_peptides_table(sub_section, peptides_table, peptides_headers, "DIA-NN")

    # Protein Quantification Table
    protein_table, protein_headers = create_protein_table(
        diann_report, sample_df, file_df
    )
    draw_protein_table(sub_section, protein_table, protein_headers, "DIA-NN")


# Draw: Peptides Quantification Table
def draw_peptides_table(sub_section, table_data, headers, report_type):

    draw_config = {
        "id": "peptides_quantification_table",
        "title": "Peptides Quantification Table",
        "save_file": False,
        "sort_rows": False,
        "only_defined_headers": True,
        "col1_header": "PeptideID",
        "no_violin": True,
    }

    # only use the first 50 lines for the table
    display_rows = 50
    table_html = table.plot(
        dict(itertools.islice(table_data.items(), display_rows)),
        headers=headers,
        pconfig=draw_config,
    )

    if report_type == "DIA-NN":
        description_text = """
            This plot shows the quantification information of peptides in the final result (DIA-NN report).
            """
        helptext_text = """
            The quantification information of peptides is obtained from the DIA-NN output file. 
            The table shows the quantitative level and distribution of peptides in different study variables, 
            run and peptiforms. The distribution show all the intensity values in a bar plot above and below 
            the average intensity for all the fractions, runs and peptiforms.

            * BestSearchScore: It is equal to min(1 - Q.Value) for DIA-NN datasets.
            * Average Intensity: Average intensity of each peptide sequence across all conditions (0 or NA ignored).
            * Peptide intensity in each condition (Eg. `CT=Mixture;CN=UPS1;QY=0.1fmol`).
            """
    elif report_type == "mzIdentML":
        description_text = """
            This plot shows the quantification information of peptides in the final result (mzIdentML).
            """
        helptext_text = """
            The quantification information of peptides is obtained from the mzIdentML. 
            The table shows the quantitative level and distribution of peptides in different study variables, 
            run and peptiforms. The distribution show all the intensity values in a bar plot above and below 
            the average intensity for all the fractions, runs and peptiforms.

            * BestSearchScore: It is equal to max(search_engine_score) for mzIdentML datasets.
            * Average Intensity: Average intensity of each peptide sequence (0 or NA ignored).
            """
    else:
        description_text = ""
        helptext_text = ""

    add_sub_section(
        sub_section=sub_section,
        plot=table_html,
        order=1,
        description=description_text,
        helptext=helptext_text,
    )


# Draw: Protein Quantification Table
def draw_protein_table(sub_sections, table_data, headers, report_type):

    draw_config = {
        "id": "protein_quant_result",
        "title": "Protein Quantification Table",
        "save_file": False,
        "sort_rows": False,
        "only_defined_headers": True,
        "col1_header": "ProteinID",
        "no_violin": True,
    }

    display_rows = 50
    table_html = table.plot(
        dict(itertools.islice(table_data.items(), display_rows)),
        headers=headers,
        pconfig=draw_config,
    )

    if report_type == "DIA-NN":
        description_text = """
            This plot shows the quantification information of proteins in the final result (DIA-NN report).
            """
        helptext_text = """
            The quantification information of proteins is obtained from the DIA-NN output file.
            The table shows the quantitative level and distribution of proteins in different study variables and run.

            * Peptides_Number: The number of peptides for each protein.
            * Average Intensity: Average intensity of each protein across all conditions (0 or NA ignored).
            * Protein intensity in each condition (Eg. `CT=Mixture;CN=UPS1;QY=0.1fmol`): Summarize intensity of peptides.
            """
    elif report_type == "mzIdentML":
        description_text = """
            This plot shows the quantification information of proteins in the final result (mzIdentML).
            """
        helptext_text = """
            The quantification information of proteins is obtained from the mzIdentML.
            The table shows the quantitative level and distribution of proteins in different study variables and run.

            * Peptides_Number: The number of peptides for each protein.
            * Average Intensity: Average intensity of each protein(0 or NA ignored).
            """
    else:
        description_text = ""
        helptext_text = ""

    add_sub_section(
        sub_section=sub_sections,
        plot=table_html,
        order=2,
        description=description_text,
        helptext=helptext_text,
    )

# DIA-NN: HeatMap
def cal_dia_heatmap(report_df):

    # "Contaminants" & "Peptide Intensity"
    pep_intensity = heatmap_cont_pep_intensity(report_df)

    # missed tryptic cleavages: there is no available data

    return pep_intensity


def heatmap_cont_pep_intensity(report_df):

    df = report_df[
        [
            "Run",
            "Protein.Names",
            "Precursor.Quantity",
            "RT",
            "Predicted.RT",
            "Precursor.Charge",
            "Normalisation.Factor",
            "RT.Stop",
            "RT.Start",
        ]
    ].copy()

    # TODO "CON"?
    df["is_contaminant"] = df["Protein.Names"].str.startswith("CON", na=False)

    # 3. "Charge"
    charge = dict()
    for raw_file, group in df[["Run", "Precursor.Charge"]].groupby("Run"):
        charge[raw_file] = group["Precursor.Charge"].value_counts()[2] / len(group)
    charge_median = np.median(list(charge.values()))
    heatmap_charge = dict(
        zip(
            charge.keys(),
            list(map(lambda v: float(1 - np.abs(v - charge_median)), charge.values())),
        )
    )

    heatmap_dict = {}
    for run, group in df.groupby("Run"):

        # 1. "Contaminants"
        cont_intensity_sum = group[group["is_contaminant"]]["Precursor.Quantity"].sum()
        if np.isnan(cont_intensity_sum) or cont_intensity_sum == 0:
            contaminant = 1
        else:
            intensity_sum = group["Precursor.Quantity"].sum()
            contaminant = 1 - cont_intensity_sum / intensity_sum

        # 2. "Peptide Intensity"
        pep_median = np.nanmedian(group["Precursor.Quantity"].to_numpy())
        pep_intensity = float(np.fmin(1.0, pep_median / (2**23)))

        # 4. "RT Alignment"
        rt_alignment = max(0.0, 1 - float(np.mean(np.abs(group["RT"] - group["Predicted.RT"]))))

        # 5. "ID rate over RT"
        ids_rate_over_rt = qual_uniform(group["RT"])

        # 6. Normalization Factor MAD
        def mean_abs_dev(x):
            mean_x = x.mean()
            return float(1 - np.mean(np.abs(x - mean_x)))

        norm_factor_mad = mean_abs_dev(group["Normalisation.Factor"])

        # 7. Peak Width = RT.Stop - RT.Start
        peak_width = max(0.0, 1 - float(np.mean(group["RT.Stop"] - group["RT.Start"])))

        # All Dict
        heatmap_dict[run] = {
            "Contaminants": contaminant,
            "Peptide Intensity": pep_intensity,
            "Charge": heatmap_charge.get(run, 0),
            "RT Alignment": rt_alignment,
            "ID rate over RT": ids_rate_over_rt,
            "Norm Factor": norm_factor_mad,
            "Peak Width": peak_width,
        }

    return heatmap_dict

def cal_feature_avg_rt(report_data, col):

    # RT: the retention time (RT) of the PSM in minutes
    sub_df = report_data[[col, "Run", "RT"]].copy()

    # RT bin
    sub_df["RT_bin"] = pd.cut(sub_df["RT"], bins=DEFAULT_BINS)
    sub_df["RT_bin_mid"] = sub_df["RT_bin"].apply(lambda x: x.mid)
    result = sub_df.groupby(["Run", "RT_bin_mid"], observed=False)[col].mean().reset_index()
    result[col] = result[col].fillna(0)

    plot_dict = {
        str(run): group.set_index("RT_bin_mid")[col].to_dict()
        for run, group in result.groupby("Run")
    }

    return plot_dict


# DIA-NN: Lowess (Loess)
def cal_rt_irt_loess(report_df, frac=0.3, data_bins: int = DEFAULT_BINS):

    if len(report_df) > 1000000:
        log.warning(f"Dataset too large ({len(report_df)} rows). Skipping LOWESS computation.")
        return None

    log.info("Start compute loess...")
    df = report_df.copy()

    # bin
    x_min, x_max = df["iRT"].min(), df["iRT"].max()
    bins = np.linspace(x_min, x_max, data_bins)

    plot_dict = dict()
    for run, group in df.groupby("Run"):

        group_sorted = group.sort_values("iRT")
        x = group_sorted["iRT"].values
        y = group_sorted["RT"].values

        # lowess
        smoothed = lowess(y, x, frac=frac)
        smoothed_x = smoothed[:, 0]
        smoothed_y = smoothed[:, 1]

        bin_indices = np.digitize(smoothed_x, bins)
        binned_dict = dict()
        for i in range(1, len(bins)):
            mask = bin_indices == i
            if np.any(mask):
                x_bin_mean = float(smoothed_x[mask].mean())
                y_bin_mean = float(smoothed_y[mask].mean())
                binned_dict[x_bin_mean] = y_bin_mean

        plot_dict[run] = binned_dict

    return plot_dict

# DIA-NN: Peptides Quantification Table
def create_peptides_table(report_df, sample_df, file_df):

    # Validation: remove rows with 0 or NA Precursor.Normalised values
    report_data = report_df[report_df["Precursor.Normalised"] > 0].copy()
    report_data = drop_empty_row(report_data, ["Protein.Names", "Stripped.Sequence"])

    report_data["BestSearchScore"] = 1 - report_data["Q.Value"]

    table_dict = dict()
    for sequence_protein, group in report_data.groupby(["Stripped.Sequence", "Protein.Names"]):

        table_dict[sequence_protein] = {
            "ProteinName": sequence_protein[1],
            "PeptideSequence": sequence_protein[0],
            "BestSearchScore": group["BestSearchScore"].min(),
            "Average Intensity": np.log10(group["Precursor.Normalised"].mean()),
        }

    headers = {
        "ProteinName": {
            "title": "Protein Name",
            "description": "Name/Identifier(s) of the protein (group)",
            "minrange": "200",
        },
        "PeptideSequence": {"title": "Peptide Sequence"},
        "BestSearchScore": {"title": "Best Search Score", "format": "{:,.4f}"},
        "Average Intensity": {
            "title": "Average Intensity",
            "description": "Average intensity across all conditions",
            "format": "{:,.4f}",
        },
    }

    if not sample_df.empty and not file_df.empty:

        sample_cond_df = pd.merge(
            sample_df[["Sample", "MSstats_Condition"]],
            file_df[["Sample", "Spectra_Filepath"]],
            on="Sample",
        )
        sample_cond_df["Run"] = sample_cond_df["Spectra_Filepath"].apply(
            lambda x: os.path.splitext(x)[0]
        )

        cond_report_data = pd.merge(
            report_data[["Stripped.Sequence", "Protein.Names", "Precursor.Normalised", "Run"]],
            sample_cond_df[["Run", "MSstats_Condition"]].drop_duplicates(),
            on="Run",
        )

        for sequence_protein, group in cond_report_data.groupby(
            ["Stripped.Sequence", "Protein.Names"]
        ):

            condition_data = dict()
            for condition, sub_group in group.groupby("MSstats_Condition"):
                condition_data[str(condition)] = np.log10(sub_group["Precursor.Normalised"].mean())

            table_dict[sequence_protein].update(condition_data)

        for exp_condition in sample_df["MSstats_Condition"].drop_duplicates():

            headers[str(exp_condition)] = {
                "title": str(exp_condition),
                "description": "MSstats Condition",
                "format": "{:,.4f}",
            }

    result_dict = {i: v for i, (_, v) in enumerate(table_dict.items(), start=1)}

    return result_dict, headers


# DIA-NN: Protein Quantification Table
def create_protein_table(report_df, sample_df, file_df):

    # Validation: remove rows with 0 or NA Precursor.Normalised values
    report_data = report_df[report_df["Precursor.Normalised"] > 0].copy()
    report_data = drop_empty_row(report_data, ["Protein.Names", "Stripped.Sequence"])

    table_dict = dict()
    for protein_name, group in report_data.groupby("Protein.Names"):

        table_dict[protein_name] = {
            "ProteinName": protein_name,
            "Peptides_Number": group["Stripped.Sequence"].nunique(),
            "Average Intensity": np.log10(group["Precursor.Normalised"].mean()),
        }

    headers = {
        "ProteinName": {
            "title": "Protein Name",
            "description": "Name/Identifier(s) of the protein (group)",
        },
        "Peptides_Number": {
            "title": "Number of Peptides",
            "description": "Number of peptides per proteins",
            "format": "{:,.0f}",
        },
        "Average Intensity": {
            "title": "Average Intensity",
            "description": "Average intensity across all conditions",
            "format": "{:,.4f}",
        },
    }

    if not sample_df.empty and not file_df.empty:

        sample_cond_df = pd.merge(
            sample_df[["Sample", "MSstats_Condition"]],
            file_df[["Sample", "Spectra_Filepath"]],
            on="Sample",
        )
        sample_cond_df["Run"] = sample_cond_df["Spectra_Filepath"].apply(
            lambda x: os.path.splitext(x)[0]
        )

        cond_report_data = pd.merge(
            report_data[["Stripped.Sequence", "Protein.Names", "Precursor.Normalised", "Run"]],
            sample_cond_df[["Run", "MSstats_Condition"]].drop_duplicates(),
            on="Run",
        )

        for protein_name, group in cond_report_data.groupby("Protein.Names"):

            condition_data = dict()
            for condition, sub_group in group.groupby("MSstats_Condition"):
                condition_data[str(condition)] = np.log10(sub_group["Precursor.Normalised"].mean())

            table_dict[protein_name].update(condition_data)

        for exp_condition in sample_df["MSstats_Condition"].drop_duplicates():

            headers[str(exp_condition)] = {
                "title": str(exp_condition),
                "description": "MSstats Condition",
                "format": "{:,.4f}",
            }

    result_dict = {i: v for i, (_, v) in enumerate(table_dict.items(), start=1)}

    return result_dict, headers

