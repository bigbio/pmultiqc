from multiqc.plots import (
    bargraph,
    linegraph,
    box,
    scatter,
    table
)
from collections import OrderedDict

from ..core.section_groups import add_sub_section
from ..common.common_plots import remove_subtitle


# MaxQuant parameters table
def draw_parameters(sub_section, parameter_table):

    draw_config = {
        "namespace": "",
        "id": "parameters",
        "title": "Parameters",
        "save_file": False,
        "sort_rows": False,
        "only_defined_headers": True,
        "col1_header": "No.",
        "no_violin": True,
    }

    headers = {
        "parameter": {"title": "Parameter"},
        "value": {"title": "Value"}
    }

    table_html = table.plot(
        data=parameter_table,
        headers=headers,
        pconfig=draw_config
    )

    add_sub_section(
        sub_section=sub_section,
        plot=table_html,
        order=1,
        description="This table presents the parameters used in MaxQuant.",
        helptext="""
            MaxQuant parameters, extracted from parameters.txt, summarizes the settings used for the MaxQuant analysis. 
            Key parameters are MaxQuant version, Re-quantify, Match-between-runs and mass search tolerances. 
            A list of protein database files is also provided, allowing to track database completeness 
            and database version information (if given in the filename).
            """
    )

def draw_intensity_box(sub_section, distribution_box, fig_type):

    if distribution_box[1]:
        boxplot_label = ["Sample", "Contaminants"]
    else:
        boxplot_label = ["Sample"]
        distribution_box = distribution_box[:1]

    # 'intensity'
    if fig_type == "intensity":

        draw_config = {
            "id": "intensity_distribution_box",
            "cpswitch": False,
            "cpswitch_c_active": False,
            "title": "Intensity Distribution",
            "tt_decimals": 2,
            "data_labels": boxplot_label,
            "xlab": "log2(Intensity)",
        }

        box_html = box.plot(
            list_of_data_by_sample=distribution_box,
            pconfig=draw_config
        )

        box_html = remove_subtitle(box_html)

        add_sub_section(
            sub_section=sub_section,
            plot=box_html,
            order=3,
            description="",
            helptext="""
                Intensity boxplots by experimental groups. Groups are user-defined during MaxQuant configuration. 
                This plot displays a (customizable) threshold line for the desired mean intensity of proteins. 
                Groups which underperform here, are likely to also suffer from a worse MS/MS id rate and higher 
                contamination due to the lack of total protein loaded/detected. If possible, all groups should 
                show a high and consistent amount of total protein. 
                
                The height of the bar correlates to the number of proteins with non-zero abundance.
                """
        )

    # 'LFQ intensity'
    elif fig_type == "lfq_intensity":

        draw_config = {
            "id": "lfq_intensity_distribution_box",
            "cpswitch": False,
            "cpswitch_c_active": False,
            "title": "LFQ Intensity Distribution",
            "tt_decimals": 2,
            "data_labels": boxplot_label,
            "xlab": "log2(Intensity)",
        }

        box_html = box.plot(
            list_of_data_by_sample=distribution_box,
            pconfig=draw_config
        )

        box_html = remove_subtitle(box_html)

        add_sub_section(
            sub_section=sub_section,
            plot=box_html,
            order=4,
            description="Label-free quantification (LFQ) intensity boxplots by experimental groups.",
            helptext="""
                Label-free quantification (LFQ) intensity boxplots by experimental groups. Groups are user-defined 
                during MaxQuant configuration. This plot displays a (customizable) threshold line for the desired 
                mean of LFQ intensity of proteins. Raw files which underperform in *Raw* intensity, are likely to 
                show an *increased* mean here, since only high-abundance proteins are recovered and quantifyable by 
                MaxQuant in this Raw file. The remaining proteins are likely to receive an LFQ value of 0 (i.e. do 
                not contribute to the distribution).
                
                The height of the bar correlates to the number of proteins with non-zero abundance.
                """
        )

    # 'peptide intensity'
    elif fig_type == "peptide_intensity":

        draw_config = {
            "id": "peptide_intensity_distribution_box",
            "cpswitch": False,
            "cpswitch_c_active": False,
            "title": "Peptide Intensity Distribution",
            "tt_decimals": 2,
            "data_labels": boxplot_label,
            "xlab": "log2(Intensity)",
        }

        box_html = box.plot(
            list_of_data_by_sample=distribution_box,
            pconfig=draw_config
        )

        box_html = remove_subtitle(box_html)

        add_sub_section(
            sub_section=sub_section,
            plot=box_html,
            order=5,
            description="""
                Peptide precursor intensity per Raw file from evidence.txt WITHOUT match-between-runs evidence.
                """,
            helptext="""
                Peptide precursor intensity per Raw file from evidence.txt WITHOUT match-between-runs evidence. 
                Low peptide intensity usually goes hand in hand with low MS/MS identifcation rates and unfavourable 
                signal/noise ratios, which makes signal detection harder. Also instrument acquisition time increases 
                for trapping instruments. 
                Failing to reach the intensity threshold is usually due to unfavorable column conditions, inadequate 
                column loading or ionization issues. If the study is not a dilution series or pulsed SILAC experiment, 
                we would expect every condition to have about the same median log-intensity (of 2<sup>%1.1f</sup>). 
                The relative standard deviation (RSD) gives an indication about reproducibility across files and should 
                be below 5%%.
                """
        )

# MaxQuant Fig 6: PCA
def draw_pg_pca(sub_section, pca_data, fig_type):

    # fig_type: 'raw_intensity' or 'lfq_intensity'
    if fig_type == "raw_intensity":
        fig_id = "pca_of_raw_intensity"
        fig_title = "PCA of Raw Intensity"
        plot_order = 6

    if fig_type == "lfq_intensity":
        fig_id = "pca_of_lfq_intensity"
        fig_title = "PCA of LFQ Intensity"
        plot_order = 7

    draw_config = {
        "id": fig_id,
        "cpswitch": False,
        "cpswitch_c_active": False,
        "title": fig_title,
        "xlab": "PC #1",
        "ylab": "PC #2",
    }

    scatter_html = scatter.plot(
        data=pca_data,
        pconfig=draw_config
    )

    scatter_html = remove_subtitle(scatter_html)

    add_sub_section(
        sub_section=sub_section,
        plot=scatter_html,
        order=plot_order,
        description="""
            [Excludes Contaminants] Principal components plots of experimental groups (as defined 
            during MaxQuant configuration).
            """,
        helptext="""
            This plot is shown only if more than one experimental group was defined. 
            If LFQ was activated in MaxQuant, an additional PCA plot for LFQ intensities is shown. 
            Similarly, if iTRAQ/TMT reporter intensities are detected. 
            Since experimental groups and Raw files do not necessarily correspond 1:1, this plot 
            cannot use the abbreviated Raw file names, but instead must rely on automatic shortening 
            of group names.
            """
    )

# Peptide ID Count
def draw_evidence_peptide_id_count(sub_section, peptide_id_count_data):

    if peptide_id_count_data["title_value"]:
        fig_title = "Peptide ID Count" + " [" + peptide_id_count_data["title_value"] + "]"
    else:
        fig_title = "Peptide ID Count"

    draw_config = {
        "id": "peptide_id_count",
        "cpswitch": True,
        "title": fig_title,
        "tt_decimals": 0,
        "ylab": "Count",
    }

    bar_html = bargraph.plot(
        data=peptide_id_count_data["plot_data"],
        cats=peptide_id_count_data["cats"],
        pconfig=draw_config,
    )

    bar_html = remove_subtitle(bar_html)

    add_sub_section(
        sub_section=sub_section,
        plot=bar_html,
        order=4,
        description="""
            [Excludes Contaminants] Number of unique (i.e. not counted twice) peptide sequences including modifications (after FDR) per Raw file.
            """,
        helptext="""
            If MBR was enabled, three categories ('Genuine (Exclusive)', 'Genuine + Transferred', 'Transferred (Exclusive)'
            are shown, so the user can judge the gain that MBR provides. Peptides in the 'Genuine + Transferred' category 
            were identified within the Raw file by MS/MS, but at the same time also transferred to this Raw file using MBR. 
            This ID transfer can be correct (e.g. in case of different charge states), or incorrect -- see MBR-related 
            metrics to tell the difference. 
            Ideally, the 'Genuine + Transferred' category should be rather small, the other two should be large.

            If MBR would be switched off, you can expect to see the number of peptides corresponding to 'Genuine (Exclusive)' + 'Genuine + Transferred'. 
            In general, if the MBR gain is low and the MBR scores are bad (see the two MBR-related metrics),
            MBR should be switched off for the Raw files which are affected (could be a few or all). 
            """
    )

# ProteinGroups Count
def draw_evidence_protein_group_count(sub_section, protein_group_count_data):

    if protein_group_count_data["title_value"]:
        fig_title = (
            "ProteinGroups Count" + " [" + protein_group_count_data["title_value"] + "]"
        )
    else:
        fig_title = "ProteinGroups Count"

    draw_config = {
        "id": "protein_group_count",
        "cpswitch": True,
        "title": fig_title,
        "tt_decimals": 0,
        "ylab": "Count",
    }

    bar_html = bargraph.plot(
        data=protein_group_count_data["plot_data"],
        cats=protein_group_count_data["cats"],
        pconfig=draw_config,
    )

    bar_html = remove_subtitle(bar_html)

    add_sub_section(
        sub_section=sub_section,
        plot=bar_html,
        order=3,
        description="""
            [Excludes Contaminants] Number of Protein groups (after FDR) per Raw file.
            """,
        helptext="""
            If MBR was enabled, three categories ('Genuine (Exclusive)', 'Genuine + Transferred', 'Transferred (Exclusive)' 
            are shown, so the user can judge the gain that MBR provides. Here, 'Transferred (Exclusive)' means that this protein group 
            has peptide evidence which originates only from transferred peptide IDs. The quantification is (of course) always from the 
            local Raw file. 
            Proteins in the 'Genuine + Transferred' category have peptide evidence from within the Raw file by MS/MS, but at the same time 
            also peptide IDs transferred to this Raw file using MBR were used. It is not unusual to see the 'Genuine + Transferred' category be the 
            rather large, since a protein group usually has peptide evidence from both sources. 
            To see of MBR worked, it is better to look at the two MBR-related metrics.

            If MBR would be switched off, you can expect to see the number of protein groups corresponding to 'Genuine (Exclusive)' + 'Genuine + Transferred'. 
            In general, if the MBR gain is low and the MBR scores are bad (see the two MBR-related metrics), 
            MBR should be switched off for the Raw files which are affected (could be a few or all).
            """
    )

# Peak width over RT
def draw_evidence_peak_width_rt(sub_section, peak_rt_data):

    draw_config = {
        "id": "peak_width_over_RT",
        "cpswitch": False,
        "cpswitch_c_active": False,
        "title": "Peak width over RT",
        "ymin": 0,
        "tt_decimals": 3,
        "ylab": "Retention length [median]",
        "xlab": "Retention time [min]",
        "showlegend": True,
    }

    linegraph_html = linegraph.plot(
        data=peak_rt_data,
        pconfig=draw_config
    )

    linegraph_html = remove_subtitle(linegraph_html)

    add_sub_section(
        sub_section=sub_section,
        plot=linegraph_html,
        order=2,
        description="Distribution of widths of peptide elution peaks, derived from the evidence table.",
        helptext="""
            The distribution of the widths of peptide elution peaks, 
            derived from the evidence table and excluding potential contaminants, 
            is one parameter of optimal and reproducible chromatographic separation.
            """
    )

# Mass Error [ppm] boxplot
def draw_mass_error_box(sub_section, mass_error_data):

    max_abs_mass_error = max(
        abs(x) for values in mass_error_data.values() for x in values
    )

    if max_abs_mass_error <= 10:
        xmax_value = 10
        xmin_value = -10
    else:
        xmax_value = None
        xmin_value = None

    draw_config = {
        "id": "uncalibrated_mass_error_box",
        "cpswitch": False,
        "cpswitch_c_active": False,
        "title": "Uncalibrated Mass Error",
        "tt_decimals": 2,
        "xlab": "Mass Error [ppm]",
        "xmax": xmax_value,
        "xmin": xmin_value,
    }

    box_html = box.plot(
        list_of_data_by_sample=mass_error_data,
        pconfig=draw_config
    )

    box_html = remove_subtitle(box_html)

    add_sub_section(
        sub_section=sub_section,
        plot=box_html,
        order=5,
        description="[Excludes Contaminants] Mass accurary before calibration.",
        helptext="""
            Mass error of the uncalibrated mass-over-charge value of the precursor ion in comparison 
            to the predicted monoisotopic mass of the identified peptide sequence. 
            """
    )

# Summary Table
def draw_maxquant_summary_table(
        sub_section,
        msms_spectra,
        identified_msms,
        identified_peptides,
        protein_dict
    ):

    coverage = (identified_msms / msms_spectra) * 100

    summary_table = dict()
    summary_table[msms_spectra] = {
        "#Identified MS2 Spectra": identified_msms,
        "%Identified MS2 Spectra": coverage,
    }

    if identified_peptides:
        summary_table[msms_spectra]["#Peptides Identified"] = identified_peptides

    if protein_dict:
        if protein_dict["num_proteins_identified"]:
            summary_table[msms_spectra]["#Proteins Identified"] = protein_dict["num_proteins_identified"]
        if protein_dict["num_proteins_quantified"]:
            summary_table[msms_spectra]["#Proteins Quantified"] = protein_dict["num_proteins_quantified"]

    headers = OrderedDict()
    headers = {
        "#Identified MS2 Spectra": {
            "title": "#Identified MS2 Spectra",
            "description": "Total number of MS/MS spectra identified",
            "format": "{:,.0f}",
        },
        "%Identified MS2 Spectra": {
            "title": "%Identified MS2 Spectra",
            "description": "Percentage of Identified MS/MS Spectra",
            "suffix": "%",
            "format": "{:,.2f}",
        },
    }

    headers["#Identified MS2 Spectra"] = {
        "description": "Total number of MS/MS spectra identified",
        "format": "{:,.0f}",
    }
    headers["%Identified MS2 Spectra"] = {
        "description": "Percentage of Identified MS/MS Spectra",
        "format": "{:,.2f}",
        "suffix": "%",
    }

    # Create table plot
    pconfig = {
        "id": "identification_summary_table",
        "title": "Summary Table",
        "save_file": False,
        "raw_data_fn": "multiqc_summary_table_table",
        "sort_rows": False,
        "only_defined_headers": False,
        "col1_header": "#MS2 Spectra",
        "scale": "Set1",
    }

    table_html = table.plot(summary_table, headers, pconfig)

    add_sub_section(
        sub_section=sub_section,
        plot=table_html,
        order=1,
        description="This table shows the MaxQuant summary statistics.",
        helptext="""
            This table presents summary statistics generated by MaxQuant. 
            "#MS2 Spectra" is derived from msmsScans.txt (or msScans.txt); 
            "#Identified MS2 Spectra" and "#Peptides Identified" are derived from evidence.txt; 
            "#Proteins Identified" and "#Proteins Quantified" are derived from proteinGroups.txt.
            """
    )

# Number of Peptides identified Per Protein
def draw_maxquant_num_pep_pro(sub_section, num_pep_per_protein):

    data_labels = [
            {
                "name": "Frequency",
                "ylab": "Frequency",
                "tt_suffix": "",
                "tt_decimals": 0,
            },
            {
                "name": "Percentage",
                "ylab": "Percentage [%]",
                "tt_suffix": "%",
                "tt_decimals": 2,
            },
        ]

    pconfig = {
        "id": "number_of_peptides_per_proteins",
        "cpswitch": False,
        "title": "Number of Peptides identified Per Protein",
        "data_labels": data_labels,
    }

    bar_html = bargraph.plot(
        data=num_pep_per_protein,
        cats=["Frequency", "Percentage"],
        pconfig=pconfig
    )

    bar_html = remove_subtitle(bar_html)
    
    add_sub_section(
        sub_section=sub_section,
        plot=bar_html,
        order=1,
        description="This plot shows the number of peptides per protein in the MaxQuant data.",
        helptext="""
            This statistic is extracted from the proteinGroups.txt file. Proteins supported by more peptide 
            identifications can constitute more confident results.
            """
    )


# Search Engine Scores
def draw_maxquant_scores(sub_section, maxquant_scores):

    pconfig = {
        "id": "summary_of_andromeda_scores",
        "cpswitch": False,
        "title": "Summary of Andromeda Scores",
        "ylab": "Counts",
        "tt_suffix": "",
        "tt_decimals": 0,
        "data_labels": maxquant_scores["data_labels"],
    }

    bar_html = bargraph.plot(
        data=maxquant_scores["plot_data"],
        pconfig=pconfig
    )

    bar_html = remove_subtitle(bar_html)

    add_sub_section(
        sub_section=sub_section,
        plot=bar_html,
        order=1,
        description="",
        helptext="""
                This statistic is extracted from msms.txt. Andromeda score for the best associated MS/MS spectrum.
                """
    )

# TopN
def draw_msms_scans_top_n(sub_section, top_n_data):

    draw_config = {
        "id": "top_n",
        "cpswitch": True,
        "title": "TopN",
        "stacking": "group",
        "tt_decimals": 0,
        "ylab": "Highest Scan Event",
    }

    bar_html = bargraph.plot(
        data=top_n_data["plot_data"],
        cats=top_n_data["cats"],
        pconfig=draw_config
    )

    bar_html = remove_subtitle(bar_html)

    add_sub_section(
        sub_section=sub_section,
        plot=bar_html,
        order=9,
        description='This metric somewhat summarizes "TopN over RT"',
        helptext="""
            Reaching TopN on a regular basis indicates that all sections of the LC gradient 
            deliver a sufficient number of peptides to keep the instrument busy. This metric somewhat summarizes "TopN over RT".
            """
    )

# TopN over RT
def draw_msms_scans_top_over_rt(sub_section, top_over_rt_data):

    draw_config = {
        "id": "topn_over_rt",
        "cpswitch": False,
        "cpswitch_c_active": False,
        "title": "TopN over RT",
        "ymin": 0,
        "tt_decimals": 2,
        "ylab": "Highest N [median per RT bin]",
        "xlab": "Retention time [min]",
        "showlegend": True,
    }

    linegraph_html = linegraph.plot(
        data=top_over_rt_data,
        pconfig=draw_config
    )

    linegraph_html = remove_subtitle(linegraph_html)

    add_sub_section(
        sub_section=sub_section,
        plot=linegraph_html,
        order=8,
        description="TopN over retention time.",
        helptext="""
            TopN over retention time. Similar to ID over RT, this metric reflects the complexity of the sample 
            at any point in time. Ideally complexity should be made roughly equal (constant) by choosing a proper (non-linear) LC gradient. 
            See [Moruz 2014, DOI: 10.1002/pmic.201400036](https://pubmed.ncbi.nlm.nih.gov/24700534/) for details.
            """
    )

# Ion Injection Time over RT
def draw_msms_scans_ion_injec_time_rt(sub_section, ion_injec_time_rt_data):

    draw_config = {
        "id": "ion_injection_time_over_rt",
        "cpswitch": False,
        "cpswitch_c_active": False,
        "title": "Ion Injection Time over RT",
        "ymin": 0,
        "tt_decimals": 2,
        "ylab": "Ion injection time [ms]",
        "xlab": "Retention time [min]",
        "showlegend": True,
    }

    linegraph_html = linegraph.plot(
        data=ion_injec_time_rt_data,
        pconfig=draw_config
    )

    linegraph_html = remove_subtitle(linegraph_html)

    add_sub_section(
        sub_section=sub_section,
        plot=linegraph_html,
        order=7,
        description="",
        helptext="""
            Ion injection time score - should be as low as possible to allow fast cycles. Correlated with peptide intensity. 
            Note that this threshold needs customization depending on the instrument used (e.g., ITMS vs. FTMS).
            """
    )
