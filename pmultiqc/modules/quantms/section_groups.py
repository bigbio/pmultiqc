from multiqc.core.special_case_modules.custom_content import MultiqcModule, CcDict
from multiqc import report, config
from multiqc.types import Anchor, SectionId, Section
from multiqc.plots.plot import Plot


def create_ordered_group_modules(grouped_plot_data: list[dict]):

    for mod_cfg in grouped_plot_data:
        mod_id = mod_cfg["id"]
        mod_anchor = Anchor(mod_id)
        mod_config = mod_cfg.get("config", {})
        mod_title = mod_config.get("section_name", mod_id.replace("_", " ").title())

        group_module = MultiqcModule(
            id=mod_id,
            anchor=mod_anchor,
            cc_dict=CcDict(
                config={
                    "id": mod_id,
                    "section_name": mod_title,
                    "description": mod_config.get("description", f"Module for {mod_title}"),
                },
                data={}
            )
        )

        sorted_sections = sorted(mod_cfg.get("sections", []), key=lambda s: s.get("order") if s.get("order") is not None else 999)

        for i, sect in enumerate(sorted_sections):

            plot_obj = sect.get("plot")

            if isinstance(plot_obj, Plot):
                section_id = sect.get("id", f"{mod_id}_{plot_obj.pconfig.id}")
                section_name = sect.get("name", plot_obj.pconfig.title)
            else:
                section_id = sect.get("id", f"{mod_id}_section_{i}")
                section_name = sect.get("name", section_id.replace("_", " ").title())

            section = Section(
                name=section_name or "Unnamed Section",
                anchor=Anchor(section_id),
                id=SectionId(section_id),
                description=sect.get("description", ""),
                module=group_module.name,
                module_anchor=group_module.anchor,
                module_info=group_module.info,
                helptext=sect.get("helptext", ""),
                content_before_plot=sect.get("content_before_plot", ""),
                content=sect.get("content", ""),
                print_section=True
            )

            if isinstance(plot_obj, Plot):
                section.plot_anchor = plot_obj.anchor
                report.plot_by_id[plot_obj.anchor] = plot_obj
            elif isinstance(plot_obj, str):
                section.plot = plot_obj

            group_module.sections.append(section)

        report.modules.append(group_module)


def add_group_modules(groups_dict):

    section_group = list()

    group_configs = [
        {
            "id": "experiment_setup",
            "key": "experiment_sub_section",
            "name": "Experiment Setup",
            "description": ""
        },
        {
            "id": "summary_and_heatmap",
            "key": "summary_sub_section",
            "name": "Summary and HeatMap",
            "description": ""
        },
        {
            "id": "identification_summary",
            "key": "identification_sub_section",
            "name": "Identification Summary",
            "description": ""
        },
        {
            "id": "search_engine_scores",
            "key": "search_engine_sub_section",
            "name": "Search Engine Scores",
            "description": ""
        },
        {
            "id": "contaminants",
            "key": "contaminants_sub_section",
            "name": "Contaminants",
            "description": ""
        },
        {
            "id": "quantification_analysis",
            "key": "quantification_sub_section",
            "name": "Quantification Analysis",
            "description": ""
        },
        {
            "id": "ms1_analysis",
            "key": "ms1_sub_section",
            "name": "MS1 Analysis",
            "description": ""
        },
        {
            "id": "ms2_and_spectral_stats",
            "key": "ms2_sub_section",
            "name": "MS2 and Spectral Stats",
            "description": ""
        },
        {
            "id": "time_and_mass_error_trends",
            "key": "time_mass_sub_section",
            "name": "Time and Mass Error Trends",
            "description": ""
        },
    ]

    for group in group_configs:
        sections = groups_dict.get(group["key"])
        if sections:
            section_group.append({
                "id": group["id"],
                "config": {
                    "section_name": group["name"],
                    "description": group["description"]
                },
                "sections": sections,
            })

    create_ordered_group_modules(section_group)

    config.report_section_order = {
        "pmultiqc": {"order": 19},
        "experiment_setup": {"order": 18},
        "summary_and_heatmap": {"order": 17},
        "identification_summary": {"order": 16},
        "search_engine_scores": {"order": 15},
        "contaminants": {"order": 14},
        "quantification_analysis": {"order": 13},
        "ms1_analysis": {"order": 12},
        "ms2_and_spectral_stats": {"order": 11},
        "time_and_mass_error_trends": {"order": 10},
        "bigbio-quantms-summary": {"order": 9},
        "nf-core-quantms-methods-description": {"order": 8},
        "methods_description": {"order": 7},
        "bigbio-quantms-methods-description": {"order": 6},
        "software_versions": {"order": 5},
        "multiqc_software_versions": {"order": 4},
        "nf_core_quantms_software": {"order": 3},
        "nf-core-quantms-summary": {"order": 2},
        "workflow_summary": {"order": 1},
    }
