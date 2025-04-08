""" MultiQC functions to plot a table """

from collections import defaultdict, OrderedDict
import logging
import textwrap
from multiqc import config, report
from multiqc.utils import mqc_colour
from multiqc.plots import table_object, violin
from multiqc.plots.table_object import TableConfig
from multiqc.plots.table_object import ValueT
from typing import Dict

logger = logging.getLogger(__name__)

letters = "abcdefghijklmnopqrstuvwxyz"


def plot(data, headers=None, pconfig=None, max_value=0.0):

    if headers is None:
        headers = []
    if pconfig is None:
        pconfig = {}

    plot_anchor = pconfig["anchor"]
    del pconfig["anchor"]
    table_pconfig = TableConfig(**pconfig)

    # Make a DataTable object
    dt = table_object.DataTable.create(
        data=data,
        table_id=pconfig["id"],
        table_anchor=plot_anchor,
        pconfig=table_pconfig,
        headers=headers,
    )

    # s_names = set()
    # for d in dt.raw_data:
    #     for s_name in d.keys():
    #         s_names.add(s_name)
    s_names = set(dt.sections[0].rows_by_sgroup.keys())

    # Make a violin plot if we have lots of samples
    if len(s_names) >= config.max_table_rows and pconfig.no_violin is not True:
        logger.debug("Plotting violin instead of table, {} samples".format(len(s_names)))
        warning = (
            '<p class="text-muted"><span class="glyphicon glyphicon-exclamation-sign" '
            'title="A violin plot has been generated instead because of the large number of samples. '
            'See http://multiqc.info/docs/#tables--beeswarm-plots"'
            ' data-toggle="tooltip"></span> Showing {} samples.</p>'.format(len(s_names))
        )
        return warning + violin.plot(data, headers, pconfig)
    else:
        return make_table(dt, max_value)


def make_table(dt, max_value):
    # table_id = dt.pconfig.get("id", "table_{}".format("".join(random.sample(letters, 4))))
    table_id = dt.id
    table_id = report.save_htmlid(table_id)
    t_headers = OrderedDict()
    t_modal_headers = OrderedDict()
    t_rows = OrderedDict()
    t_rows_empty = OrderedDict()
    raw_vals: Dict[str, Dict[str, ValueT]] = defaultdict(lambda: dict())
    empty_cells = dict()
    hidden_cols = 1
    table_title = dt.pconfig.title
    if table_title is None:
        table_title = table_id.replace("_", " ").title()

    fixed_col = [
        "PeptideSequence",
        "ProteinName",
        "BestSearchScore",
        "Average Intensity",
        "Peptides_Number",
        "Average Spectrum Counting",
    ]
    # for idx, k, header in dt.get_headers_in_order():
    for _, k, header in dt.get_headers_in_order():

        rid = header.rid
        # Build the table header cell
        shared_key = ""
        if header.shared_key is not None:
            shared_key = " data-shared-key={}".format(header.shared_key)

        hide = ""
        muted = ""
        checked = ' checked="checked"'
        if header.hidden:
            hide = "hidden"
            muted = " text-muted"
            checked = ""
            hidden_cols += 1

        data_attr = 'data-dmax="{}" data-dmin="{}" data-namespace="{}" {}'.format(
            header.dmax, header.dmin, header.namespace, shared_key
        )

        # join with zero width white space for line break
        cont = "&#8203;".join(textwrap.wrap(header.title, 30))

        cell_contents = '<span class="mqc_table_tooltip" title="{}: {}">{}</span>'.format(
            header.namespace, header.description, cont
        )

        if k not in fixed_col and "distribution" in k:
            cont = header.title.replace("_distribution", "")
            cont = "&#8203;".join(textwrap.wrap(cont, 30))
            cell_contents = '<span class="mqc_table_tooltip" title="{}: {}">{}</span>'.format(
                header.namespace, header.description, cont
            )
            t_headers[rid] = (
                '<th id="header_{rid}" class="{rid} {h} col-condition-sparkline" {da}>{c}</th>'.format(
                    rid=rid, h=hide, da=data_attr, c=cell_contents
                )
            )
        elif k not in fixed_col:
            t_headers[rid] = (
                '<th id="header_{rid}" class="{rid} {h} col-condition" {da}>{c}</th>'.format(
                    rid=rid, h=hide, da=data_attr, c=cell_contents
                )
            )
        else:
            t_headers[rid] = '<th id="header_{rid}" class="{rid} {h}" {da}>{c}</th>'.format(
                rid=rid, h=hide, da=data_attr, c=cell_contents, t=table_id
            )
        empty_cells[rid] = '<td class="data-coloured {rid} {h}"></td>'.format(rid=rid, h=hide)

        # Build the modal table row
        t_modal_headers[rid] = (
            """
        <tr class="{rid}{muted}" style="background-color: rgba({col}, 0.15);">
            <td class="sorthandle ui-sortable-handle">||</span></td>
            <td style="text-align:center;">
                <input class="mqc_table_col_visible" type="checkbox" {checked} value="{rid}" data-target="#{tid}">
            </td>
            <td>{name}</td>
            <td>{title}</td>
            <td>{desc}</td>
            <td>{col_id}</td>
            <td>{sk}</td>
        </tr>""".format(
                rid=rid,
                muted=muted,
                checked=checked,
                tid=table_id,
                col=header.colour,
                name=header.namespace,
                title=header.title,
                desc=header.description,
                col_id="<code>{}</code>".format(k),
                sk=header.shared_key,
            )
        )

        # Make a colour scale
        if not header.scale:
            c_scale = None
        else:
            c_scale = mqc_colour.mqc_colour_scale(header.scale, header.dmin, header.dmax, id=dt.id)

        # Collect conditional formatting config
        cond_formatting_rules = {}
        if header.cond_formatting_rules:
            cond_formatting_rules[rid] = header.cond_formatting_rules
        cond_formatting_rules.update(config.table_cond_formatting_rules)

        cond_formatting_colours = header.cond_formatting_colours
        cond_formatting_colours.extend(config.table_cond_formatting_colours)

        # Add the data table cells
        # for s_name, samp in dt.raw_data[idx].items():
        # Multiqc 1.22 --> 1.26
        for _, raw_data in dt.sections[0].rows_by_sgroup.items():
            s_name = raw_data[0].sample
            samp = raw_data[0].raw_data

            if k in samp:
                val: ValueT = samp[k]
                # valstr: str = dt.formatted_data[idx][s_name][k]
                valstr: str = str(samp[k])
                # kname = "{}_{}".format(header["namespace"], rid)
                # dt.raw_vals[s_name][kname] = val
                raw_vals[s_name][f"{header.namespace}_{rid}"] = val

                if c_scale and c_scale.name not in c_scale.qualitative_scales:
                    dmin = header.dmin
                    dmax = header.dmax
                    if dmin is not None and dmax is not None and dmax != dmin:
                        try:
                            val_float = float(val)
                        except ValueError:
                            percentage = 0.0
                        else:
                            percentage = ((val_float - dmin) / (dmax - dmin)) * 100
                            # Treat 0 as 0-width and make bars width of absolute value
                            if header.bars_zero_centrepoint:
                                dmax = max(abs(dmin), abs(dmax))
                                dmin = 0
                                percentage = ((abs(val_float) - dmin) / (dmax - dmin)) * 100
                            percentage = min(percentage, 100)
                            percentage = max(percentage, 0)
                    else:
                        percentage = 0.0
                else:
                    percentage = 100.0

                # This is horrible, but Python locale settings are worse
                if config.thousandsSep_format is None:
                    config.thousandsSep_format = '<span class="mqc_small_space"></span>'
                if config.decimalPoint_format is None:
                    config.decimalPoint_format = "."

                valstr = valstr.replace(".", "DECIMAL").replace(",", "THOUSAND")
                valstr = valstr.replace("DECIMAL", config.decimalPoint_format).replace(
                    "THOUSAND", config.thousandsSep_format
                )

                suffix = header.suffix
                if suffix:
                    # Add a space before the suffix, but not as an actual character, so ClipboardJS would copy
                    # the whole value without the space. Also, remove &nbsp; that we don't want ClipboardJS to copy.
                    suffix = suffix.replace("&nbsp;", " ").strip()
                    valstr += "<span class='mqc_small_space'></span>" + suffix

                # Conditional formatting
                # Build empty dict for cformatting matches
                cmatches = {}
                for cfc in cond_formatting_colours:
                    for cfck in cfc:
                        cmatches[cfck] = False
                # Find general rules followed by column-specific rules
                for cfk in ["all_columns", rid, table_id]:
                    if cfk in cond_formatting_rules:
                        # Loop through match types
                        for ftype in cmatches.keys():
                            # Loop through array of comparison types
                            for cmp in cond_formatting_rules[cfk].get(ftype, []):
                                try:
                                    # Each comparison should be a dict with single key: val
                                    if (
                                        "s_eq" in cmp
                                        and str(cmp["s_eq"]).lower() == str(val).lower()
                                    ):
                                        cmatches[ftype] = True
                                    if (
                                        "s_contains" in cmp
                                        and str(cmp["s_contains"]).lower() in str(val).lower()
                                    ):
                                        cmatches[ftype] = True
                                    if (
                                        "s_ne" in cmp
                                        and str(cmp["s_ne"]).lower() != str(val).lower()
                                    ):
                                        cmatches[ftype] = True
                                    if "eq" in cmp and float(cmp["eq"]) == float(val):
                                        cmatches[ftype] = True
                                    if "ne" in cmp and float(cmp["ne"]) != float(val):
                                        cmatches[ftype] = True
                                    if "gt" in cmp and float(cmp["gt"]) < float(val):
                                        cmatches[ftype] = True
                                    if "lt" in cmp and float(cmp["lt"]) > float(val):
                                        cmatches[ftype] = True
                                except:
                                    logger.warning(
                                        "Not able to apply table conditional formatting to '{}' ({})".format(
                                            val, cmp
                                        )
                                    )
                # Apply HTML in order of config keys
                badge_col = None
                for cfc in cond_formatting_colours:
                    for cfck in cfc:  # should always be one, but you never know
                        if cmatches[cfck]:
                            badge_col = cfc[cfck]
                if badge_col is not None:
                    valstring = '<span class="badge" style="background-color:{}">{}</span>'.format(
                        badge_col, valstr
                    )

                # Determine background color based on scale. Only relevant for hashable values. If value is for some
                # reason a dict or a list, it's not hashable and the logic determining the color will not work.
                hashable = True
                try:
                    hash(val)
                except TypeError:
                    hashable = False
                    print(
                        f"Value {val} is not hashable for table {dt.id}, column {k}, sample {s_name}"
                    )

                # Categorical backgorund colours supplied
                if val in header.bgcols.keys():
                    col = 'style="background-color:{} !important;"'.format(header.bgcols[val])
                    if s_name not in t_rows:
                        t_rows[s_name] = dict()
                    t_rows[s_name][rid] = '<td class="{rid} {h}" {c}>{v}</td>'.format(
                        rid=rid, h=hide, c=col, v=valstr
                    )

                # Build table cell background colour bar
                elif hashable and header.scale:
                    if c_scale is not None:
                        col = " background-color:{} !important;".format(c_scale.get_colour(val))
                    else:
                        col = ""
                    bar_html = '<span class="bar" style="width:{}%;{}"></span>'.format(
                        percentage, col
                    )
                    val_html = '<span class="val">{}</span>'.format(valstr)
                    wrapper_html = '<div class="wrapper">{}</div>'.format(val_html)

                    if s_name not in t_rows:
                        t_rows[s_name] = dict()
                    if "_distribution" in rid:
                        if valstr == "":
                            t_rows[s_name][rid] = (
                                "<td class=\"data-sparkline col-condition-sparkline\" data-sparkline='{v}'></td>".format(
                                    rid=rid, h=hide, v=valstr
                                )
                            )
                        else:
                            valstr.replace(",", "&#44;")
                            # valstring = ", ".join(valstring.split(" ;")[0].split(" ")) + " ;" + str(valstring.split(" ;")[1])
                            t_rows[s_name][rid] = (
                                "<td class=\"data-sparkline col-condition-sparkline\" data-sparkline='{v}'></td>".format(
                                    rid=rid, h=hide, v=valstr
                                )
                            )
                    elif header.title in fixed_col:
                        if "Average_Intensity-1" in rid:
                            t_rows[s_name][rid] = (
                                '<td class="data-coloured Average_Intensity {h}">{c}</td>'.format(
                                    h=hide, c=wrapper_html
                                )
                            )
                        else:
                            t_rows[s_name][rid] = (
                                '<td class="data-coloured {rid} {h}">{c}</td>'.format(
                                    rid=rid, h=hide, c=wrapper_html
                                )
                            )
                    else:
                        t_rows[s_name][rid] = (
                            '<td class="data-coloured col-condition {rid} {h}">{c}</td>'.format(
                                rid=rid, h=hide, c=wrapper_html
                            )
                        )

                # Scale / background colours are disabled
                else:
                    if s_name not in t_rows:
                        t_rows[s_name] = dict()
                    t_rows[s_name][rid] = '<td class="{rid} {h}">{v}</td>'.format(
                        rid=rid, h=hide, v=valstr
                    )

                # Is this cell hidden or empty?
                if s_name not in t_rows_empty:
                    t_rows_empty[s_name] = dict()
                t_rows_empty[s_name][rid] = header.hidden or str(val).strip() == ""

        # Remove header if we don't have any filled cells for it
        if sum([len(rows) for rows in t_rows.values()]) == 0:
            if header.get("hidden", False) is True:
                hidden_cols -= 1
            t_headers.pop(rid, None)
            t_modal_headers.pop(rid, None)
            logger.debug("Removing header {} from table, as no data".format(k))

    #
    # Put everything together
    #

    # Buttons above the table
    html = ""
    if not config.simple_output:

        # Copy Table Button
        html += """
        <button type="button" class="mqc_table_copy_btn btn btn-default btn-sm" data-clipboard-target="#{tid}">
            <span class="glyphicon glyphicon-copy"></span> Copy table
        </button>
        """.format(
            tid=table_id
        )

        # Configure Columns Button
        if len(t_headers) > 1:
            html += """
            <button type="button" class="mqc_table_configModal_btn btn btn-default btn-sm" data-toggle="modal" data-target="#{tid}_configModal">
                <span class="glyphicon glyphicon-th"></span> Configure Columns
            </button>
            """.format(
                tid=table_id
            )

        # Sort By Highlight button
        html += """
        <button type="button" class="mqc_table_sortHighlight btn btn-default btn-sm" data-target="#{tid}" data-direction="desc" style="display:none;">
            <span class="glyphicon glyphicon-sort-by-attributes-alt"></span> Sort by highlight
        </button>
        """.format(
            tid=table_id
        )

        # Scatter Plot Button
        if len(t_headers) > 1:
            html += """
            <button type="button" class="mqc_table_makeScatter btn btn-default btn-sm" data-toggle="modal" data-target="#tableScatterModal" data-table="#{tid}">
                <span class="glyphicon glyphicon glyphicon-stats"></span> Plot
            </button>
            """.format(
                tid=table_id
            )

        # "Showing x of y columns" text
        row_visibilities = [all(t_rows_empty[s_name].values()) for s_name in t_rows_empty]
        visible_rows = [x for x in row_visibilities if not x]

        # Visible rows
        t_showing_rows_txt = 'Showing <sup id="{tid}_numrows" class="mqc_table_numrows">{nvisrows}</sup>/<sub>{nrows}</sub> rows'.format(
            tid=table_id, nvisrows=len(visible_rows), nrows=len(t_rows)
        )

        # How many columns are visible?
        ncols_vis = (len(t_headers) + 1) - hidden_cols
        t_showing_cols_txt = ""
        if len(t_headers) > 1:
            t_showing_cols_txt = ' and <sup id="{tid}_numcols" class="mqc_table_numcols">{ncols_vis}</sup>/<sub>{ncols}</sub> columns'.format(
                tid=table_id, ncols_vis=ncols_vis, ncols=len(t_headers)
            )

        # Build table header text
        html += """
        <small id="{tid}_numrows_text" class="mqc_table_numrows_text">{rows}{cols}.</small>
        """.format(
            tid=table_id, rows=t_showing_rows_txt, cols=t_showing_cols_txt
        )
        if "peptide" in table_title:
            html += """
            <button type="button" class="btn btn-default btn-sm" new-data=1 data-action="set_data" data-target="quantification_of_peptides" id="peptide-distribution-button">
            <span class="glyphicon glyphicon glyphicon-stats"></span> Show replicates</button>
            """
        else:
            html += """
            <button type="button" class="btn btn-default btn-sm" new-data=1 data-action="set_data" data-target="quantification_of_proteins" id="protein-distribution-button">
            <span class="glyphicon glyphicon glyphicon-stats"></span> Show replicates</button>
            """

    # Build the table itself
    collapse_class = "mqc-table-collapse" if len(t_rows) > 10 and config.collapse_tables else ""
    html += """
        <div id="{tid}_container" class="mqc_table_container">
            <div class="table-responsive mqc-table-responsive {cc}">
                <table id="{tid}" class="table table-condensed mqc_table" data-title="{title}" data-dmax="{dmax}" >
        """.format(
        tid=table_id, title=table_title, cc=collapse_class, dmax=max_value
    )

    # Build the header row
    col1_header = dt.pconfig.col1_header

    html += '<thead><tr><th class="rowheader">{}</th>{}</tr></thead>'.format(
        col1_header, "".join(t_headers.values())
    )

    # Build the table body
    html += "<tbody>"
    t_row_keys = t_rows.keys()
    if dt.pconfig.sort_rows:
        t_row_keys = sorted(t_row_keys)
    for s_name in t_row_keys:
        # Hide the row if all cells are empty or hidden
        row_hidden = ' style="display:none"' if all(t_rows_empty[s_name].values()) else ""
        html += "<tr{}>".format(row_hidden)
        # Sample name row header
        # Wrap with zero width space character for line breaks that is not visible later
        content = "&#8203;".join(textwrap.wrap(s_name, 40))
        html += '<th class="rowheader" data-original-sn="{sn}">{sn}</th>'.format(sn=content)
        for k in t_headers:
            html += t_rows[s_name].get(k, empty_cells[k])
        html += "</tr>"
    html += "</tbody></table></div>"
    if len(t_rows) > 10 and config.collapse_tables:
        html += '<div class="mqc-table-expand"><span class="glyphicon glyphicon-chevron-down" aria-hidden="true"></span></div>'
    html += "</div>"

    # Build the bootstrap modal to customise columns and order
    if not config.simple_output:
        html += """
    <!-- MultiQC Table Columns Modal -->
    <div class="modal fade" id="{tid}_configModal" tabindex="-1">
        <div class="modal-dialog modal-lg">
            <div class="modal-content">
            <div class="modal-header">
                <button type="button" class="close" data-dismiss="modal" aria-label="Close"><span aria-hidden="true">&times;</span></button>
                <h4 class="modal-title">{title}: Columns</h4>f
            </div>
            <div class="modal-body">
                <p>Uncheck the tick box to hide columns. Click and drag the handle on the left to change order.</p>
                <p>
                    <button class="btn btn-default btn-sm mqc_configModal_bulkVisible" data-target="#{tid}" data-action="showAll">Show All</button>
                    <button class="btn btn-default btn-sm mqc_configModal_bulkVisible" data-target="#{tid}" data-action="showNone">Show None</button>
                </p>
                <table class="table mqc_table mqc_sortable mqc_configModal_table" id="{tid}_configModal_table" data-title="{title}">
                <thead>
                    <tr>
                    <th class="sorthandle" style="text-align:center;">Sort</th>
                    <th style="text-align:center;">Visible</th>
                    <th>Group</th>
                    <th>Column</th>
                    <th>Description</th>
                    <th>ID</th>
                    <th>Scale</th>
                    </tr>
                </thead>
                <tbody>
                    {trows}
                </tbody>
            </table>
        </div>
        <div class="modal-footer"> <button type="button" class="btn btn-default" data-dismiss="modal">Close</button> </div>
    </div> </div> </div>""".format(
            tid=table_id, title=table_title, trows="".join(t_modal_headers.values())
        )

    # Save the raw values to a file if requested
    if dt.pconfig.save_file:
        fn = dt.pconfig.get("raw_data_fn", "multiqc_{}".format(table_id))
        report.write_data_file(raw_vals, fn)
        report.saved_raw_data[fn] = dt.raw_vals

    return html
