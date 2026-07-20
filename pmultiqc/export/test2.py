
import json
###find exact key names & types
json_path = "/home/timo/Desktop/BA/data/mq_data/jsons/multiqc_data.json"

with open(json_path, "r") as f:
    data = json.load(f)

pep_plot = data.get("report_plot_data", {}).get("peptide_id_count", {})
datasets = pep_plot.get("datasets", [])

if datasets:
    ds = datasets[0]
    print("--- Checking All Keys inside Dataset for Data Arrays ---")
    
    # Check common Plotly array locations
    for key in ["x", "y", "series", "data", "cats"]:
        if key in ds:
            val = ds[key]
            print(f"Found '{key}' ({type(val)}): length = {len(val) if hasattr(val, '__len__') else 'N/A'}")
            print(f"  -> Preview: {str(val)[:200]}\n")
            
    # If none of those hit, let's look at everything that is a list or dict
    print("--- Full Dump of Remaining Keys with Types ---")
    for k, v in ds.items():
        if k not in ["x", "y", "series", "data", "cats", "samples", "layout", "trace_params"]:
            print(f"Key: {k} | Type: {type(v)} | Value Preview: {str(v)[:100]}")