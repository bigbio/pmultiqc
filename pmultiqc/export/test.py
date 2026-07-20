# test_batch_datasets.py
import os
import logging

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

from pmultiqc.export.mzqc_writer import write_mzqc

def run_batch_test():
    # Targeted folder containing only your json datasets
    base_dir = "/home/timo/Desktop/BA/data/mq_data/jsons"
    target_filename = "multiqc_data.json"
    
    if not os.path.exists(base_dir):
        logger.error(f"The target directory does not exist: {base_dir}")
        return

    logger.info(f"Scanning target directory: {base_dir} for multiqc_data.json files...")
    
    found_reports = []
    for root, dirs, files in os.walk(base_dir):
        if target_filename in files:
            found_reports.append(os.path.join(root, target_filename))

    logger.info(f"Found {len(found_reports)} '{target_filename}' files to process.")

    success_count = 0
    
    for json_file_path in found_reports:
        path_parts = json_file_path.split(os.sep)
        
        try:
            jsons_index = path_parts.index("jsons")
            # If "jsons" is the immediate parent of the file, use the folder above it ("mq_data")
            if jsons_index == len(path_parts) - 2:
                dataset_name = path_parts[jsons_index - 1] # Grabs 'mq_data'
            else:
                dataset_name = path_parts[jsons_index + 1] # Grabs subfolder name
        except (ValueError, IndexError):
            dataset_name = path_parts[-2] if len(path_parts) > 1 else "unknown_dataset"

        output_dir = os.path.dirname(json_file_path)
        input_dir = output_dir 

        logger.info(f"Processing {dataset_name} (File path: {json_file_path})")
        
        # Build the exact token ID
        current_job_id = f"{dataset_name}_quality"

        try:
            write_mzqc(
                output_dir=output_dir,
                input_dir=input_dir,
                job_id=current_job_id,
                workflow="maxquant",
                json_path=json_file_path
            )
            
            # --- Robust Verification Check ---
            # We will verify using BOTH the exact job_id and the json fallback name 
            expected_file_a = os.path.join(output_dir, f"{current_job_id}.mzQC")
            expected_file_b = f"{json_file_path}_quality.mzQC"
            
            if os.path.exists(expected_file_a):
                logger.info(f" -> Successfully generated: {expected_file_a}")
                success_count += 1
            elif os.path.exists(expected_file_b):
                logger.info(f" -> Successfully generated (via fallback path): {expected_file_b}")
                success_count += 1
            else:
                logger.error(f" -> Failed to locate expected output file. Checked:\n    1: {expected_file_a}\n    2: {expected_file_b}")
                
        except Exception as e:
            logger.error(f" -> Critical error processing {dataset_name}: {e}", exc_info=True)

if __name__ == "__main__":
    run_batch_test()