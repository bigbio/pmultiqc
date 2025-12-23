"""
Disk cleanup script for pmultiqc service
Removes old job directories and corrupted files to free up disk space
"""

import os
import shutil
import time
import uuid as uuid_lib
import zipfile


def cleanup_job_directories(base_dir, max_age_days=7):
    """Clean up old job directories"""
    if not os.path.exists(base_dir):
        print(f"Directory {base_dir} does not exist")
        return 0

    total_cleaned = 0
    cutoff_time = time.time() - (max_age_days * 24 * 3600)

    for item in os.listdir(base_dir):
        item_path = os.path.join(base_dir, item)
        # Security: Validate that it's a directory and has UUID format (8-4-4-4-12 format)
        if os.path.isdir(item_path) and len(item) == 36:
            # Additional validation: Check if it looks like a valid UUID
            try:
                uuid_lib.UUID(item)  # This will raise ValueError if not a valid UUID
            except ValueError:
                print(f"Skipping non-UUID directory: {item}")
                continue
            
            try:
                mtime = os.path.getmtime(item_path)
                if mtime < cutoff_time:
                    size = sum(
                        os.path.getsize(os.path.join(dirpath, filename))
                        for dirpath, dirnames, filenames in os.walk(item_path)
                        for filename in filenames
                    )
                    print(f"Removing old job directory: {item} ({size / (1024**3):.2f} GB)")
                    shutil.rmtree(item_path)
                    total_cleaned += size
            except Exception as e:
                print(f"Error cleaning up {item}: {e}")

    return total_cleaned


def cleanup_corrupted_zips(base_dir):
    """Clean up corrupted zip files"""
    if not os.path.exists(base_dir):
        return 0

    total_cleaned = 0

    for root, dirs, files in os.walk(base_dir):
        for file in files:
            if file.endswith(".zip"):
                file_path = os.path.join(root, file)
                try:
                    # Check if zip is valid
                    with zipfile.ZipFile(file_path, "r") as zipf:
                        zipf.testzip()
                except zipfile.BadZipFile:
                    size = os.path.getsize(file_path)
                    print(f"Removing corrupted zip: {file_path} ({size / (1024**3):.2f} GB)")
                    os.remove(file_path)
                    total_cleaned += size
                except Exception as e:
                    print(f"Error checking zip {file_path}: {e}")

    return total_cleaned


def main():
    # Configuration
    upload_dir = "local_uploads"
    output_dir = "local_outputs"
    html_reports_dir = "local_html_reports"

    print("pmultiqc Disk Cleanup Script")
    print("=" * 40)

    # Clean up old job directories
    print("\n1. Cleaning up old job directories...")
    upload_cleaned = cleanup_job_directories(upload_dir)
    output_cleaned = cleanup_job_directories(output_dir)
    html_cleaned = cleanup_job_directories(html_reports_dir)

    # Clean up corrupted zip files
    print("\n2. Cleaning up corrupted zip files...")
    zip_cleaned = cleanup_corrupted_zips(output_dir)

    # Summary
    total_cleaned = upload_cleaned + output_cleaned + html_cleaned + zip_cleaned
    print("\n" + "=" * 40)
    print("Cleanup Summary:")
    print(f"Upload directories cleaned: {upload_cleaned / (1024**3):.2f} GB")
    print(f"Output directories cleaned: {output_cleaned / (1024**3):.2f} GB")
    print(f"HTML reports cleaned: {html_cleaned / (1024**3):.2f} GB")
    print(f"Corrupted zips cleaned: {zip_cleaned / (1024**3):.2f} GB")
    print(f"Total space freed: {total_cleaned / (1024**3):.2f} GB")

    # Show current disk usage
    print("\nCurrent disk usage:")
    for dir_name in [upload_dir, output_dir, html_reports_dir]:
        if os.path.exists(dir_name):
            total_size = sum(
                os.path.getsize(os.path.join(dirpath, filename))
                for dirpath, dirnames, filenames in os.walk(dir_name)
                for filename in filenames
            )
            print(f"{dir_name}: {total_size / (1024**3):.2f} GB")


if __name__ == "__main__":
    main()