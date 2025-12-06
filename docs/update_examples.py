import argparse
import json
import os
import shutil
import subprocess
import time
import urllib
import urllib.request
from ftplib import FTP
from pathlib import Path
from urllib.parse import urlparse

from tqdm import tqdm

from pmultiqc.modules.common.file_utils import extract_files


def download_file(url, save_path, max_retries=3, backoff_factor=2):
    parsed_url = urlparse(url)
    filename = os.path.basename(parsed_url.path)
    local_filepath = os.path.join(save_path, filename)
    Path(save_path).mkdir(parents=True, exist_ok=True)

    ftp_url = url.replace("https://", "ftp://", 1)
    ftp_parsed = urlparse(ftp_url)
    ftp_host = ftp_parsed.hostname
    ftp_path = ftp_parsed.path

    # -------- Try FTP first --------
    for attempt in range(1, max_retries + 1):
        try:
            print(f"[Attempt {attempt}] Trying FTP download: {ftp_host}{ftp_path}")
            ftp = FTP(ftp_host, timeout=30)
            ftp.login()
            total_size = ftp.size(ftp_path)

            with open(local_filepath, "wb") as f, tqdm(
                total=total_size, unit="B", unit_scale=True, desc=filename
            ) as pbar:

                def callback(data):
                    f.write(data)
                    pbar.update(len(data))

                ftp.retrbinary(f"RETR {ftp_path}", callback)

            ftp.quit()
            print(f"✅ FTP download complete: {filename}")
            return  # success
        except Exception as ftp_err:
            print(f"⚠️ FTP attempt {attempt} failed: {ftp_err}")
            if attempt < max_retries:
                sleep_time = backoff_factor**attempt
                print(f"Retrying in {sleep_time} seconds...")
                time.sleep(sleep_time)
            else:
                print("FTP download failed after maximum retries. Trying HTTPS...")

    # -------- Try HTTPS fallback --------
    for attempt in range(1, max_retries + 1):
        try:
            print(f"[Attempt {attempt}] Trying HTTPS download: {url}")

            def reporthook(block_num, block_size, total_size):
                if reporthook.pbar is None:
                    reporthook.pbar = tqdm(
                        total=total_size, unit="B", unit_scale=True, desc=filename
                    )
                downloaded = block_num * block_size
                reporthook.pbar.update(downloaded - reporthook.pbar.n)

            reporthook.pbar = None

            # Try with default SSL context first
            try:
                urllib.request.urlretrieve(url, local_filepath, reporthook)
            except urllib.error.URLError as ssl_err:
                import ssl

                # Check if this is specifically an SSL certificate error
                if isinstance(ssl_err.reason, ssl.SSLError) and "CERTIFICATE_VERIFY_FAILED" in str(
                    ssl_err.reason
                ):
                    print(
                        "⚠️ SSL certificate verification failed, trying with SSL context that bypasses verification..."
                    )
                    print(
                        "🔓 WARNING: This bypasses SSL certificate verification and reduces security!"
                    )

                    # Store original opener to restore later
                    original_opener = urllib.request._opener

                    try:
                        # Create SSL context that bypasses certificate verification
                        ssl_context = ssl.create_default_context()
                        ssl_context.check_hostname = False
                        ssl_context.verify_mode = ssl.CERT_NONE

                        # Create opener with custom SSL context
                        opener = urllib.request.build_opener(
                            urllib.request.HTTPSHandler(context=ssl_context)
                        )
                        urllib.request.install_opener(opener)

                        # Retry with SSL context that bypasses verification
                        urllib.request.urlretrieve(url, local_filepath, reporthook)
                    finally:
                        # Restore original opener to limit scope of changes
                        if original_opener:
                            urllib.request.install_opener(original_opener)
                        else:
                            urllib.request.install_opener(urllib.request.build_opener())
                else:
                    raise ssl_err
            if reporthook.pbar:
                reporthook.pbar.close()

            print(f"✅ HTTPS download complete: {filename}")
            return  # success
        except Exception as http_err:
            print(f"⚠️ HTTPS attempt {attempt} failed: {http_err}")
            if attempt < max_retries:
                sleep_time = backoff_factor**attempt
                print(f"Retrying in {sleep_time} seconds...")
                time.sleep(sleep_time)
            else:
                print("❌ All download attempts failed.")
                raise http_err


def delete_old_examples(folder_path):
    for filename in os.listdir(folder_path):
        file_path = os.path.join(folder_path, filename)

        if os.path.isfile(file_path):
            os.remove(file_path)
            print(f"Deleted file: {file_path}")

        elif os.path.isdir(file_path):
            shutil.rmtree(file_path)
            print(f"Deleted folder: {file_path}")


def run_pmultiqc(download_path, report_path, plugin_type):

    if plugin_type == "maxquant":
        command = ["multiqc", "--maxquant_plugin", download_path, "-o", report_path]

    elif plugin_type == "proteobench":
        command = ["multiqc", "--proteobench_plugin", download_path, "-o", report_path]

    elif plugin_type == "mzid":
        command = ["multiqc", "--mzid_plugin", download_path, "-o", report_path]

    elif plugin_type == "dia" or plugin_type == "tmt" or plugin_type == "lfq":
        command = [
            "multiqc",
            "--quantms_plugin",
            download_path,
            "--config",
            os.path.join(download_path, "multiqc_config.yml"),
            "-o",
            report_path,
        ]
    elif plugin_type == "diann":
        command = ["multiqc", "--diann_plugin", download_path, "-o", report_path]

    subprocess.run(command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)



def new_examples(config_file, project_accession=None):

    if not os.path.exists(config_file):
        print(f"Error: config file {config_file} not exist!")
        return

    with open(config_file, "r") as f:
        config = json.load(f)

    projects = config.get("projects", [])

    # Filter to specific project if accession is provided
    if project_accession:
        projects = [p for p in projects if p.get("accession") == project_accession]
        if not projects:
            print(f"❌ Project with accession '{project_accession}' not found in config!")
            return
        print(f"✅ Processing project: {project_accession}")

    for project in projects:
        accession = project.get("accession")
        urls = project.get("urls")
        report_path = project.get("path")
        file_type = project.get("file_type")

        print(f"\n{'='*60}")
        print(f"Processing project: {accession}")
        print(f"Report path: {report_path}")
        print(f"File type: {file_type}")
        print(f"{'='*60}\n")

        download_path = "./data_temp"

        if os.path.exists(download_path):
            shutil.rmtree(download_path)
        Path(download_path).mkdir(parents=True, exist_ok=True)

        try:
            for url in urls:
                print(f"Downloading: {url}")
                download_file(url, download_path)

            extract_files(download_path)

            Path(report_path).mkdir(parents=True, exist_ok=True)
            delete_old_examples(report_path)

            print("Running pmultiqc...")
            run_pmultiqc(download_path, report_path, file_type)

            print(f"✅ Successfully processed project: {accession}")
        except Exception as e:
            print(f"❌ Error processing project {accession}: {e}")
            raise e
        finally:
            # Clean up download directory
            if os.path.exists(download_path):
                shutil.rmtree(download_path)


def main():
    parser = argparse.ArgumentParser(description="Update pmultiqc examples")
    parser.add_argument(
        "--project",
        type=str,
        help="Process only a specific project by accession (e.g., 'LFQ_PXD007683')",
    )
    parser.add_argument(
        "--config",
        type=str,
        default="docs/config.json",
        help="Path to config.json file (default: docs/config.json)",
    )
    args = parser.parse_args()

    new_examples(args.config, args.project)


if __name__ == "__main__":
    main()