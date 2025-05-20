import os
import json
import subprocess
import shutil
import time
import urllib
import zipfile
import gzip
import tarfile
from ftplib import FTP
from pathlib import Path
from urllib.parse import urlparse

from tqdm import tqdm


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
                sleep_time = backoff_factor ** attempt
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
                    reporthook.pbar = tqdm(total=total_size, unit="B", unit_scale=True, desc=filename)
                downloaded = block_num * block_size
                reporthook.pbar.update(downloaded - reporthook.pbar.n)

            reporthook.pbar = None

            urllib.request.urlretrieve(url, local_filepath, reporthook)
            if reporthook.pbar:
                reporthook.pbar.close()

            print(f"✅ HTTPS download complete: {filename}")
            return  # success
        except Exception as http_err:
            print(f"⚠️ HTTPS attempt {attempt} failed: {http_err}")
            if attempt < max_retries:
                sleep_time = backoff_factor ** attempt
                print(f"Retrying in {sleep_time} seconds...")
                time.sleep(sleep_time)
            else:
                print("❌ All download attempts failed.")
                raise http_err


def extract_zip(file_path, extract_to):
    with zipfile.ZipFile(file_path, "r") as zip_ref:
        zip_ref.extractall(extract_to)
        print(f"Extracted {file_path} to {extract_to}")


def extract_gz(file_path, extract_to):
    with gzip.open(file_path, "rb") as f_in:
        out_path = os.path.join(extract_to, os.path.basename(file_path).replace(".gz", ""))
        with open(out_path, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)
            print(f"Extracted {file_path} to {out_path}")


def extract_tar(file_path, extract_to):
    with tarfile.open(file_path, "r:*") as tar_ref:
        tar_ref.extractall(extract_to)
        print(f"Extracted {file_path} to {extract_to}")


def extract_files(folder_path):
    for root, _, files in os.walk(folder_path):
        for file in files:
            file_path = os.path.join(root, file)
            # *.zip
            if file.endswith(".zip"):
                extract_zip(file_path, root)
            # *.gz
            elif file.endswith(".gz") and not file.endswith(".tar.gz"):
                extract_gz(file_path, root)
            # .tar, .tar.gz, .tar.bz2
            elif (
                file.endswith(".tar")
                or file.endswith(".tar.gz")
                or file.endswith(".tgz")
                or file.endswith(".tar.bz2")
            ):
                extract_tar(file_path, root)


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
        command = ["multiqc", "--parse_maxquant", download_path, "-o", report_path]

    elif plugin_type == "mzid":
        command = ["multiqc", "--mzid_plugin", download_path, "-o", report_path]

    elif plugin_type == "dia" or plugin_type == "tmt" or plugin_type == "lfq":
        command = [
            "multiqc",
            download_path,
            "--config",
            os.path.join(download_path, "multiqc_config.yml"),
            "-o",
            report_path,
        ]
        
    elif plugin_type == "diann":
        command = ["multiqc", download_path, "-o", report_path]

    subprocess.run(command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)


def new_examples(config_file):

    if not os.path.exists(config_file):
        print(f"config file {config_file} not exist!")
        return

    with open(config_file, "r") as f:
        config = json.load(f)

    for project in config.get("projects"):

        urls = project.get("urls")
        report_path = project.get("path")
        file_type = project.get("file_type")

        download_path = "./data"
        Path(download_path).mkdir(parents=True, exist_ok=True)

        for url in urls:
            download_file(url, download_path)

        extract_files(download_path)

        Path(report_path).mkdir(parents=True, exist_ok=True)
        delete_old_examples(report_path)

        run_pmultiqc(download_path, report_path, file_type)

        shutil.rmtree(download_path)


new_examples("docs/config.json")