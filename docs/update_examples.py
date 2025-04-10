import os
import json
import subprocess
import shutil
import zipfile
import gzip
import tarfile
from pathlib import Path


def download_file(url, save_path):
    print(f"Downloading {url} to {save_path} ...")
    subprocess.run(
        ["wget", "-nv", "-P", save_path, url],
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    print(f"File downloaded successfully: {os.path.basename(url)}")


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