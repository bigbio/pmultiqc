#!/usr/bin/env python
"""
PMultiQC Service API
A Flask-based web service for generating PMultiQC reports from uploaded data files.
"""

import os
import zipfile
import shutil
import uuid
import logging
import threading
import sys
import requests
import json
import time
import subprocess
import hashlib
import traceback
from datetime import datetime
from typing import Dict, Any, List, Optional

from flask import Flask, request, jsonify, send_file, render_template, redirect
from flask_cors import CORS
from werkzeug.utils import secure_filename

# Try to import yaml, but make it optional
try:
    import yaml
    YAML_AVAILABLE = True
except ImportError:
    YAML_AVAILABLE = False

# Configure logging for Kubernetes
def setup_logging():
    """Configure logging for Kubernetes environment."""
    # Get log level from environment or default to INFO
    log_level = os.environ.get('LOG_LEVEL', 'INFO').upper()
    
    # Create a formatter that includes timestamp, level, and module
    formatter = logging.Formatter(
        '%(asctime)s [%(levelname)s] %(name)s: %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    
    # Create console handler that writes to stdout
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(getattr(logging, log_level))
    console_handler.setFormatter(formatter)
    
    # Configure root logger
    root_logger = logging.getLogger()
    root_logger.setLevel(getattr(logging, log_level))
    
    # Remove existing handlers to avoid duplicates
    for handler in root_logger.handlers[:]:
        root_logger.removeHandler(handler)
    
    # Add our console handler
    root_logger.addHandler(console_handler)
    
    # Create application logger
    logger = logging.getLogger(__name__)
    logger.info(f"Logging configured with level: {log_level}")
    
    return logger

# Initialize logging
logger = setup_logging()

app = Flask(__name__, template_folder='templates', static_folder='templates')
CORS(app)

# Configuration with environment variable and YAML support
def load_yaml_config(config_path: str = None) -> Dict[str, Any]:
    """Load configuration from YAML file if available."""
    if not YAML_AVAILABLE:
        logger.warning("PyYAML not available, skipping YAML config file")
        return {}
    
    # Try to find config file
    config_paths = [
        config_path,  # Explicitly provided path
        os.environ.get('CONFIG_FILE'),  # From environment variable
        'config.yaml',  # Current directory
        'config.yml',   # Alternative extension
        '/etc/pmultiqc/config.yaml',  # System-wide config
        os.path.join(os.path.dirname(__file__), 'config.yaml')  # Next to app.py
    ]
    
    for path in config_paths:
        if path and os.path.exists(path):
            try:
                with open(path, 'r') as f:
                    config = yaml.safe_load(f)
                logger.info(f"Loaded configuration from: {path}")
                return config
            except Exception as e:
                logger.error(f"Error loading config from {path}: {e}")
    
    logger.info("No YAML config file found, using environment variables and defaults")
    return {}

def get_config_value(env_var: str, yaml_key: str, default_value: str, yaml_config: Dict[str, Any]) -> str:
    """Get configuration value from environment variable, YAML config, or use default."""
    # Priority: Environment variable > YAML config > default
    env_value = os.environ.get(env_var)
    if env_value is not None:
        logger.debug(f"Using environment variable {env_var}: {env_value}")
        return env_value
    
    yaml_value = yaml_config.get(yaml_key)
    if yaml_value is not None:
        logger.debug(f"Using YAML config {yaml_key}: {yaml_value}")
        return str(yaml_value)
    
    logger.debug(f"Using default value for {yaml_key}: {default_value}")
    return default_value

def get_config_int(env_var: str, yaml_key: str, default_value: int, yaml_config: Dict[str, Any]) -> int:
    """Get integer configuration value from environment variable, YAML config, or use default."""
    # Priority: Environment variable > YAML config > default
    env_value = os.environ.get(env_var)
    if env_value is not None:
        try:
            value = int(env_value)
            logger.debug(f"Using environment variable {env_var}: {value}")
            return value
        except ValueError:
            logger.warning(f"Invalid integer value for {env_var}: {env_value}, trying YAML config")
    
    yaml_value = yaml_config.get(yaml_key)
    if yaml_value is not None:
        try:
            value = int(yaml_value)
            logger.debug(f"Using YAML config {yaml_key}: {value}")
            return value
        except ValueError:
            logger.warning(f"Invalid integer value for {yaml_key} in YAML: {yaml_value}, using default")
    
    logger.debug(f"Using default value for {yaml_key}: {default_value}")
    return default_value

# Load YAML configuration
yaml_config = load_yaml_config()

# Configuration
app.config['MAX_CONTENT_LENGTH'] = get_config_int('MAX_CONTENT_LENGTH', 'max_content_length', 2 * 1024 * 1024 * 1024, yaml_config)  # 2GB max file size
app.config['UPLOAD_FOLDER'] = get_config_value('UPLOAD_FOLDER', 'upload_folder', '/tmp/pmultiqc_uploads', yaml_config)
app.config['OUTPUT_FOLDER'] = get_config_value('OUTPUT_FOLDER', 'output_folder', '/tmp/pmultiqc_outputs', yaml_config)
app.config['HTML_REPORTS_FOLDER'] = get_config_value('HTML_REPORTS_FOLDER', 'html_reports_folder', '/tmp/pmultiqc_html_reports', yaml_config)
app.config['BASE_URL'] = get_config_value('BASE_URL', 'base_url', 'http://localhost:5000', yaml_config)

# Log configuration values
logger.info(f"Configuration loaded:")
logger.info(f"  MAX_CONTENT_LENGTH: {app.config['MAX_CONTENT_LENGTH']} bytes ({app.config['MAX_CONTENT_LENGTH'] / (1024*1024*1024):.1f} GB)")
logger.info(f"  UPLOAD_FOLDER: {app.config['UPLOAD_FOLDER']}")
logger.info(f"  OUTPUT_FOLDER: {app.config['OUTPUT_FOLDER']}")
logger.info(f"  HTML_REPORTS_FOLDER: {app.config['HTML_REPORTS_FOLDER']}")
logger.info(f"  BASE_URL: {app.config['BASE_URL']}")
if yaml_config:
    logger.info(f"  YAML config keys: {list(yaml_config.keys())}")

# Ensure directories exist
os.makedirs(app.config['UPLOAD_FOLDER'], exist_ok=True)
os.makedirs(app.config['OUTPUT_FOLDER'], exist_ok=True)
os.makedirs(app.config['HTML_REPORTS_FOLDER'], exist_ok=True)

# Job status tracking
job_status_dict = {}

# Allowed file extensions
ALLOWED_EXTENSIONS = {'zip'}

# Files with completion time exceeding 24 hours will be cleaned up.
EXPIRATION_SECONDS = 24 * 3600
CLEANUP_JOB_SECONDS = 3600

def allowed_file(filename):
    """Check if the uploaded file has an allowed extension."""
    return '.' in filename and filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS

def current_time():
    now = int(time.time())
    dt = datetime.fromtimestamp(now)
    return now, dt.strftime("%Y-%m-%d %H:%M:%S")

def get_pride_files(accession: str) -> List[Dict]:
    """
    Fetch files from PRIDE API for a given accession.
    
    Args:
        accession: PRIDE accession (e.g., PXD039077)
    
    Returns:
        List of file dictionaries from PRIDE API
    """
    try:
        url = f"https://www.ebi.ac.uk/pride/ws/archive/v3/projects/{accession}/files"
        params = {"page": 0, "pageSize": 100}  # Get up to 100 files
        
        logger.info(f"Fetching files from PRIDE for accession: {accession}")
        response = requests.get(url, params=params, timeout=30)
        response.raise_for_status()
        
        files = response.json()
        logger.info(f"Retrieved {len(files)} files from PRIDE for {accession}")
        return files
        
    except requests.RequestException as e:
        logger.error(f"Error fetching PRIDE files for {accession}: {e}")
        raise Exception(f"Failed to fetch files from PRIDE: {e}")

def filter_search_files(files: List[Dict]) -> List[Dict]:
    """
    Filter files to only include search engine output files.
    
    Args:
        files: List of file dictionaries from PRIDE API
    
    Returns:
        List of search engine output files
    """
    search_files = []
    
    for file_info in files:
        # Check if it's a search engine output file
        if (file_info.get('fileCategory', {}).get('value') == 'SEARCH' or
            'report' in file_info.get('fileName', '').lower() or
            'evidence' in file_info.get('fileName', '').lower() or
            'peptides' in file_info.get('fileName', '').lower() or
            'proteingroups' in file_info.get('fileName', '').lower() or
            'msms' in file_info.get('fileName', '').lower()):
            
            search_files.append(file_info)
    
    logger.info(f"Found {len(search_files)} search engine output files")
    return search_files

def download_pride_file(file_info: Dict, download_dir: str) -> str:
    """
    Download a file from PRIDE.
    
    Args:
        file_info: File information from PRIDE API
        download_dir: Directory to save the file
    
    Returns:
        Path to the downloaded file
    """
    try:
        # Get FTP URL (preferred over Aspera for reliability)
        ftp_url = None
        for location in file_info.get('publicFileLocations', []):
            if location.get('accession') == 'PRIDE:0000469':  # FTP Protocol
                ftp_url = location.get('value')
                break
        
        if not ftp_url:
            raise Exception(f"No FTP URL found for file {file_info.get('fileName')}")
        
        # Convert FTP URL to HTTP URL for easier download
        http_url = ftp_url.replace('ftp://', 'https://')
        
        filename = file_info.get('fileName')
        file_path = os.path.join(download_dir, filename)
        
        logger.info(f"Downloading {filename} from {http_url}")
        
        # Download file with progress tracking
        response = requests.get(http_url, stream=True, timeout=300)
        response.raise_for_status()
        
        with open(file_path, 'wb') as f:
            for chunk in response.iter_content(chunk_size=8192):
                if chunk:
                    f.write(chunk)
        
        logger.info(f"Successfully downloaded {filename} ({os.path.getsize(file_path)} bytes)")
        return file_path
        
    except Exception as e:
        logger.error(f"Error downloading file {file_info.get('fileName')}: {e}")
        raise Exception(f"Failed to download {file_info.get('fileName')}: {e}")

def detect_pride_input_type(files: List[Dict]) -> str:
    """
    Detect input type based on PRIDE files.
    
    Args:
        files: List of file dictionaries from PRIDE API
    
    Returns:
        Input type string
    """
    file_names = [f.get('fileName', '').lower() for f in files]
    
    # Check for MaxQuant files
    maxquant_files = ['evidence.txt', 'msms.txt', 'peptides.txt', 'proteingroups.txt']
    if any(any(mq_file in fname for mq_file in maxquant_files) for fname in file_names):
        return "maxquant"
    
    # Check for DIANN files
    diann_files = ['report.tsv', 'report.parquet']
    if any(any(diann_file in fname for diann_file in diann_files) for fname in file_names):
        return "diann"
    
    # Check for quantms files
    quantms_files = ['ms_info.parquet', '.mztab']
    if any(any(qt_file in fname for qt_file in quantms_files) for fname in file_names):
        return "quantms"
    
    # Check for mzIdentML files
    mzid_files = ['.mzid', '.mzml']
    if any(any(mzid_file in fname for mzid_file in mzid_files) for fname in file_names):
        return "mzidentml"
    
    return "unknown"

# List all files in the extracted directory
def scan_files(upload_path):
    files = []
    file_paths = {}

    for root, _, filenames in os.walk(upload_path):
        for filename in filenames:
            full_path = os.path.join(root, filename)
            name_lower = filename.lower()
            files.append(name_lower)
            file_paths[name_lower] = full_path

    return files, file_paths

def detect_input_type(upload_path: str) -> str:
    """
    Detect the type of input data based on the contents of the extracted zip file.
    
    Returns:
        str: One of 'maxquant', 'diann', 'quantms', 'mzidentml', or 'unknown'
    """
    try:
        files, file_paths = scan_files(upload_path)

        # Check for quantms files
        quantms_files = ["_ms_info.parquet", ".mzTab"]
        if "multiqc_config.yml" in files and any(any(substr in f for substr in quantms_files) for f in files):
            return "quantms", file_paths["multiqc_config.yml"]

        # Check for MaxQuant files
        maxquant_files = ['evidence.txt', 'msms.txt', 'peptides.txt', 'proteinGroups.txt']
        if any(f in files for f in maxquant_files):
            return "maxquant", None
        
        # Check for DIANN files
        diann_files = ['report.tsv', 'report.parquet']
        if any(f in files for f in diann_files) or any('diann_report' in f for f in files):
            return "diann", None
        
        # Check for mzIdentML + mzML files
        mzid_files = ['.mzid', '.mzml']
        if any(any(f.endswith(ext) for ext in mzid_files) for f in files):
            return "mzidentml", None
        
        return "unknown", None
    
    except Exception as e:
        logger.error(f"Error detecting input type: {e}")
        return "unknown", None

def run_pmultiqc(input_path: str, output_path: str, input_type: str, pmultiqc_config: str) -> Dict[str, Any]:
    """
    Run PMultiQC on the input data and generate a report.
    
    Args:
        input_path: Path to the input data directory
        output_path: Path where the output report should be saved
        input_type: Type of input data ('maxquant', 'diann', 'quantms', 'mzidentml')
    
    Returns:
        Dict containing success status and any error messages
    """
    try:
        # Set up MultiQC arguments based on input type
        args = [
            'multiqc',
            input_path,
            '-o', output_path,
            '--force'
        ]
        
        # Add type-specific arguments
        if input_type == 'maxquant':
            args.extend(['--parse_maxquant'])
        elif input_type == "quantms":
            args.extend(["--config", pmultiqc_config])
        elif input_type == 'diann':
            # DIANN files are handled automatically by pmultiqc
            pass
        elif input_type == 'mzidentml':
            args.extend(['--mzid_plugin'])
        
        # Run MultiQC with PMultiQC plugin
        logger.info(f"Running PMultiQC with args: {args}")
        
        # Run multiqc as a subprocess with real-time output
        try:
            # Run the multiqc command with real-time output
            process = subprocess.Popen(
                args, 
                stdout=subprocess.PIPE, 
                stderr=subprocess.PIPE, 
                text=True, 
                bufsize=1,
                universal_newlines=True
            )
            
            # Stream output in real-time
            output_lines = []
            error_lines = []
            
            while True:
                output = process.stdout.readline()
                error = process.stderr.readline()
                
                if output:
                    output_lines.append(output.strip())
                    logger.info(f"MultiQC: {output.strip()}")
                
                if error:
                    error_lines.append(error.strip())
                    logger.warning(f"MultiQC: {error.strip()}")
                
                # Check if process has finished
                if process.poll() is not None:
                    break
            
            # Get any remaining output
            remaining_output, remaining_error = process.communicate()
            if remaining_output:
                output_lines.extend(remaining_output.strip().split('\n'))
                logger.info(f"MultiQC remaining output: {remaining_output.strip()}")
            if remaining_error:
                error_lines.extend(remaining_error.strip().split('\n'))
                logger.warning(f"MultiQC remaining error: {remaining_error.strip()}")
            
            if process.returncode == 0:
                return {
                    'success': True, 
                    'message': 'Report generated successfully',
                    'output': output_lines,
                    'errors': error_lines
                }
            else:
                return {
                    'success': False, 
                    'message': f'MultiQC failed with return code {process.returncode}',
                    'output': output_lines,
                    'errors': error_lines
                }
                
        except subprocess.CalledProcessError as e:
            logger.error(f"MultiQC failed with return code {e.returncode}")
            logger.error(f"MultiQC stdout: {e.stdout}")
            logger.error(f"MultiQC stderr: {e.stderr}")
            return {
                'success': False, 
                'message': f'MultiQC failed: {e.stderr}',
                'output': e.stdout.split('\n') if e.stdout else [],
                'errors': e.stderr.split('\n') if e.stderr else []
            }
            
    except Exception as e:
        logger.error(f"Error running PMultiQC: {e}")
        return {'success': False, 'message': str(e), 'output': [], 'errors': []}

def create_zip_report(output_path: str, zip_path: str) -> bool:
    """
    Create a zip file containing the PMultiQC report.
    
    Args:
        output_path: Path to the MultiQC output directory
        zip_path: Path where the zip file should be created
    
    Returns:
        bool: True if successful, False otherwise
    """
    try:
        logger.info(f"Creating zip report from {output_path} to {zip_path}")
        
        # Check if output directory exists and has files
        if not os.path.exists(output_path):
            logger.error(f"Output path does not exist: {output_path}")
            return False
        
        files_found = []
        for root, dirs, files in os.walk(output_path):
            for file in files:
                files_found.append(os.path.join(root, file))
        
        logger.info(f"Found {len(files_found)} files to zip")
        
        if not files_found:
            logger.error("No files found to zip")
            return False
        
        # Create zip file with better error handling
        with zipfile.ZipFile(zip_path, 'w', zipfile.ZIP_DEFLATED) as zipf:
            for file_path in files_found:
                try:
                    arcname = os.path.relpath(file_path, output_path)
                    logger.info(f"Adding {file_path} as {arcname}")
                    zipf.write(file_path, arcname)
                except Exception as e:
                    logger.error(f"Error adding file {file_path} to zip: {e}")
                    # Continue with other files instead of failing completely
                    continue
        
        # Check if zip file was created and has content
        if os.path.exists(zip_path) and os.path.getsize(zip_path) > 0:
            logger.info(f"Successfully created zip file: {zip_path} ({os.path.getsize(zip_path)} bytes)")
            return True
        else:
            logger.error(f"Zip file was not created or is empty: {zip_path}")
            return False
            
    except Exception as e:
        logger.error(f"Error creating zip report: {e}")
        return False

def copy_html_report_for_online_viewing(output_path: str, job_id: str, accession: str = None) -> Optional[Dict[str, str]]:
    """
    Copy HTML reports to public directory for online viewing.
    
    Args:
        output_path: Path to the MultiQC output directory
        job_id: Job ID for the report
        accession: PRIDE accession number (optional, for PRIDE datasets)
    
    Returns:
        dict: Dictionary mapping report names to their URLs, or None if failed
    """
    try:
        logger.info(f"Starting copy_html_report_for_online_viewing: output_path={output_path}, job_id={job_id}, accession={accession}")
        
        # Log the contents of the output directory first
        logger.info(f"Contents of output directory {output_path}:")
        if os.path.exists(output_path):
            for root, dirs, files in os.walk(output_path):
                level = root.replace(output_path, '').count(os.sep)
                indent = ' ' * 2 * level
                logger.info(f"{indent}{os.path.basename(root)}/")
                subindent = ' ' * 2 * (level + 1)
                for file in files:
                    logger.info(f"{subindent}{file}")
        else:
            logger.error(f"Output directory does not exist: {output_path}")
            return None
        
        # Look for all HTML report files
        html_reports = []
        for root, dirs, files in os.walk(output_path):
            for file in files:
                if file == 'multiqc_report.html':
                    html_report_path = os.path.join(root, file)
                    # Get the relative path from output_path to determine report name
                    rel_path = os.path.relpath(root, output_path)
                    # Use "." as root if the report is directly in output_path
                    if rel_path == '.':
                        rel_path = 'root'
                    html_reports.append((rel_path, html_report_path))
                    logger.info(f"Found HTML report: {html_report_path} (relative path: {rel_path})")
        
        if not html_reports:
            logger.warning(f"No HTML reports found in {output_path}")
            logger.info("Searching for any .html files in the output directory...")
            for root, dirs, files in os.walk(output_path):
                for file in files:
                    if file.endswith('.html'):
                        logger.info(f"Found HTML file: {os.path.join(root, file)}")
            return None
        
        # Create a unique hash for the report URL
        hash_input = f"{job_id}_{accession}" if accession else job_id
        report_hash = hashlib.md5(hash_input.encode()).hexdigest()[:12]
        logger.info(f"Generated report hash: {report_hash} (from input: {hash_input})")
        
        # Create public directory for this report
        public_dir = os.path.join(app.config['HTML_REPORTS_FOLDER'], report_hash)
        os.makedirs(public_dir, exist_ok=True)
        logger.info(f"Created public directory: {public_dir}")
        
        # Create a job file for cleanup tracking
        job_file_path = os.path.join(public_dir, '.job_id')
        with open(job_file_path, 'w') as f:
            f.write(job_id)
        
        # Copy all files from output directory to public directory
        files_copied = 0
        for root, dirs, files in os.walk(output_path):
            for file in files:
                src_path = os.path.join(root, file)
                rel_path = os.path.relpath(src_path, output_path)
                dst_path = os.path.join(public_dir, rel_path)
                
                # Ensure destination directory exists
                os.makedirs(os.path.dirname(dst_path), exist_ok=True)
                
                # Copy file
                shutil.copy2(src_path, dst_path)
                files_copied += 1
                logger.info(f"Copied {src_path} to {dst_path}")
        
        logger.info(f"Total files copied: {files_copied}")
        
        # Create URLs for each HTML report
        html_urls = {}
        for report_name, html_report_path in html_reports:
            # Find the relative path to the HTML report for the URL
            html_rel_path = os.path.relpath(html_report_path, output_path)
            # Use base URL to generate absolute URLs
            base_url = app.config['BASE_URL'].rstrip('/')
            html_url = f'{base_url}/reports/{report_hash}/{html_rel_path}'
            html_urls[report_name] = html_url
            logger.info(f"HTML report for '{report_name}' available at: {html_url}")
        
        logger.info(f"Generated HTML URLs: {html_urls}")
        return html_urls
        
    except Exception as e:
        logger.error(f"Error copying HTML reports for PRIDE: {e}")
        logger.error(f"Traceback: {traceback.format_exc()}")
        return None

def save_job_data_to_disk(job_id: str, job_data: dict, output_dir: str):
    """
    Save job data to disk for future retrieval when job is removed from memory.
    
    Args:
        job_id: The job ID
        job_data: The complete job data dictionary
        output_dir: The output directory for this job
    """
    try:
        # Save console output
        console_output = job_data.get('console_output', [])
        if console_output:
            console_output_file = os.path.join(output_dir, 'console_output.txt')
            with open(console_output_file, 'w', encoding='utf-8') as f:
                for line in console_output:
                    f.write(line + '\n')
            logger.info(f"Saved {len(console_output)} console output lines for job {job_id}")
        
        # Save console errors
        console_errors = job_data.get('console_errors', [])
        if console_errors:
            console_errors_file = os.path.join(output_dir, 'console_errors.txt')
            with open(console_errors_file, 'w', encoding='utf-8') as f:
                for line in console_errors:
                    f.write(line + '\n')
            logger.info(f"Saved {len(console_errors)} console error lines for job {job_id}")
        
        # Save job info
        job_info = {
            'job_id': job_id,
            'status': job_data.get('status'),
            'input_type': job_data.get('input_type'),
            'accession': job_data.get('accession'),
            'results': job_data.get('results', []),
            'finished_at': job_data.get('finished_at'),
            'files_processed': job_data.get('files_processed'),
            'total_files': job_data.get('total_files'),
            'files_downloaded': job_data.get('files_downloaded'),
            'html_report_urls': job_data.get('html_report_urls')
        }
        
        job_info_file = os.path.join(output_dir, 'job_info.json')
        with open(job_info_file, 'w', encoding='utf-8') as f:
            json.dump(job_info, f, indent=2)
        logger.info(f"Saved job info for job {job_id} to {job_info_file}")
        
    except Exception as e:
        logger.error(f"Error saving job data to disk for job {job_id}: {e}")

# Scheduled cleanup of the tmp folder
def cleanup_expired_jobs():

    now = current_time()[0]
    to_delete_jobids = []

    for job_id, info in job_status_dict.items():

        finished_at_seconds = info.get('finished_at_seconds')

        if finished_at_seconds is None:
            continue

        if now - finished_at_seconds > EXPIRATION_SECONDS:
            cleanup_job_files(job_id, now)
            to_delete_jobids.append(job_id)

    for job_id in to_delete_jobids:
        job_data = job_status_dict.pop(job_id, None)
        if job_data:
            logger.info(f"Job ID {job_id} has been deleted from job_status_dict at: {now} (was {job_data.get('status', 'unknown')})")
        else:
            logger.warning(f"Job ID {job_id} was already missing from job_status_dict at: {now}")

    threading.Timer(CLEANUP_JOB_SECONDS, cleanup_expired_jobs).start()

def cleanup_job_files(job_id, now_time):
    try:
        upload_dir = os.path.join(app.config['UPLOAD_FOLDER'], job_id)
        output_dir = os.path.join(app.config['OUTPUT_FOLDER'], job_id)

        for directory in [upload_dir, output_dir]:
            if os.path.exists(directory):
                shutil.rmtree(directory)
                logger.info(f"Removed directory {directory} for job {job_id} at {now_time}")
        
        # Clean up HTML reports for PRIDE datasets
        if os.path.exists(app.config['HTML_REPORTS_FOLDER']):
            for report_dir in os.listdir(app.config['HTML_REPORTS_FOLDER']):
                report_path = os.path.join(app.config['HTML_REPORTS_FOLDER'], report_dir)
                if os.path.isdir(report_path):
                    job_file_path = os.path.join(report_path, '.job_id')
                    if os.path.exists(job_file_path):
                        with open(job_file_path, 'r') as f:
                            stored_job_id = f.read().strip()
                            if stored_job_id == job_id:
                                shutil.rmtree(report_path)
                                logger.info(f"Removed HTML report directory {report_path} for job {job_id} at {now_time}")

    except Exception as e:
        logger.error(f"Error cleaning files for job {job_id}: {e}")

if __name__ == '__main__':
    cleanup_expired_jobs()


@app.route('/', methods=['GET'])
def index():
    """Main page with upload interface."""
    return render_template('index.html')

@app.route('/health', methods=['GET'])
def health_check():
    """Health check endpoint."""
    return jsonify({
        'status': 'healthy',
        'timestamp': datetime.now().isoformat(),
        'service': 'pmultiqc-service'
    })


@app.route('/generate-report', methods=['POST'])
def generate_report():
    """
    Generate a PMultiQC report from uploaded data.
    
    Expected form data:
    - file: ZIP file containing the data to analyze
    
    Returns:
    - JSON response with job ID and status
    """
    try:
        # Check if file was uploaded
        if 'file' not in request.files:
            return jsonify({'error': 'No file provided'}), 400
        
        file = request.files['file']
        if file.filename == '':
            return jsonify({'error': 'No file selected'}), 400
        
        if not allowed_file(file.filename):
            return jsonify({'error': 'Only ZIP files are allowed'}), 400
        
        # Generate unique job ID
        job_id = str(uuid.uuid4())
        
        # Create job directories
        upload_dir = os.path.join(app.config['UPLOAD_FOLDER'], job_id)
        output_dir = os.path.join(app.config['OUTPUT_FOLDER'], job_id)
        os.makedirs(upload_dir, exist_ok=True)
        os.makedirs(output_dir, exist_ok=True)
        
        # Save uploaded file
        filename = secure_filename(file.filename)
        zip_path = os.path.join(upload_dir, filename)
        file.save(zip_path)
        
        # Extract the zip file
        extract_path = os.path.join(upload_dir, 'extracted')
        os.makedirs(extract_path, exist_ok=True)
        
        with zipfile.ZipFile(zip_path, 'r') as zip_ref:
            zip_ref.extractall(extract_path)
        
        # Detect input type
        input_type, quantms_config = detect_input_type(extract_path)
        logger.info(f"Detected input type: {input_type}")
        
        if input_type == 'unknown':
            return jsonify({
                'error': 'Could not detect input type. Please ensure your ZIP file contains valid data files.'
            }), 400
        
        # Run PMultiQC
        result = run_pmultiqc(extract_path, output_dir, input_type, quantms_config)
        
        if not result['success']:
            return jsonify({
                'error': f'Failed to generate report: {result["message"]}'
            }), 500
        
        # Create zip report
        zip_report_path = os.path.join(output_dir, f'pmultiqc_report_{job_id}.zip')
        if not create_zip_report(output_dir, zip_report_path):
            logger.error(f"Failed to create zip report for job {job_id}")
            return jsonify({
                'error': 'Failed to create zip report. Please try again.'
            }), 500
        
        logger.info(f"Successfully completed job {job_id}")

        # Copy HTML report to public directory for online viewing
        try:
            html_urls = copy_html_report_for_online_viewing(output_dir, job_id)
            logger.info(f"Generated html_urls for sync job {job_id}: {html_urls}")
        except Exception as e:
            logger.error(f"Failed to copy HTML reports for sync job {job_id}: {e}")
            html_urls = None

        job_update = {
            'job_id': job_id,
            'status': 'completed',
            'input_type': input_type,
            'finished_at_seconds': current_time()[0],
            'finished_at': current_time()[1]
        }
        
        if html_urls:
            job_update['html_report_urls'] = html_urls

        job_status_dict[job_id] = job_update

        # Save job data to disk for future retrieval
        save_job_data_to_disk(job_id, job_update, output_dir)

        base_url = app.config['BASE_URL'].rstrip('/')
        response_data = {
            'job_id': job_id,
            'status': 'completed',
            'input_type': input_type,
            'download_url': f'{base_url}/download-report/{job_id}',
            'message': 'Report generated successfully',
            'finished_at': current_time()[1]
        }
        
        if html_urls:
            response_data['html_report_urls'] = html_urls

        return jsonify(response_data)
        
    except Exception as e:
        logger.error(f"Error in generate_report: {e}")
        return jsonify({'error': 'Internal server error'}), 500

@app.route('/upload-large', methods=['POST'])
def upload_large_file():
    """
    Upload large files using streaming to handle files > 1GB.
    
    Expected form data:
    - file: ZIP file containing the data to analyze
    
    Returns:
    - JSON response with job ID and status
    """
    try:
        # Check if file was uploaded
        if 'file' not in request.files:
            return jsonify({'error': 'No file provided'}), 400
        
        file = request.files['file']
        if file.filename == '':
            return jsonify({'error': 'No file selected'}), 400
        
        if not allowed_file(file.filename):
            return jsonify({'error': 'Only ZIP files are allowed'}), 400
        
        # Generate unique job ID
        job_id = str(uuid.uuid4())
        
        # Create job directories
        upload_dir = os.path.join(app.config['UPLOAD_FOLDER'], job_id)
        output_dir = os.path.join(app.config['OUTPUT_FOLDER'], job_id)
        os.makedirs(upload_dir, exist_ok=True)
        os.makedirs(output_dir, exist_ok=True)
        
        # Save uploaded file using streaming
        filename = secure_filename(file.filename)
        zip_path = os.path.join(upload_dir, filename)
        
        logger.info(f"Starting large file upload for job {job_id}: {filename}")
        
        # Stream the file to disk
        with open(zip_path, 'wb') as f:
            while True:
                chunk = file.stream.read(8192)  # 8KB chunks
                if not chunk:
                    break
                f.write(chunk)
        
        file_size = os.path.getsize(zip_path)
        logger.info(f"Upload completed for job {job_id}: {file_size} bytes")
        
        # Extract the zip file
        extract_path = os.path.join(upload_dir, 'extracted')
        os.makedirs(extract_path, exist_ok=True)
        
        logger.info(f"Extracting ZIP file for job {job_id}")
        with zipfile.ZipFile(zip_path, 'r') as zip_ref:
            zip_ref.extractall(extract_path)
        
        # Detect input type
        input_type, quantms_config = detect_input_type(extract_path)
        logger.info(f"Detected input type: {input_type}")
        
        if input_type == 'unknown':
            return jsonify({
                'error': 'Could not detect input type. Please ensure your ZIP file contains valid data files.'
            }), 400
        
        # Run PMultiQC
        result = run_pmultiqc(extract_path, output_dir, input_type, quantms_config)
        
        if not result['success']:
            return jsonify({
                'error': f'Failed to generate report: {result["message"]}'
            }), 500
        
        # Create zip report
        zip_report_path = os.path.join(output_dir, f'pmultiqc_report_{job_id}.zip')
        if not create_zip_report(output_dir, zip_report_path):
            logger.error(f"Failed to create zip report for job {job_id}")
            return jsonify({
                'error': 'Failed to create zip report. Please try again.'
            }), 500
        
        logger.info(f"Successfully completed large file job {job_id}")
        
        # Copy HTML report to public directory for online viewing
        try:
            html_urls = copy_html_report_for_online_viewing(output_dir, job_id)
            logger.info(f"Generated html_urls for large file job {job_id}: {html_urls}")
        except Exception as e:
            logger.error(f"Failed to copy HTML reports for large file job {job_id}: {e}")
            html_urls = None

        job_update = {
            'job_id': job_id,
            'status': 'completed',
            'input_type': input_type,
            'finished_at_seconds': current_time()[0],
            'finished_at': current_time()[1]
        }
        
        if html_urls:
            job_update['html_report_urls'] = html_urls

        job_status_dict[job_id] = job_update

        # Save job data to disk for future retrieval
        save_job_data_to_disk(job_id, job_update, output_dir)

        base_url = app.config['BASE_URL'].rstrip('/')
        response_data = {
            'job_id': job_id,
            'status': 'completed',
            'input_type': input_type,
            'download_url': f'{base_url}/download-report/{job_id}',
            'message': 'Report generated successfully',
            'file_size': file_size,
            'finished_at': current_time()[1]
        }
        
        if html_urls:
            response_data['html_report_urls'] = html_urls

        return jsonify(response_data)
        
    except Exception as e:
        logger.error(f"Error in upload_large_file: {e}")
        return jsonify({'error': 'Internal server error'}), 500

def run_pmultiqc_with_progress(input_path: str, output_path: str, input_type: str, pmultiqc_config: str, job_id: str) -> Dict[str, Any]:
    """
    Run PMultiQC with real-time progress updates stored in job_status_dict.
    """
    try:
        # Set up MultiQC arguments based on input type
        args = [
            'multiqc',
            input_path,
            '-o', output_path,
            '--force'
        ]
        
        # Add type-specific arguments
        if input_type == 'maxquant':
            args.extend(['--parse_maxquant'])
        elif input_type == "quantms":
            args.extend(["--config", pmultiqc_config])
        elif input_type == 'diann':
            # DIANN files are handled automatically by pmultiqc
            pass
        elif input_type == 'mzidentml':
            args.extend(['--mzid_plugin'])
        
        # Run MultiQC with PMultiQC plugin
        logger.info(f"Running PMultiQC with args: {args}")
        
        # Initialize console output in job status
        job_status_dict[job_id]['console_output'] = []
        job_status_dict[job_id]['console_errors'] = []
        
        # Run multiqc as a subprocess with real-time output
        try:
            # Run the multiqc command with real-time output
            process = subprocess.Popen(
                args, 
                stdout=subprocess.PIPE, 
                stderr=subprocess.PIPE, 
                text=True, 
                bufsize=1,
                universal_newlines=True
            )
            
            # Stream output in real-time
            output_lines = []
            error_lines = []
            
            while True:
                output = process.stdout.readline()
                error = process.stderr.readline()
                
                if output:
                    output_line = output.strip()
                    output_lines.append(output_line)
                    job_status_dict[job_id]['console_output'].append(output_line)
                    logger.info(f"MultiQC: {output_line}")
                
                if error:
                    error_line = error.strip()
                    error_lines.append(error_line)
                    job_status_dict[job_id]['console_errors'].append(error_line)
                    logger.warning(f"MultiQC: {error_line}")
                
                # Check if process has finished
                if process.poll() is not None:
                    break
            
            # Get any remaining output
            remaining_output, remaining_error = process.communicate()
            if remaining_output:
                for line in remaining_output.strip().split('\n'):
                    if line.strip():
                        output_lines.append(line.strip())
                        job_status_dict[job_id]['console_output'].append(line.strip())
                        logger.info(f"MultiQC remaining output: {line.strip()}")
            if remaining_error:
                for line in remaining_error.strip().split('\n'):
                    if line.strip():
                        error_lines.append(line.strip())
                        job_status_dict[job_id]['console_errors'].append(line.strip())
                        logger.warning(f"MultiQC remaining error: {line.strip()}")
            
            if process.returncode == 0:
                return {
                    'success': True, 
                    'message': 'Report generated successfully',
                    'output': output_lines,
                    'errors': error_lines
                }
            else:
                return {
                    'success': False, 
                    'message': f'MultiQC failed with return code {process.returncode}',
                    'output': output_lines,
                    'errors': error_lines
                }
                
        except subprocess.CalledProcessError as e:
            logger.error(f"MultiQC failed with return code {e.returncode}")
            logger.error(f"MultiQC stdout: {e.stdout}")
            logger.error(f"MultiQC stderr: {e.stderr}")
            return {
                'success': False, 
                'message': f'MultiQC failed: {e.stderr}',
                'output': e.stdout.split('\n') if e.stdout else [],
                'errors': e.stderr.split('\n') if e.stderr else []
            }
            
    except Exception as e:
        logger.error(f"Error running PMultiQC: {e}")
        return {'success': False, 'message': str(e), 'output': [], 'errors': []}

def process_job_async(job_id: str, extract_path: str, output_dir: str, input_type: str, quantms_config: str):
    """Process a job asynchronously in a separate thread."""
    try:
        job_status_dict[job_id] = {
            'status': 'running_pmultiqc',
            'progress': 0,
            'console_output': [],
            'console_errors': [],
            'processing_stage': 'Starting pmultiqc analysis...'
        }
        
        # Run PMultiQC with real-time progress
        logger.info(f"Starting async processing for job {job_id}")
        result = run_pmultiqc_with_progress(extract_path, output_dir, input_type, quantms_config, job_id)
        
        if not result['success']:
            job_status_dict[job_id].update({
                'status': 'failed',
                'error': result['message'],
                'finished_at_seconds': current_time()[0],
                'finished_at': current_time()[1]
            })
            return
        
        job_status_dict[job_id].update({
            'status': 'creating_report', 
            'progress': 75,
            'processing_stage': 'Creating zip report...'
        })
        
        # Create zip report
        zip_report_path = os.path.join(output_dir, f'pmultiqc_report_{job_id}.zip')
        if not create_zip_report(output_dir, zip_report_path):
            job_status_dict[job_id].update({
                'status': 'failed',
                'error': 'Failed to create zip report',
                'finished_at_seconds': current_time()[0],
                'finished_at': current_time()[1],
            })
            return
        
        # Copy HTML report to public directory for online viewing
        try:
            html_urls = copy_html_report_for_online_viewing(output_dir, job_id)
            logger.info(f"Generated html_urls for job {job_id}: {html_urls}")
        except Exception as e:
            logger.error(f"Failed to copy HTML reports for job {job_id}: {e}")
            html_urls = None
        
        job_update = {
            'status': 'completed',
            'progress': 100,
            'input_type': input_type,
            'processing_stage': 'Analysis completed successfully!',
            'finished_at_seconds': current_time()[0],
            'finished_at': current_time()[1]
        }
        
        # Only add html_report_urls if we successfully generated them
        if html_urls:
            job_update['html_report_urls'] = html_urls
        
        job_status_dict[job_id].update(job_update)
        
        # Save job data to disk for future retrieval
        save_job_data_to_disk(job_id, job_status_dict[job_id], output_dir)
        
        logger.info(f"Successfully completed async job {job_id}")
        
    except Exception as e:
        logger.error(f"Error in async processing for job {job_id}: {e}")
        job_status_dict[job_id].update({
            'status': 'failed',
            'error': str(e),
            'finished_at_seconds': current_time()[0],
            'finished_at': current_time()[1]
        })

@app.route('/debug-upload', methods=['POST'])
def debug_upload():
    """Debug endpoint to check what's being received."""
    logger.info("Debug upload endpoint called")
    logger.info(f"Request files: {list(request.files.keys())}")
    logger.info(f"Request form: {list(request.form.keys())}")
    logger.info(f"Content type: {request.content_type}")
    return jsonify({
        'files': list(request.files.keys()),
        'form': list(request.form.keys()),
        'content_type': request.content_type
    })

@app.route('/debug-reports', methods=['GET'])
def debug_reports():
    """Debug endpoint to check HTML reports directory."""
    try:
        reports_dir = app.config['HTML_REPORTS_FOLDER']
        debug_info = {
            'reports_directory': reports_dir,
            'directory_exists': os.path.exists(reports_dir),
            'contents': {}
        }
        
        if os.path.exists(reports_dir):
            for item in os.listdir(reports_dir):
                item_path = os.path.join(reports_dir, item)
                if os.path.isdir(item_path):
                    debug_info['contents'][item] = {}
                    for root, dirs, files in os.walk(item_path):
                        rel_path = os.path.relpath(root, item_path)
                        if rel_path == '.':
                            rel_path = 'root'
                        debug_info['contents'][item][rel_path] = {
                            'directories': dirs,
                            'files': files
                        }
        
        return jsonify(debug_info)
        
    except Exception as e:
        logger.error(f"Error in debug_reports: {e}")
        return jsonify({'error': str(e)}), 500

@app.route('/debug-job/<job_id>', methods=['GET'])
def debug_job(job_id):
    """Debug endpoint to check job status and data."""
    try:
        # Check current time and show all available jobs
        now = current_time()
        debug_info = {
            'current_time': now[1],
            'current_timestamp': now[0],
            'total_jobs_in_memory': len(job_status_dict),
            'all_job_ids': list(job_status_dict.keys()),
            'requested_job_id': job_id,
            'job_found_in_memory': job_id in job_status_dict
        }
        
        if job_id in job_status_dict:
            job_data = job_status_dict[job_id]
            
            # Calculate how long ago the job finished
            finished_at = job_data.get('finished_at_seconds')
            if finished_at:
                age_seconds = now[0] - finished_at
                debug_info['job_age_seconds'] = age_seconds
                debug_info['job_age_human'] = f"{age_seconds // 3600}h {(age_seconds % 3600) // 60}m {age_seconds % 60}s"
            
            # Also check if output directory exists
            output_dir = os.path.join(app.config['OUTPUT_FOLDER'], job_id)
            debug_info.update({
                'job_data': job_data,
                'output_directory': output_dir,
                'output_exists': os.path.exists(output_dir),
                'output_contents': []
            })
            
            if os.path.exists(output_dir):
                for root, dirs, files in os.walk(output_dir):
                    for file in files:
                        debug_info['output_contents'].append(os.path.join(root, file))
            
            return jsonify(debug_info)
        else:
            # Check if files exist even though job is not in memory
            output_dir = os.path.join(app.config['OUTPUT_FOLDER'], job_id)
            zip_report_path = os.path.join(output_dir, f'pmultiqc_report_{job_id}.zip')
            
            debug_info.update({
                'output_directory': output_dir,
                'output_exists': os.path.exists(output_dir),
                'zip_report_exists': os.path.exists(zip_report_path),
                'fallback_would_trigger': os.path.exists(zip_report_path)
            })
            
            return jsonify(debug_info)
            
    except Exception as e:
        logger.error(f"Error in debug_job: {e}")
        return jsonify({'error': str(e)}), 500

@app.route('/upload-async', methods=['POST'])
def upload_async():
    """
    Upload and process files asynchronously for very large files.
    
    Expected form data:
    - files: ZIP file containing the data to analyze
    
    Returns:
    - JSON response with job ID and status
    """
    try:
        # Debug: Log all form data and files
        logger.debug(f"Request files: {list(request.files.keys())}")
        logger.debug(f"Request form: {list(request.form.keys())}")
        
        # Check if file was uploaded
        if 'files' not in request.files:
            logger.error(f"No 'files' field found. Available fields: {list(request.files.keys())}")
            return jsonify({'error': 'No file provided'}), 400
        
        file = request.files['files']
        if file.filename == '':
            return jsonify({'error': 'No file selected'}), 400
        
        if not allowed_file(file.filename):
            return jsonify({'error': 'Only ZIP files are allowed'}), 400
        
        # Generate unique job ID
        job_id = str(uuid.uuid4())
        
        # Create job directories
        upload_dir = os.path.join(app.config['UPLOAD_FOLDER'], job_id)
        output_dir = os.path.join(app.config['OUTPUT_FOLDER'], job_id)
        os.makedirs(upload_dir, exist_ok=True)
        os.makedirs(output_dir, exist_ok=True)
        
        # Save uploaded file using streaming
        filename = secure_filename(file.filename)
        zip_path = os.path.join(upload_dir, filename)
        
        logger.info(f"Starting async upload for job {job_id}: {filename}")
        job_status_dict[job_id] = {'status': 'uploading', 'progress': 0}
        
        # Stream the file to disk
        with open(zip_path, 'wb') as f:
            while True:
                chunk = file.stream.read(8192)  # 8KB chunks
                if not chunk:
                    break
                f.write(chunk)
        
        file_size = os.path.getsize(zip_path)
        logger.info(f"Upload completed for job {job_id}: {file_size} bytes")
        
        # Extract the zip file
        extract_path = os.path.join(upload_dir, 'extracted')
        os.makedirs(extract_path, exist_ok=True)
        
        job_status_dict[job_id] = {'status': 'extracting', 'progress': 25}
        logger.info(f"Extracting ZIP file for job {job_id}")
        
        with zipfile.ZipFile(zip_path, 'r') as zip_ref:
            zip_ref.extractall(extract_path)
        
        # Detect input type
        input_type, quantms_config = detect_input_type(extract_path)
        logger.info(f"Detected input type: {input_type}")
        
        if input_type == 'unknown':
            job_status_dict[job_id] = {
                'status': 'failed',
                'error': 'Could not detect input type',
                'finished_at_seconds': current_time()[0],
                'finished_at': current_time()[1]
            }
            return jsonify({
                'error': 'Could not detect input type. Please ensure your ZIP file contains valid data files.'
            }), 400
        
        # Start async processing
        job_status_dict[job_id] = {'status': 'queued', 'progress': 50}
        thread = threading.Thread(
            target=process_job_async,
            args=(job_id, extract_path, output_dir, input_type, quantms_config)
        )
        thread.daemon = True
        thread.start()
        
        return jsonify({
            'job_id': job_id,
            'status': 'queued',
            'message': 'File uploaded successfully. Processing started.',
            'status_url': f'/job-status/{job_id}',
            'file_size': file_size
        })
        
    except Exception as e:
        logger.error(f"Error in upload_async: {e}")
        return jsonify({'error': 'Internal server error'}), 500

@app.route('/download-report/<job_id>', methods=['GET'])
def download_report(job_id):
    """
    Download the generated PMultiQC report.
    
    Args:
        job_id: The job ID returned from generate-report
    
    Returns:
        ZIP file containing the PMultiQC report
    """
    try:
        # Validate job_id format
        try:
            uuid.UUID(job_id)
        except ValueError:
            return jsonify({'error': 'Invalid job ID'}), 400
        
        output_dir = os.path.join(app.config['OUTPUT_FOLDER'], job_id)
        zip_report_path = os.path.join(output_dir, f'pmultiqc_report_{job_id}.zip')
        
        if not os.path.exists(zip_report_path):
            return jsonify({'error': 'Report not found'}), 404
        
        return send_file(
            zip_report_path,
            as_attachment=True,
            download_name=f'pmultiqc_report_{job_id}.zip',
            mimetype='application/zip'
        )
        
    except Exception as e:
        logger.error(f"Error in download_report: {e}")
        return jsonify({'error': 'Internal server error'}), 500

@app.route('/reports/<report_hash>/info', methods=['GET'])
def get_report_info(report_hash):
    """
    Get basic information about a report without loading large data files.
    """
    try:
        logger.info(f"Report info request: hash={report_hash}")
        
        # Validate report_hash format
        if not report_hash or len(report_hash) != 12 or not all(c in '0123456789abcdef' for c in report_hash):
            return jsonify({'error': 'Invalid report hash'}), 400
        
        report_dir = os.path.join(app.config['HTML_REPORTS_FOLDER'], report_hash)
        
        if not os.path.exists(report_dir):
            return jsonify({'error': 'Report not found'}), 404
        
        # Get basic info about the report
        info = {
            'report_hash': report_hash,
            'report_directory': report_dir,
            'subdirectories': [],
            'total_files': 0,
            'total_size_bytes': 0,
            'large_files': []
        }
        
        for root, dirs, files in os.walk(report_dir):
            rel_path = os.path.relpath(root, report_dir)
            if rel_path != '.':
                info['subdirectories'].append(rel_path)
            
            for file in files:
                if file == '.job_id':
                    continue
                    
                file_path = os.path.join(root, file)
                file_size = os.path.getsize(file_path)
                info['total_files'] += 1
                info['total_size_bytes'] += file_size
                
                # Track large files (>1MB)
                if file_size > 1024 * 1024:
                    rel_file_path = os.path.relpath(file_path, report_dir)
                    info['large_files'].append({
                        'path': rel_file_path,
                        'size_bytes': file_size,
                        'size_human': f"{file_size / (1024*1024):.1f} MB"
                    })
        
        # Add human readable total size
        if info['total_size_bytes'] > 1024 * 1024:
            info['total_size_human'] = f"{info['total_size_bytes'] / (1024*1024):.1f} MB"
        else:
            info['total_size_human'] = f"{info['total_size_bytes'] / 1024:.1f} KB"
        
        return jsonify(info)
        
    except Exception as e:
        logger.error(f"Error getting report info: {e}")
        return jsonify({'error': 'Internal server error'}), 500

@app.route('/reports/<report_hash>/test', methods=['GET'])
def test_report_access(report_hash):
    """
    Simple test page to verify report directory access without loading large files.
    """
    try:
        # Validate report_hash format
        if not report_hash or len(report_hash) != 12 or not all(c in '0123456789abcdef' for c in report_hash):
            return "Invalid report hash", 400
        
        report_dir = os.path.join(app.config['HTML_REPORTS_FOLDER'], report_hash)
        
        if not os.path.exists(report_dir):
            return "Report not found", 404
        
        # Generate a simple test HTML page
        html_content = f"""
        <!DOCTYPE html>
        <html>
        <head>
            <title>Report Test - {report_hash}</title>
            <style>
                body {{ font-family: Arial, sans-serif; margin: 20px; }}
                .info {{ background: #e7f3ff; padding: 10px; border-radius: 5px; margin: 10px 0; }}
                .file-list {{ background: #f5f5f5; padding: 10px; border-radius: 5px; max-height: 300px; overflow-y: auto; }}
            </style>
        </head>
        <body>
            <h1>Report Test Page</h1>
            <div class="info">
                <strong>Report Hash:</strong> {report_hash}<br>
                <strong>Directory:</strong> {report_dir}<br>
                <strong>Status:</strong> Directory exists and is accessible
            </div>
            
            <h2>Available Files:</h2>
            <div class="file-list">
        """
        
        # List files in the report directory
        for root, dirs, files in os.walk(report_dir):
            rel_path = os.path.relpath(root, report_dir)
            if rel_path != '.':
                html_content += f"<strong>{rel_path}/</strong><br>"
            
            for file in files:
                if file == '.job_id':
                    continue
                file_path = os.path.join(root, file)
                file_size = os.path.getsize(file_path)
                rel_file_path = os.path.relpath(file_path, report_dir)
                
                if file == 'multiqc_report.html':
                    html_content += f'&nbsp;&nbsp;<a href="/reports/{report_hash}/{rel_file_path}" target="_blank"><strong>{file}</strong></a> ({file_size} bytes)<br>'
                else:
                    html_content += f"&nbsp;&nbsp;{file} ({file_size} bytes)<br>"
        
        html_content += """
            </div>
            <div class="info">
                <strong>Note:</strong> This is a test page to verify report access. 
                Click on the multiqc_report.html links above to view the actual reports.
            </div>
        </body>
        </html>
        """
        
        return html_content, 200, {'Content-Type': 'text/html'}
        
    except Exception as e:
        logger.error(f"Error in test_report_access: {e}")
        return f"Error: {str(e)}", 500

@app.route('/reports/<report_hash>/<path:filename>', methods=['GET'])
def serve_html_report(report_hash, filename):
    """
    Serve HTML reports for online viewing.
    
    Args:
        report_hash: Hash identifier for the report
        filename: Filename to serve (e.g., multiqc_report.html)
    
    Returns:
        HTML file or other assets from the MultiQC report
    """
    try:
        start_time = time.time()
        logger.info(f"HTML report request: hash={report_hash}, filename={filename}")
        
        # Validate report_hash format (12 character hex)
        if not report_hash or len(report_hash) != 12 or not all(c in '0123456789abcdef' for c in report_hash):
            logger.error(f"Invalid report hash format: {report_hash}")
            return jsonify({'error': 'Invalid report hash'}), 400
        
        # Construct file path
        file_path = os.path.join(app.config['HTML_REPORTS_FOLDER'], report_hash, filename)
        logger.info(f"Looking for file at: {file_path}")
        
        # List contents of the report directory for debugging
        report_dir = os.path.join(app.config['HTML_REPORTS_FOLDER'], report_hash)
        if os.path.exists(report_dir):
            logger.info(f"Report directory exists. Contents:")
            for root, dirs, files in os.walk(report_dir):
                level = root.replace(report_dir, '').count(os.sep)
                indent = ' ' * 2 * level
                logger.info(f"{indent}{os.path.basename(root)}/")
                subindent = ' ' * 2 * (level + 1)
                for file in files:
                    logger.info(f"{subindent}{file}")
        else:
            logger.error(f"Report directory does not exist: {report_dir}")
            return jsonify({'error': 'Report directory not found'}), 404
        
        # Security check: ensure the file is within the HTML reports directory
        real_path = os.path.realpath(file_path)
        html_reports_real = os.path.realpath(app.config['HTML_REPORTS_FOLDER'])
        
        if not real_path.startswith(html_reports_real):
            logger.error(f"Security violation: {real_path} not within {html_reports_real}")
            return jsonify({'error': 'Access denied'}), 403
        
        if not os.path.exists(file_path):
            logger.error(f"File not found: {file_path}")
            return jsonify({'error': f'File not found: {filename}'}), 404
        
        # Determine MIME type based on file extension
        mime_type = 'text/html'  # Default for HTML files
        if filename.endswith('.css'):
            mime_type = 'text/css'
        elif filename.endswith('.js'):
            mime_type = 'application/javascript'
        elif filename.endswith('.png'):
            mime_type = 'image/png'
        elif filename.endswith('.jpg') or filename.endswith('.jpeg'):
            mime_type = 'image/jpeg'
        elif filename.endswith('.svg'):
            mime_type = 'image/svg+xml'
        elif filename.endswith('.ico'):
            mime_type = 'image/x-icon'
        elif filename.endswith('.json'):
            mime_type = 'application/json'
        elif filename.endswith('.txt'):
            mime_type = 'text/plain'
        elif filename.endswith('.parquet'):
            mime_type = 'application/octet-stream'
        elif filename.endswith('.log'):
            mime_type = 'text/plain'
        
        # Get file size for logging
        file_size = os.path.getsize(file_path)
        logger.info(f"Serving file: {file_path} ({file_size} bytes) with MIME type: {mime_type}")
        
        # Log timing for large files
        if file_size > 1024 * 1024:  # Files larger than 1MB
            logger.warning(f"Large file being served: {filename} ({file_size / (1024*1024):.1f} MB)")
        
        # Add cache headers for static assets
        if filename.endswith(('.css', '.js', '.png', '.jpg', '.jpeg', '.svg', '.ico')):
            response = send_file(file_path, mimetype=mime_type, 
                               max_age=3600,  # Cache for 1 hour
                               as_attachment=False)
        else:
            response = send_file(file_path, mimetype=mime_type, as_attachment=False)
        
        # Log completion time
        elapsed = time.time() - start_time
        logger.info(f"File served in {elapsed:.2f}s: {filename}")
        
        return response
        
    except Exception as e:
        logger.error(f"Error serving HTML report: {e}")
        return jsonify({'error': 'Internal server error'}), 500

@app.route('/pride/<report_hash>/<path:filename>', methods=['GET'])
def serve_pride_html_report_redirect(report_hash, filename):
    """
    Backward compatibility: redirect old PRIDE URLs to new reports URLs.
    """
    new_url = f'/reports/{report_hash}/{filename}'
    logger.info(f"Redirecting PRIDE URL /pride/{report_hash}/{filename} to {new_url}")
    return redirect(new_url, code=301)

@app.route('/job-status/<job_id>', methods=['GET'])
def job_status(job_id):
    """
    Get the status of a job.
    
    Args:
        job_id: The job ID to check
    
    Returns:
        JSON with job status information
    """
    try:
        # Validate job_id format
        try:
            uuid.UUID(job_id)
        except ValueError:
            return jsonify({'error': 'Invalid job ID'}), 400
        
        # Check if job exists in tracking dictionary
        logger.info(f"Checking job {job_id}. job_status_dict contains {len(job_status_dict)} jobs: {list(job_status_dict.keys())}")
        if job_id in job_status_dict:
            job_data = job_status_dict[job_id]
            logger.info(f"Job {job_id} found in job_status_dict with status: {job_data.get('status')} and keys: {list(job_data.keys())}")
            return jsonify(job_data)
        else:
            logger.warning(f"Job {job_id} NOT found in job_status_dict")
        
        # Check if job is completed (for backward compatibility)
        output_dir = os.path.join(app.config['OUTPUT_FOLDER'], job_id)
        zip_report_path = os.path.join(output_dir, f'pmultiqc_report_{job_id}.zip')
        
        if os.path.exists(zip_report_path):
            logger.info(f"Job {job_id} not in job_status_dict but zip file exists")
            
            # Try to find HTML reports for this job
            html_reports_found = {}
            reports_dir = app.config['HTML_REPORTS_FOLDER']
            
            # Check all report directories for this job
            if os.path.exists(reports_dir):
                for report_hash in os.listdir(reports_dir):
                    report_path = os.path.join(reports_dir, report_hash)
                    job_file_path = os.path.join(report_path, '.job_id')
                    
                    if os.path.exists(job_file_path):
                        try:
                            with open(job_file_path, 'r') as f:
                                stored_job_id = f.read().strip()
                                
                            if stored_job_id == job_id:
                                # Found reports for this job - scan for HTML files
                                for root, dirs, files in os.walk(report_path):
                                    for file in files:
                                        if file == 'multiqc_report.html':
                                            rel_path = os.path.relpath(root, report_path)
                                            if rel_path == '.':
                                                rel_path = 'root'
                                            html_rel_path = os.path.relpath(os.path.join(root, file), report_path)
                                            # Use base URL to generate absolute URLs
                                            base_url = app.config['BASE_URL'].rstrip('/')
                                            html_url = f'{base_url}/reports/{report_hash}/{html_rel_path}'
                                            html_reports_found[rel_path] = html_url
                                            logger.info(f"Found HTML report for job {job_id}: {html_url}")
                                break
                        except Exception as e:
                            logger.error(f"Error reading job file {job_file_path}: {e}")
            
            # Try to recover additional job data from output directory
            console_output = []
            console_errors = []
            results = []
            input_type = None
            accession = None
            
            # Look for console output files
            console_output_file = os.path.join(output_dir, 'console_output.txt')
            console_errors_file = os.path.join(output_dir, 'console_errors.txt')
            job_info_file = os.path.join(output_dir, 'job_info.json')
            
            if os.path.exists(console_output_file):
                try:
                    with open(console_output_file, 'r', encoding='utf-8') as f:
                        console_output = [line.strip() for line in f.readlines() if line.strip()]
                    logger.info(f"Recovered {len(console_output)} console output lines for job {job_id}")
                except Exception as e:
                    logger.error(f"Error reading console output file: {e}")
            
            if os.path.exists(console_errors_file):
                try:
                    with open(console_errors_file, 'r', encoding='utf-8') as f:
                        console_errors = [line.strip() for line in f.readlines() if line.strip()]
                    logger.info(f"Recovered {len(console_errors)} console error lines for job {job_id}")
                except Exception as e:
                    logger.error(f"Error reading console errors file: {e}")
            
            if os.path.exists(job_info_file):
                try:
                    with open(job_info_file, 'r', encoding='utf-8') as f:
                        job_info = json.load(f)
                        input_type = job_info.get('input_type')
                        accession = job_info.get('accession')
                        results = job_info.get('results', [])
                    logger.info(f"Recovered job info for job {job_id}: input_type={input_type}, accession={accession}, results={len(results)}")
                except Exception as e:
                    logger.error(f"Error reading job info file: {e}")
            
            # Fallback response with recovered data
            base_url = app.config['BASE_URL'].rstrip('/')
            response = {
                'job_id': job_id,
                'status': 'completed',
                'progress': 100,
                'download_url': f'{base_url}/download-report/{job_id}',
                'console_output': console_output,
                'console_errors': console_errors,
                'results': results
            }
            
            if input_type:
                response['input_type'] = input_type
            if accession:
                response['accession'] = accession
            
            if html_reports_found:
                response['html_report_urls'] = html_reports_found
                logger.info(f"Added HTML reports to fallback response: {html_reports_found}")
            
            return jsonify(response)
        else:
            return jsonify({
                'job_id': job_id,
                'status': 'not_found'
            }), 404
        
    except Exception as e:
        logger.error(f"Error in job_status: {e}")
        return jsonify({'error': 'Internal server error'}), 500
    
@app.route('/cleanup/<job_id>', methods=['DELETE'])
def cleanup_job(job_id):
    """
    Clean up job files.
    
    Args:
        job_id: The job ID to clean up
    
    Returns:
        JSON response indicating cleanup status
    """
    try:
        # Validate job_id format
        try:
            uuid.UUID(job_id)
        except ValueError:
            return jsonify({'error': 'Invalid job ID'}), 400
        
        upload_dir = os.path.join(app.config['UPLOAD_FOLDER'], job_id)
        output_dir = os.path.join(app.config['OUTPUT_FOLDER'], job_id)
        
        # Remove directories if they exist
        if os.path.exists(upload_dir):
            shutil.rmtree(upload_dir)
        
        if os.path.exists(output_dir):
            shutil.rmtree(output_dir)
        
        return jsonify({
            'job_id': job_id,
            'status': 'cleaned_up',
            'message': 'Job files removed successfully'
        })
        
    except Exception as e:
        logger.error(f"Error in cleanup_job: {e}")
        return jsonify({'error': 'Internal server error'}), 500

def process_pride_job_async(job_id: str, accession: str, output_dir: str):
    """Process a PRIDE job asynchronously in a separate thread."""
    try:
        job_status_dict[job_id] = {
            'status': 'fetching_files',
            'progress': 10,
        }
        
        # Fetch files from PRIDE
        logger.info(f"Starting PRIDE processing for job {job_id}, accession {accession}")
        files = get_pride_files(accession)
        
        job_status_dict[job_id] = {'status': 'filtering_files', 'progress': 20}
        
        # Filter search engine output files
        search_files = filter_search_files(files)
        
        if not search_files:
            job_status_dict[job_id] = {
                'status': 'failed',
                'error': 'No search engine output files found in PRIDE submission',
                'finished_at_seconds': current_time()[0],
                'finished_at': current_time()[1]
            }
            return
        
        job_status_dict[job_id] = {'status': 'downloading_files', 'progress': 30}
        
        # Create download directory
        download_dir = os.path.join(app.config['UPLOAD_FOLDER'], job_id, 'pride_downloads')
        os.makedirs(download_dir, exist_ok=True)
        
        # Download files
        downloaded_files = []
        total_files = len(search_files)
        
        for i, file_info in enumerate(search_files):
            try:
                file_path = download_pride_file(file_info, download_dir)
                downloaded_files.append(file_path)
                progress = 30 + int((i + 1) / total_files * 40)  # 30-70% for downloads
                job_status_dict[job_id] = {'status': 'downloading_files', 'progress': progress}
            except Exception as e:
                logger.warning(f"Failed to download {file_info.get('fileName')}: {e}")
                # Continue with other files
        
        if not downloaded_files:
            job_status_dict[job_id] = {
                'status': 'failed',
                'error': 'Failed to download any files from PRIDE',
                'finished_at_seconds': current_time()[0],
                'finished_at': current_time()[1]
            }
            return
        
        job_status_dict[job_id] = {'status': 'extracting_files', 'progress': 70}
        
        # Process each downloaded file separately
        all_results = []
        total_processed = 0
        
        for i, zip_file in enumerate(downloaded_files):
            if zip_file.lower().endswith('.zip'):
                try:
                    logger.info(f"Processing {zip_file}")
                    
                    # Create separate directory for this file
                    file_name = os.path.splitext(os.path.basename(zip_file))[0]
                    file_extract_dir = os.path.join(download_dir, f'extracted_{file_name}')
                    file_output_dir = os.path.join(output_dir, f'report_{file_name}')
                    
                    os.makedirs(file_extract_dir, exist_ok=True)
                    os.makedirs(file_output_dir, exist_ok=True)
                    
                    # Extract this ZIP file
                    logger.info(f"Extracting {zip_file}")
                    with zipfile.ZipFile(zip_file, 'r') as zip_ref:
                        zip_ref.extractall(file_extract_dir)
                    
                    # Detect input type for this file
                    input_type, quantms_config = detect_input_type(file_extract_dir)
                    logger.info(f"Detected input type for {file_name}: {input_type}")
                    
                    if input_type == 'unknown':
                        logger.warning(f"Could not detect input type for {file_name}")
                        continue
                    
                    # Run PMultiQC on this file
                    result = run_pmultiqc(file_extract_dir, file_output_dir, input_type, quantms_config)
                    
                    if result['success']:
                        # Create zip report for this file
                        zip_report_path = os.path.join(file_output_dir, f'pmultiqc_report_{file_name}.zip')
                        if create_zip_report(file_output_dir, zip_report_path):
                            all_results.append({
                                'file_name': file_name,
                                'input_type': input_type,
                                'report_path': zip_report_path,
                                'output': result.get('output', []),
                                'errors': result.get('errors', [])
                            })
                            total_processed += 1
                            logger.info(f"Successfully processed {file_name}")
                        else:
                            logger.error(f"Failed to create zip report for {file_name}")
                    else:
                        logger.error(f"Failed to process {file_name}: {result['message']}")
                        
                except Exception as e:
                    logger.error(f"Error processing {zip_file}: {e}")
        
        if not all_results:
            job_status_dict[job_id] = {
                'status': 'failed',
                'error': 'Failed to process any files from PRIDE',
                'finished_at_seconds': current_time()[0],
                'finished_at': current_time()[1]
            }
            return
        
        # Create combined zip report with all results
        combined_zip_path = os.path.join(output_dir, f'pmultiqc_report_{job_id}.zip')
        if create_zip_report(output_dir, combined_zip_path):
            # Copy HTML report to public directory for PRIDE datasets
            try:
                html_urls = copy_html_report_for_online_viewing(output_dir, job_id, accession)
                logger.info(f"Generated html_urls for PRIDE job {job_id}: {html_urls}")
            except Exception as e:
                logger.error(f"Failed to copy HTML reports for PRIDE job {job_id}: {e}")
                html_urls = None
            
            # Combine all console outputs
            all_output = []
            all_errors = []
            for result in all_results:
                all_output.extend(result.get('output', []))
                all_errors.extend(result.get('errors', []))
            
            job_update = {
                'status': 'completed',
                'progress': 100,
                'input_type': 'multiple',
                'accession': accession,
                'files_processed': total_processed,
                'total_files': len(downloaded_files),
                'results': all_results,
                'console_output': all_output,
                'console_errors': all_errors,
                'finished_at_seconds': current_time()[0],
                'finished_at': current_time()[1]
            }
            
            # Only add html_report_urls if we successfully generated them
            if html_urls:
                job_update['html_report_urls'] = html_urls
            
            # Store job data with extra verification
            job_status_dict[job_id] = job_update
            
            # Save job data to disk for future retrieval
            save_job_data_to_disk(job_id, job_update, output_dir)
            
            # Verify the job was actually stored
            if job_id in job_status_dict:
                stored_data = job_status_dict[job_id]
                logger.info(f"Successfully completed PRIDE job {job_id} with {total_processed} files processed")
                logger.info(f"PRIDE job {job_id} stored in job_status_dict with keys: {list(job_update.keys())}")
                logger.info(f"Verification: job is in dict with status: {stored_data.get('status')} and html_report_urls: {bool(stored_data.get('html_report_urls'))}")
                logger.info(f"job_status_dict now contains {len(job_status_dict)} jobs: {list(job_status_dict.keys())}")
                
                # Force a small delay to ensure any race conditions are avoided
                import time
                time.sleep(0.1)
                
                # Double-check after delay
                if job_id in job_status_dict:
                    logger.info(f"Double-check: PRIDE job {job_id} still in job_status_dict after delay")
                else:
                    logger.error(f"ERROR: PRIDE job {job_id} disappeared from job_status_dict after delay!")
            else:
                logger.error(f"ERROR: Failed to store PRIDE job {job_id} in job_status_dict!")
        else:
            job_status_dict[job_id] = {
                'status': 'failed',
                'error': 'Failed to create combined zip report',
                'finished_at_seconds': current_time()[0],
                'finished_at': current_time()[1]
            }
        
    except Exception as e:
        logger.error(f"Error in PRIDE processing for job {job_id}: {e}")
        job_status_dict[job_id] = {
            'status': 'failed',
            'error': str(e),
            'finished_at_seconds': current_time()[0],
            'finished_at': current_time()[1]
        }

@app.route('/process-pride', methods=['POST'])
def process_pride():
    """Process PRIDE dataset asynchronously."""
    try:
        data = request.get_json()
        if not data or 'accession' not in data:
            return jsonify({'error': 'Missing accession parameter'}), 400
        
        accession = data['accession'].strip().upper()
        
        # Validate accession format
        if not accession.startswith('PXD'):
            return jsonify({'error': 'Invalid PRIDE accession format. Should start with PXD'}), 400
        
        # Generate unique job ID
        job_id = str(uuid.uuid4())
        
        # Create output directory
        output_dir = os.path.join(app.config['OUTPUT_FOLDER'], job_id)
        os.makedirs(output_dir, exist_ok=True)
        
        # Initialize job status
        job_status_dict[job_id] = {
            'status': 'queued',
            'accession': accession,
            'progress': 0,
            'started_at': current_time()[1],
            'finished_at': None,
            'error': None,
            'input_type': None,
            'files_processed': 0,
            'total_files': 0,
            'files_downloaded': 0,
            'console_output': [],
            'console_errors': [],
            'results': [],
            'html_report_urls': None
        }
        
        # Start async processing
        thread = threading.Thread(
            target=process_pride_job_async,
            args=(job_id, accession, output_dir)
        )
        thread.daemon = True
        thread.start()
        
        logger.info(f"Started PRIDE processing for accession {accession} with job ID {job_id}")
        
        return jsonify({
            'job_id': job_id,
            'accession': accession,
            'status': 'queued',
            'message': 'PRIDE processing started successfully',
            'status_url': f'/job-status/{job_id}'
        }), 200
        
    except Exception as e:
        logger.error(f"Error starting PRIDE processing: {e}")
        return jsonify({'error': f'Failed to start PRIDE processing: {str(e)}'}), 500

@app.route('/submit', methods=['GET'])
def submit_pride():
    """Direct submission page for PRIDE datasets via URL parameter."""
    accession = request.args.get('accession', '').strip().upper()
    
    if not accession:
        # No accession provided, show the main page
        return render_template('index.html')
    
    # Validate accession format
    if not accession.startswith('PXD'):
        return render_template('index.html', error=f'Invalid PRIDE accession format: {accession}. Should start with PXD.')
    
    # Render the submission page with the accession pre-filled
    return render_template('submit.html', accession=accession)

@app.route('/api/submit-pride', methods=['POST'])
def api_submit_pride():
    """API endpoint for direct PRIDE submission."""
    try:
        data = request.get_json()
        if not data or 'accession' not in data:
            return jsonify({'error': 'Missing accession parameter'}), 400
        
        accession = data['accession'].strip().upper()
        
        # Validate accession format
        if not accession.startswith('PXD'):
            return jsonify({'error': 'Invalid PRIDE accession format. Should start with PXD'}), 400
        
        # Generate unique job ID
        job_id = str(uuid.uuid4())
        
        # Create output directory
        output_dir = os.path.join(app.config['OUTPUT_FOLDER'], job_id)
        os.makedirs(output_dir, exist_ok=True)
        
        # Initialize job status
        job_status_dict[job_id] = {
            'status': 'queued',
            'accession': accession,
            'progress': 0,
            'started_at': current_time()[1],
            'finished_at': None,
            'error': None,
            'input_type': None,
            'files_processed': 0,
            'total_files': 0,
            'files_downloaded': 0,
            'console_output': [],
            'console_errors': [],
            'results': [],
            'html_report_urls': None
        }
        
        # Start async processing
        thread = threading.Thread(
            target=process_pride_job_async,
            args=(job_id, accession, output_dir)
        )
        thread.daemon = True
        thread.start()
        
        logger.info(f"Started PRIDE processing for accession {accession} with job ID {job_id}")
        
        return jsonify({
            'job_id': job_id,
            'accession': accession,
            'status': 'queued',
            'message': 'PRIDE processing started successfully',
            'status_url': f'/job-status/{job_id}',
            'redirect_url': f'/results/{job_id}'
        }), 200
        
    except Exception as e:
        logger.error(f"Error starting PRIDE processing: {e}")
        return jsonify({'error': f'Failed to start PRIDE processing: {str(e)}'}), 500

@app.route('/results/<job_id>', methods=['GET'])
def view_results(job_id):
    """View results page for a specific job."""
    if job_id not in job_status_dict:
        return render_template('index.html', error=f'Job {job_id} not found.')
    
    job_data = job_status_dict[job_id]
    return render_template('results.html', job_id=job_id, job_data=job_data)

if __name__ == '__main__':
    app.run(debug=True, host='0.0.0.0', port=5000) 