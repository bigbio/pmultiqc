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
from datetime import datetime
from typing import Dict, Any, List, Optional
from urllib.parse import urlparse

from flask import Flask, request, jsonify, send_file, render_template
from flask_cors import CORS
from werkzeug.utils import secure_filename
import time

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

app = Flask(__name__, template_folder='templates')
CORS(app)

# Configuration
app.config['MAX_CONTENT_LENGTH'] = 2 * 1024 * 1024 * 1024  # 2GB max file size
app.config['UPLOAD_FOLDER'] = '/tmp/pmultiqc_uploads'
app.config['OUTPUT_FOLDER'] = '/tmp/pmultiqc_outputs'

# Ensure directories exist
os.makedirs(app.config['UPLOAD_FOLDER'], exist_ok=True)
os.makedirs(app.config['OUTPUT_FOLDER'], exist_ok=True)

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
        
        # Run multiqc as a subprocess
        import subprocess
        
        try:
            # Run the multiqc command
            result = subprocess.run(args, capture_output=True, text=True, check=True)
            logger.info(f"MultiQC stdout: {result.stdout}")
            if result.stderr:
                logger.warning(f"MultiQC stderr: {result.stderr}")
            return {'success': True, 'message': 'Report generated successfully'}
        except subprocess.CalledProcessError as e:
            logger.error(f"MultiQC failed with return code {e.returncode}")
            logger.error(f"MultiQC stdout: {e.stdout}")
            logger.error(f"MultiQC stderr: {e.stderr}")
            return {'success': False, 'message': f'MultiQC failed: {e.stderr}'}
            
    except Exception as e:
        logger.error(f"Error running PMultiQC: {e}")
        return {'success': False, 'message': str(e)}

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
        job_status_dict.pop(job_id, None)
        logger.info(f"Job ID {job_id} has been deleted at: {now}")

    threading.Timer(CLEANUP_JOB_SECONDS, cleanup_expired_jobs).start()

def cleanup_job_files(job_id, now_time):
    try:
        upload_dir = os.path.join(app.config['UPLOAD_FOLDER'], job_id)
        output_dir = os.path.join(app.config['OUTPUT_FOLDER'], job_id)

        for directory in [upload_dir, output_dir]:
            if os.path.exists(directory):
                shutil.rmtree(directory)
                logger.info(f"Removed directory {directory} for job {job_id} at {now_time}")

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

        job_status_dict[job_id] = {
            'job_id': job_id,
            'status': 'completed',
            'input_type': input_type,
            'finished_at_seconds': current_time()[0],
            'finished_at': current_time()[1]
        }

        return jsonify({
            'job_id': job_id,
            'status': 'completed',
            'input_type': input_type,
            'download_url': f'/download-report/{job_id}',
            'message': 'Report generated successfully',
            'finished_at': current_time()[1]
        })
        
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
        job_status_dict[job_id] = {
            'job_id': job_id,
            'status': 'completed',
            'input_type': input_type,
            'finished_at_seconds': current_time()[0],
            'finished_at': current_time()[1]
        }

        return jsonify({
            'job_id': job_id,
            'status': 'completed',
            'input_type': input_type,
            'download_url': f'/download-report/{job_id}',
            'message': 'Report generated successfully',
            'file_size': file_size,
            'finished_at': current_time()[1]
        })
        
    except Exception as e:
        logger.error(f"Error in upload_large_file: {e}")
        return jsonify({'error': 'Internal server error'}), 500

def process_job_async(job_id: str, extract_path: str, output_dir: str, input_type: str, quantms_config: str):
    """Process a job asynchronously in a separate thread."""
    try:
        job_status_dict[job_id] = {
            'status': 'processing',
            'progress': 0,
        }
        
        # Run PMultiQC
        logger.info(f"Starting async processing for job {job_id}")
        result = run_pmultiqc(extract_path, output_dir, input_type, quantms_config)
        
        if not result['success']:
            job_status_dict[job_id] = {
                'status': 'failed',
                'error': result['message'],
                'finished_at_seconds': current_time()[0],
                'finished_at': current_time()[1]
            }
            return
        
        job_status_dict[job_id] = {'status': 'creating_report', 'progress': 75}
        
        # Create zip report
        zip_report_path = os.path.join(output_dir, f'pmultiqc_report_{job_id}.zip')
        if not create_zip_report(output_dir, zip_report_path):
            job_status_dict[job_id] = {
                'status': 'failed',
                'error': 'Failed to create zip report',
                'finished_at_seconds': current_time()[0],
                'finished_at': current_time()[1],
            }
            return
        
        job_status_dict[job_id] = {
            'status': 'completed',
            'progress': 100,
            'input_type': input_type,
            'finished_at_seconds': current_time()[0],
            'finished_at': current_time()[1]
        }
        logger.info(f"Successfully completed async job {job_id}")
        
    except Exception as e:
        logger.error(f"Error in async processing for job {job_id}: {e}")
        job_status_dict[job_id] = {
            'status': 'failed',
            'error': str(e),
            'finished_at_seconds': current_time()[0],
            'finished_at': current_time()[1]
        }

@app.route('/upload-async', methods=['POST'])
def upload_async():
    """
    Upload and process files asynchronously for very large files.
    
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
        if job_id in job_status_dict:
            return jsonify(job_status_dict[job_id])
        
        # Check if job is completed (for backward compatibility)
        output_dir = os.path.join(app.config['OUTPUT_FOLDER'], job_id)
        zip_report_path = os.path.join(output_dir, f'pmultiqc_report_{job_id}.zip')
        
        if os.path.exists(zip_report_path):
            return jsonify({
                'job_id': job_id,
                'status': 'completed',
                'progress': 100,
                'download_url': f'/download-report/{job_id}'
            })
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
        
        job_status_dict[job_id] = {'status': 'detecting_type', 'progress': 75}
        
        # Detect input type
        input_type = detect_pride_input_type(search_files)
        logger.info(f"Detected PRIDE input type: {input_type}")
        
        if input_type == 'unknown':
            job_status_dict[job_id] = {
                'status': 'failed',
                'error': 'Could not detect input type from PRIDE files',
                'finished_at_seconds': current_time()[0],
                'finished_at': current_time()[1]
            }
            return
        
        job_status_dict[job_id] = {'status': 'processing', 'progress': 80}
        
        # Run PMultiQC
        result = run_pmultiqc(download_dir, output_dir, input_type, None)
        
        if not result['success']:
            job_status_dict[job_id] = {
                'status': 'failed',
                'error': result['message'],
                'finished_at_seconds': current_time()[0],
                'finished_at': current_time()[1]
            }
            return
        
        job_status_dict[job_id] = {'status': 'creating_report', 'progress': 90}
        
        # Create zip report
        zip_report_path = os.path.join(output_dir, f'pmultiqc_report_{job_id}.zip')
        if not create_zip_report(output_dir, zip_report_path):
            job_status_dict[job_id] = {
                'status': 'failed',
                'error': 'Failed to create zip report',
                'finished_at_seconds': current_time()[0],
                'finished_at': current_time()[1]
            }
            return
        
        job_status_dict[job_id] = {
            'status': 'completed',
            'progress': 100,
            'input_type': input_type,
            'accession': accession,
            'files_downloaded': len(downloaded_files),
            'finished_at_seconds': current_time()[0],
            'finished_at': current_time()[1]
        }
        logger.info(f"Successfully completed PRIDE job {job_id}")
        
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
    """
    Process data from PRIDE using an accession number.
    
    Expected JSON data:
    - accession: PRIDE accession (e.g., PXD039077)
    
    Returns:
    - JSON response with job ID and status
    """
    try:
        data = request.get_json()
        if not data or 'accession' not in data:
            return jsonify({'error': 'Accession number is required'}), 400
        
        accession = data['accession'].strip().upper()
        
        # Validate accession format (basic check)
        if not accession.startswith('PXD'):
            return jsonify({'error': 'Invalid PRIDE accession format. Should start with PXD'}), 400
        
        # Generate unique job ID
        job_id = str(uuid.uuid4())
        
        # Create job directories
        upload_dir = os.path.join(app.config['UPLOAD_FOLDER'], job_id)
        output_dir = os.path.join(app.config['OUTPUT_FOLDER'], job_id)
        os.makedirs(upload_dir, exist_ok=True)
        os.makedirs(output_dir, exist_ok=True)
        
        logger.info(f"Starting PRIDE processing for accession {accession}, job {job_id}")
        
        # Start async processing
        job_status_dict[job_id] = {'status': 'queued', 'progress': 0}
        thread = threading.Thread(
            target=process_pride_job_async,
            args=(job_id, accession, output_dir)
        )
        thread.daemon = True
        thread.start()
        
        return jsonify({
            'job_id': job_id,
            'accession': accession,
            'status': 'queued',
            'message': 'PRIDE processing started successfully',
            'status_url': f'/job-status/{job_id}'
        })
        
    except Exception as e:
        logger.error(f"Error in process_pride: {e}")
        return jsonify({'error': 'Internal server error'}), 500

if __name__ == '__main__':
    app.run(debug=True, host='0.0.0.0', port=5000) 