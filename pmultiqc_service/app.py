#!/usr/bin/env python
"""
PMultiQC Service API - FastAPI Version
A FastAPI-based web service for generating PMultiQC reports from uploaded data files.
"""

import os
import logging
import sys
import uuid
import shutil
import time
import zipfile
import subprocess
import asyncio
import json
import traceback
import threading
import requests
import hashlib
from datetime import datetime
from typing import Dict, Any, List, Optional, Union

from fastapi import FastAPI, File, UploadFile, HTTPException, BackgroundTasks, Form, Request
from fastapi.openapi.utils import get_openapi
from fastapi.responses import JSONResponse, FileResponse, HTMLResponse, RedirectResponse
from fastapi.middleware.cors import CORSMiddleware
from fastapi.templating import Jinja2Templates
from fastapi.staticfiles import StaticFiles
import uvicorn

import sqlite3
from pathlib import Path

# Configuration
UPLOAD_FOLDER = os.environ.get('UPLOAD_FOLDER', '/tmp/pmultiqc_uploads')
OUTPUT_FOLDER = os.environ.get('OUTPUT_FOLDER', '/tmp/pmultiqc_outputs')
HTML_REPORTS_FOLDER = os.environ.get('HTML_REPORTS_FOLDER', '/tmp/pmultiqc_html_reports')
BASE_URL = os.environ.get('BASE_URL', 'http://localhost:5000')

# Ensure directories exist
os.makedirs(UPLOAD_FOLDER, exist_ok=True)
os.makedirs(OUTPUT_FOLDER, exist_ok=True)
os.makedirs(HTML_REPORTS_FOLDER, exist_ok=True)

# Initialize Jinja2 templates
templates = Jinja2Templates(directory="templates")

# Allowed file extensions
ALLOWED_EXTENSIONS = {'zip'}

# Database configuration - use the same directory as uploads for consistency
DATABASE_PATH = os.path.join(UPLOAD_FOLDER, 'jobs.db')

# Configure logging
def setup_logging():
    """Configure logging for Kubernetes environment."""
    log_level = os.environ.get('LOG_LEVEL', 'INFO').upper()
    
    formatter = logging.Formatter(
        '%(asctime)s [%(levelname)s] %(name)s: %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(getattr(logging, log_level))
    console_handler.setFormatter(formatter)
    
    root_logger = logging.getLogger()
    root_logger.setLevel(getattr(logging, log_level))
    
    for handler in root_logger.handlers[:]:
        root_logger.removeHandler(handler)
    
    root_logger.addHandler(console_handler)
    
    logger = logging.getLogger(__name__)
    logger.info(f"Logging configured with level: {log_level}")
    
    return logger

# Initialize logging
logger = setup_logging()

def init_database():
    """Initialize the SQLite database for job storage."""
    try:
        conn = sqlite3.connect(DATABASE_PATH)
        cursor = conn.cursor()
        
        # Create jobs table
        cursor.execute('''
            CREATE TABLE IF NOT EXISTS jobs (
                job_id TEXT PRIMARY KEY,
                status TEXT NOT NULL,
                progress INTEGER DEFAULT 0,
                input_type TEXT,
                accession TEXT,
                started_at TEXT,
                finished_at TEXT,
                error TEXT,
                files_processed INTEGER,
                total_files INTEGER,
                files_downloaded INTEGER,
                processing_stage TEXT,
                console_output TEXT,
                console_errors TEXT,
                html_report_urls TEXT,
                download_url TEXT,
                download_details TEXT,
                created_at TEXT DEFAULT CURRENT_TIMESTAMP,
                updated_at TEXT DEFAULT CURRENT_TIMESTAMP
            )
        ''')
        
        # Create index for faster lookups
        cursor.execute('CREATE INDEX IF NOT EXISTS idx_job_status ON jobs(status)')
        cursor.execute('CREATE INDEX IF NOT EXISTS idx_job_created ON jobs(created_at)')
        
        conn.commit()
        conn.close()
        logger.info(f"Database initialized at {DATABASE_PATH}")
        
    except Exception as e:
        logger.error(f"Error initializing database: {e}")
        raise

def save_job_to_db(job_id: str, job_data: dict):
    """Save or update job data in the database."""
    try:
        conn = sqlite3.connect(DATABASE_PATH)
        cursor = conn.cursor()
        
        # Convert lists to JSON strings for storage
        console_output = json.dumps(job_data.get('console_output', []))
        console_errors = json.dumps(job_data.get('console_errors', []))
        html_report_urls = json.dumps(job_data.get('html_report_urls', {}))
        download_details = json.dumps(job_data.get('download_details', {}))
        
        cursor.execute('''
            INSERT OR REPLACE INTO jobs (
                job_id, status, progress, input_type, accession, started_at, finished_at,
                error, files_processed, total_files, files_downloaded, processing_stage,
                console_output, console_errors, html_report_urls, download_url, download_details, updated_at
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, CURRENT_TIMESTAMP)
        ''', (
            job_id,
            job_data.get('status', 'unknown'),
            job_data.get('progress', 0),
            job_data.get('input_type'),
            job_data.get('accession'),
            job_data.get('started_at'),
            job_data.get('finished_at'),
            job_data.get('error'),
            job_data.get('files_processed'),
            job_data.get('total_files'),
            job_data.get('files_downloaded'),
            job_data.get('processing_stage'),
            console_output,
            console_errors,
            html_report_urls,
            job_data.get('download_url'),
            download_details
        ))
        
        conn.commit()
        conn.close()
        logger.info(f"Job {job_id} saved to database")
        
    except Exception as e:
        logger.error(f"Error saving job {job_id} to database: {e}")

def get_job_from_db(job_id: str) -> Optional[dict]:
    """Retrieve job data from the database."""
    try:
        conn = sqlite3.connect(DATABASE_PATH)
        cursor = conn.cursor()
        
        cursor.execute('SELECT * FROM jobs WHERE job_id = ?', (job_id,))
        row = cursor.fetchone()
        conn.close()
        
        if row:
            # Convert JSON strings back to lists/dicts
            job_data = {
                'job_id': row[0],
                'status': row[1],
                'progress': row[2],
                'input_type': row[3],
                'accession': row[4],
                'started_at': row[5],
                'finished_at': row[6],
                'error': row[7],
                'files_processed': row[8],
                'total_files': row[9],
                'files_downloaded': row[10],
                'processing_stage': row[11],
                'console_output': json.loads(row[12]) if row[12] else [],
                'console_errors': json.loads(row[13]) if row[13] else [],
                'html_report_urls': json.loads(row[14]) if row[14] else {},
                'download_url': row[15],
                'download_details': json.loads(row[16]) if row[16] else {},
                'created_at': row[17],
                'updated_at': row[18]
            }
            logger.info(f"Retrieved job {job_id} from database")
            return job_data
        else:
            logger.info(f"Job {job_id} not found in database")
            return None
            
    except Exception as e:
        logger.error(f"Error retrieving job {job_id} from database: {e}")
        return None

def update_job_progress(job_id: str, status: str, progress: int = None, **kwargs):
    """Update job progress and status in the database."""
    try:
        conn = sqlite3.connect(DATABASE_PATH)
        cursor = conn.cursor()
        
        # Build update query dynamically
        update_fields = ['status = ?', 'updated_at = CURRENT_TIMESTAMP']
        params = [status]
        
        if progress is not None:
            update_fields.append('progress = ?')
            params.append(progress)
        
        for key, value in kwargs.items():
            if key in ['input_type', 'accession', 'started_at', 'finished_at', 'error', 
                      'files_processed', 'total_files', 'download_url', 'processing_stage', 'files_downloaded']:
                update_fields.append(f'{key} = ?')
                params.append(value)
            elif key in ['console_output', 'console_errors', 'html_report_urls', 'download_details']:
                update_fields.append(f'{key} = ?')
                params.append(json.dumps(value))
        
        params.append(job_id)
        
        query = f'UPDATE jobs SET {", ".join(update_fields)} WHERE job_id = ?'
        cursor.execute(query, params)
        
        conn.commit()
        conn.close()
        logger.info(f"Updated job {job_id} in database: status={status}, progress={progress}")
        
    except Exception as e:
        logger.error(f"Error updating job {job_id} in database: {e}")

def cleanup_old_jobs(days_to_keep: int = 30):
    """Clean up old completed/failed jobs from the database."""
    try:
        conn = sqlite3.connect(DATABASE_PATH)
        cursor = conn.cursor()
        
        # Delete jobs older than specified days
        cursor.execute('''
            DELETE FROM jobs 
            WHERE status IN ('completed', 'failed') 
            AND updated_at < datetime('now', '-{} days')
        '''.format(days_to_keep))
        
        deleted_count = cursor.rowcount
        conn.commit()
        conn.close()
        
        if deleted_count > 0:
            logger.info(f"Cleaned up {deleted_count} old jobs from database")
            
    except Exception as e:
        logger.error(f"Error cleaning up old jobs: {e}")

# Initialize database on startup
init_database()

# Clean up old jobs on startup
cleanup_old_jobs(days_to_keep=30)

# Create FastAPI app
app = FastAPI(
    title="PMultiQC Service API",
    description="A FastAPI-based web service for generating PMultiQC reports from uploaded data files and PRIDE datasets.",
    version="1.0.0",
    contact={
        "name": "PMultiQC Team",
        "url": "https://github.com/bigbio/pmultiqc",
    },
    license_info={
        "name": "MIT",
        "url": "https://opensource.org/licenses/MIT",
    },
    docs_url=None,  # Disable automatic docs to use our custom one
    redoc_url=None,  # Disable automatic redoc to avoid conflicts
    openapi_url="/openapi.json",
    # Disable automatic OpenAPI generation to use our custom one
    openapi_tags=None
)

# Add CORS middleware
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Mount static files
app.mount("/static", StaticFiles(directory="templates"), name="static")

# Configure OpenAPI for subpath deployment
def custom_openapi():
    if app.openapi_schema:
        return app.openapi_schema
    
    openapi_schema = get_openapi(
        title=app.title,
        version=app.version,
        description=app.description,
        routes=app.routes,
    )
    
    # Set the correct server URL using BASE_URL environment variable
    # Remove trailing slash if present to avoid double slashes
    base_url = BASE_URL.rstrip('/')
    openapi_schema["servers"] = [
        {
            "url": base_url,
            "description": "Production server"
        }
    ]
    
    app.openapi_schema = openapi_schema
    return app.openapi_schema

app.openapi = custom_openapi

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

def download_pride_file(file_info: Dict, download_dir: str, job_id: str = None) -> str:
    """
    Download a file from PRIDE with detailed progress tracking.
    
    Args:
        file_info: File information from PRIDE API
        download_dir: Directory to save the file
        job_id: Job ID for progress updates
    
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
        
        # Download file with detailed progress tracking
        response = requests.get(http_url, stream=True, timeout=300)
        response.raise_for_status()
        
        # Get total file size for progress calculation
        total_size = int(response.headers.get('content-length', 0))
        downloaded_size = 0
        start_time = time.time()
        
        with open(file_path, 'wb') as f:
            for chunk in response.iter_content(chunk_size=8192):
                if chunk:
                    f.write(chunk)
                    downloaded_size += len(chunk)
                    
                    # Update progress every 1MB or every 5 seconds
                    if downloaded_size % (1024 * 1024) == 0 or (time.time() - start_time) > 5:
                        if job_id and total_size > 0:
                            # Calculate transfer speed and ETA
                            elapsed_time = time.time() - start_time
                            if elapsed_time > 0:
                                speed = downloaded_size / elapsed_time  # bytes per second
                                remaining_bytes = total_size - downloaded_size
                                eta_seconds = remaining_bytes / speed if speed > 0 else 0
                                
                                # Update job progress with detailed transfer info
                                update_job_progress(
                                    job_id, 
                                    'downloading_files', 
                                    progress=None,  # Keep existing progress
                                    processing_stage=f'Downloading {filename}...',
                                    download_details={
                                        'current_file': filename,
                                        'downloaded_bytes': downloaded_size,
                                        'total_bytes': total_size,
                                        'speed_bytes_per_sec': speed,
                                        'eta_seconds': eta_seconds,
                                        'elapsed_seconds': elapsed_time
                                    }
                                )
                        
                        start_time = time.time()  # Reset timer
        
        file_size = os.path.getsize(file_path)
        logger.info(f"Successfully downloaded {filename} ({file_size} bytes)")
        return file_path
        
    except Exception as e:
        logger.error(f"Error downloading file {file_info.get('fileName')}: {e}")
        raise Exception(f"Failed to download {file_info.get('fileName')}: {e}")

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
        public_dir = os.path.join(HTML_REPORTS_FOLDER, report_hash)
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
            base_url = BASE_URL.rstrip('/')
            html_url = f'{base_url}/reports/{report_hash}/{html_rel_path}'
            html_urls[report_name] = html_url
            logger.info(f"HTML report for '{report_name}' available at: {html_url}")
        
        logger.info(f"Generated HTML URLs: {html_urls}")
        return html_urls
        
    except Exception as e:
        logger.error(f"Error copying HTML reports for PRIDE: {e}")
        logger.error(f"Traceback: {traceback.format_exc()}")
        return None

def process_pride_job_async(job_id: str, accession: str, output_dir: str):
    """Process a PRIDE job asynchronously in a separate thread."""
    try:
        # Initialize job in database
        update_job_progress(job_id, 'fetching_files', 10, accession=accession, 
                          processing_stage='Fetching file list from PRIDE database...')
        
        # Fetch files from PRIDE
        logger.info(f"Starting PRIDE processing for job {job_id}, accession {accession}")
        files = get_pride_files(accession)
        
        update_job_progress(job_id, 'filtering_files', 20, 
                          processing_stage='Filtering search engine output files...')
        
        # Filter search engine output files
        search_files = filter_search_files(files)
        
        if not search_files:
            update_job_progress(job_id, 'failed', error='No search engine output files found in PRIDE submission', 
                              finished_at=datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
            return
        
        update_job_progress(job_id, 'downloading_files', 30, 
                          processing_stage=f'Downloading {len(search_files)} files from PRIDE...',
                          total_files=len(search_files), files_downloaded=0)
        
        # Create download directory
        download_dir = os.path.join(UPLOAD_FOLDER, job_id, 'pride_downloads')
        os.makedirs(download_dir, exist_ok=True)
        
        # Download files
        downloaded_files = []
        total_files = len(search_files)
        
        for i, file_info in enumerate(search_files):
            try:
                logger.info(f"Downloading file {i+1}/{total_files}: {file_info.get('fileName')}")
                file_path = download_pride_file(file_info, download_dir, job_id)
                downloaded_files.append(file_path)
                
                # Update progress for downloads (30-70%)
                progress = 30 + int((i + 1) / total_files * 40)
                update_job_progress(job_id, 'downloading_files', progress, 
                                  files_downloaded=i+1, total_files=total_files,
                                  processing_stage=f'Downloaded {i+1}/{total_files} files...')
                
            except Exception as e:
                logger.warning(f"Failed to download {file_info.get('fileName')}: {e}")
                # Continue with other files
        
        if not downloaded_files:
            update_job_progress(job_id, 'failed', error='Failed to download any files from PRIDE', 
                              finished_at=datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
            return
        
        update_job_progress(job_id, 'extracting_files', 70, 
                          processing_stage='Extracting and analyzing downloaded files...')
        
        # Process each downloaded file separately
        all_results = []
        total_processed = 0
        total_files_to_process = len(downloaded_files)
        
        for i, zip_file in enumerate(downloaded_files):
            if zip_file.lower().endswith('.zip'):
                try:
                    file_name = os.path.splitext(os.path.basename(zip_file))[0]
                    logger.info(f"Processing file {i+1}/{total_files_to_process}: {file_name}")
                    
                    # Update progress for processing (70-90%)
                    progress = 70 + int((i + 1) / total_files_to_process * 20)
                    update_job_progress(job_id, 'processing', progress, 
                                      files_processed=total_processed, total_files=total_files_to_process,
                                      processing_stage=f'Processing {file_name} ({i+1}/{total_files_to_process})...')
                    
                    # Create separate directory for this file
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
                    result = run_pmultiqc_with_progress(file_extract_dir, file_output_dir, input_type, quantms_config, job_id)
                    
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
                    else:
                        logger.warning(f"PMultiQC failed for {file_name}: {result.get('message')}")
                        
                except Exception as e:
                    logger.error(f"Error processing {zip_file}: {e}")
                    continue
        
        if not all_results:
            update_job_progress(job_id, 'failed', error='No files were successfully processed', 
                              finished_at=datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
            return
        
        # Create combined zip report
        combined_zip_path = os.path.join(output_dir, f'pmultiqc_report_{job_id}.zip')
        
        with zipfile.ZipFile(combined_zip_path, 'w') as combined_zip:
            for result in all_results:
                if os.path.exists(result['report_path']):
                    combined_zip.write(result['report_path'], os.path.basename(result['report_path']))
        
        # Save job data to disk for persistence
        job_data = {
            'job_id': job_id,
            'status': 'completed',
            'progress': 100,
            'accession': accession,
            'files_processed': total_processed,
            'total_files': len(downloaded_files),
            'results': all_results,
            'finished_at': datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            'download_url': f'{BASE_URL.rstrip("/")}/download-report/{job_id}'
        }
        
        save_job_data_to_disk(job_id, job_data, output_dir)
        
        # Save to database
        save_job_to_db(job_id, job_data)
        
        logger.info(f"Successfully completed PRIDE job {job_id}: {total_processed}/{len(downloaded_files)} files processed")
        
    except Exception as e:
        logger.error(f"Error in PRIDE async processing for job {job_id}: {e}")
        update_job_progress(job_id, 'failed', error=str(e), 
                          finished_at=datetime.now().strftime("%Y-%m-%d %H:%M:%S"))

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

def detect_input_type(upload_path: str) -> tuple:
    """
    Detect the type of input data based on the contents of the extracted zip file.
    
    Returns:
        tuple: (input_type, quantms_config_path) where input_type is one of 'maxquant', 'diann', 'quantms', 'mzidentml', or 'unknown'
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
        update_job_progress(job_id, 'running_pmultiqc', 0, 
                          console_output=[], console_errors=[], 
                          processing_stage='Starting pmultiqc analysis...')
        
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
                    update_job_progress(job_id, 'running_pmultiqc', 0, 
                                      console_output=output_lines, console_errors=error_lines, 
                                      processing_stage='Running PMultiQC...')
                    logger.info(f"MultiQC: {output_line}")
                
                if error:
                    error_line = error.strip()
                    error_lines.append(error_line)
                    update_job_progress(job_id, 'running_pmultiqc', 0, 
                                      console_output=output_lines, console_errors=error_lines, 
                                      processing_stage='Running PMultiQC...')
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
                        update_job_progress(job_id, 'running_pmultiqc', 0, 
                                          console_output=output_lines, console_errors=error_lines, 
                                          processing_stage='Running PMultiQC...')
                        logger.info(f"MultiQC remaining output: {line.strip()}")
            if remaining_error:
                for line in remaining_error.strip().split('\n'):
                    if line.strip():
                        error_lines.append(line.strip())
                        update_job_progress(job_id, 'running_pmultiqc', 0, 
                                          console_output=output_lines, console_errors=error_lines, 
                                          processing_stage='Running PMultiQC...')
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

def process_job_async(job_id: str, extract_path: str, output_dir: str, input_type: str, quantms_config: str):
    """Process a job asynchronously."""
    try:
        # Update job status to running
        update_job_progress(job_id, 'running_pmultiqc', 0, 
                          console_output=[], console_errors=[], 
                          processing_stage='Starting pmultiqc analysis...')
        
        # Run PMultiQC with real-time progress
        logger.info(f"Starting async processing for job {job_id}")
        result = run_pmultiqc_with_progress(extract_path, output_dir, input_type, quantms_config, job_id)
        
        if not result['success']:
            update_job_progress(job_id, 'failed', error=result['message'], 
                              finished_at=datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
            return
        
        # Update job status to creating report
        update_job_progress(job_id, 'creating_report', 75, 
                          processing_stage='Creating zip report...')
        
        # Create zip report
        zip_report_path = os.path.join(output_dir, f'pmultiqc_report_{job_id}.zip')
        if not create_zip_report(output_dir, zip_report_path):
            update_job_progress(job_id, 'failed', error='Failed to create zip report', 
                              finished_at=datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
            return
        
        # Update job to completed
        final_job_data = {
            'status': 'completed',
            'progress': 100,
            'input_type': input_type,
            'processing_stage': 'Analysis completed successfully!',
            'finished_at': datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            'download_url': f'{BASE_URL.rstrip("/")}/download-report/{job_id}'
        }
        
        # Save final job data to database
        save_job_to_db(job_id, final_job_data)
        
        logger.info(f"Successfully completed async job {job_id}")
        
    except Exception as e:
        logger.error(f"Error in async processing for job {job_id}: {e}")
        update_job_progress(job_id, 'failed', error=str(e), 
                          finished_at=datetime.now().strftime("%Y-%m-%d %H:%M:%S"))

# API Routes
@app.get("/", response_class=HTMLResponse)
async def index(request: Request):
    """Main page of index.html, with upload interface and PRIDE submission page

    Parameters:
        request: Request object
        config: Configuration object

    Returns:
        HTMLResponse: index.html template with upload interface and PRIDE submission page
    """
    return templates.TemplateResponse("index.html", {
        "request": request,
        "config": {"BASE_URL": BASE_URL}
    })

@app.get("/submit")
async def submit_pride(request: Request):
    """Direct submission page for PRIDE datasets via URL parameter.

    Parameters:
        request: Request object
        config: Configuration object

    Returns:
        HTMLResponse: submit.html template with PRIDE submission page
    """
    accession = request.query_params.get('accession', '').strip().upper()
    
    if not accession:
        # No accession provided, show the main page
        return templates.TemplateResponse("index.html", {"request": request, "config": {"BASE_URL": BASE_URL}})
    
    # Validate accession format
    if not accession.startswith('PXD'):
        return templates.TemplateResponse("index.html", {"request": request, "config": {"BASE_URL": BASE_URL}, "error": f'Invalid PRIDE accession format: {accession}. Should start with PXD.'})
    
    # Render the submission page with the accession pre-filled
    return templates.TemplateResponse("submit.html", {"request": request, "accession": accession, "config": {"BASE_URL": BASE_URL}})

@app.get("/results/{job_id}")
async def view_results(request: Request, job_id: str):
    """View results page for a specific job.

    Parameters:
        request: Request object
        job_id: Job ID

    Returns:
        HTMLResponse: results.html template with results page
    """
    try:
        # Validate job_id format
        try:
            uuid.UUID(job_id)
        except ValueError:
            return templates.TemplateResponse("index.html", {"request": request, "config": {"BASE_URL": BASE_URL}, "error": f'Invalid job ID format.'})
        
        logger.info(f"Results page request for {job_id}")
        
        # Try to get job from database first
        job_data = get_job_from_db(job_id)
        
        if job_data:
            logger.info(f"Found job {job_id} in database with status: {job_data.get('status')}")
        else:
            logger.warning(f"Job {job_id} NOT found in database")
            
            # Check if job is completed (for backward compatibility with old jobs)
            output_dir = os.path.join(OUTPUT_FOLDER, job_id)
            zip_report_path = os.path.join(output_dir, f'pmultiqc_report_{job_id}.zip')
            
            if os.path.exists(zip_report_path):
                logger.info(f"Job {job_id} not in database but zip file exists")
                
                # Try to recover job data from output directory
                console_output = []
                console_errors = []
                input_type = None
                
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
                        logger.info(f"Recovered job info for job {job_id}: input_type={input_type}")
                    except Exception as e:
                        logger.error(f"Error reading job info file: {e}")
                
                # Create job_data from recovered information
                job_data = {
                    'job_id': job_id,
                    'status': 'completed',
                    'progress': 100,
                    'input_type': input_type,
                    'accession': None,
                    'started_at': None,
                    'finished_at': None,
                    'error': None,
                    'files_processed': None,
                    'total_files': None,
                    'console_output': console_output,
                    'console_errors': console_errors,
                    'html_report_urls': None,
                    'download_url': f'{BASE_URL.rstrip("/")}/download-report/{job_id}'
                }
                
                # Save recovered job to database for future requests
                save_job_to_db(job_id, job_data)
                
                logger.info("Returning fallback response for completed job " + job_id)
            else:
                return templates.TemplateResponse("index.html", {"request": request, "config": {"BASE_URL": BASE_URL}, "error": f'Job {job_id} not found.'})
        
        # Ensure download_url is present in job_data
        if 'download_url' not in job_data:
            job_data['download_url'] = f'{BASE_URL.rstrip("/")}/download-report/{job_id}'
        
        return templates.TemplateResponse("results.html", {"request": request, "job_id": job_id, "job_data": job_data, "config": {"BASE_URL": BASE_URL}})
        
    except Exception as e:
        logger.error(f"Error in view_results for job {job_id}: {e}")
        return templates.TemplateResponse("index.html", {"request": request, "config": {"BASE_URL": BASE_URL}, "error": f'Error loading job {job_id}: {str(e)}'})

@app.get("/reports/{report_hash}/{filename:path}")
async def serve_html_report(report_hash: str, filename: str):
    """Serve HTML reports for online viewing.

    Parameters:
        report_hash: Hash identifier for the report
        filename: Filename to serve (e.g., multiqc_report.html)
    """
    try:
        start_time = time.time()
        logger.info(f"HTML report request: hash={report_hash}, filename={filename}")
        
        # Validate report_hash format (12 character hex)
        if not report_hash or len(report_hash) != 12 or not all(c in '0123456789abcdef' for c in report_hash):
            logger.error(f"Invalid report hash format: {report_hash}")
            raise HTTPException(status_code=400, detail="Invalid report hash")
        
        # Construct file path
        file_path = os.path.join(HTML_REPORTS_FOLDER, report_hash, filename)
        logger.info(f"Looking for file at: {file_path}")
        
        # List contents of the report directory for debugging
        report_dir = os.path.join(HTML_REPORTS_FOLDER, report_hash)
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
            raise HTTPException(status_code=404, detail="Report directory not found")
        
        # Security check: ensure the file is within the HTML reports directory
        real_path = os.path.realpath(file_path)
        html_reports_real = os.path.realpath(HTML_REPORTS_FOLDER)
        
        if not real_path.startswith(html_reports_real):
            logger.error(f"Security violation: {real_path} not within {html_reports_real}")
            raise HTTPException(status_code=403, detail="Access denied")
        
        if not os.path.exists(file_path):
            logger.error(f"File not found: {file_path}")
            raise HTTPException(status_code=404, detail=f"File not found: {filename}")
        
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
            response = FileResponse(
                path=file_path,
                media_type=mime_type,
                headers={"Cache-Control": "max-age=3600"}  # Cache for 1 hour
            )
        else:
            response = FileResponse(path=file_path, media_type=mime_type)
        
        # Log completion time
        elapsed = time.time() - start_time
        logger.info(f"File served in {elapsed:.2f}s: {filename}")
        
        return response
        
    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Error serving HTML report: {e}")
        raise HTTPException(status_code=500, detail="Internal server error")

@app.get("/health")
async def health_check():
    """Health check endpoint."""
    return {
        "status": "healthy",
        "timestamp": datetime.now().isoformat(),
        "service": "pmultiqc-service",
        "base_url": BASE_URL,
        "endpoints": [
            '/upload-async',
            '/generate-report', 
            '/process-pride',
            '/submit-pride',
            '/job-status/{job_id}',
            '/download-report/{job_id}',
            '/reports/{hash}/{filename}',
            '/docs',
            '/swagger.json'
        ]
    }
    
@app.post("/process-pride")
async def process_pride(request: Request):
    """
    Process PRIDE dataset asynchronously
    
    Submit a PRIDE accession number to download and process the dataset.
    The job will be processed asynchronously and you can check the status
    using the returned job_id.
    """
    data = await request.json()
    accession = data.get('accession', '').strip().upper()
    
    # Validate accession format
    if not accession.startswith('PXD'):
        raise HTTPException(status_code=400, detail="Invalid PRIDE accession format. Should start with PXD")
    
    # Generate unique job ID
    job_id = str(uuid.uuid4())
    
    # Create output directory
    output_dir = os.path.join(OUTPUT_FOLDER, job_id)
    os.makedirs(output_dir, exist_ok=True)
    
    # Initialize job status
    update_job_progress(job_id, 'queued', 0, accession=accession)
    
    logger.info(f"Started PRIDE processing for accession {accession} with job ID {job_id}")
    
    return {
        "job_id": job_id,
        "accession": accession,
        "status": "queued",
        "message": "PRIDE processing started successfully",
        "status_url": f'/job-status/{job_id}',
        "redirect_url": f'/results/{job_id}'
    }

@app.post("/upload-async")
async def upload_async(file: UploadFile = File(..., alias="files")):
    """
    Upload and process files asynchronously
    
    Upload a ZIP file containing data files for PMultiQC analysis.
    The job will be processed asynchronously and you can check the status
    using the returned job_id.
    """
    # Check file extension
    if not file.filename or not file.filename.lower().endswith('.zip'):
        raise HTTPException(status_code=400, detail="Only ZIP files are allowed")
    
    # Generate unique job ID
    job_id = str(uuid.uuid4())
    
    # Create job directories
    upload_dir = os.path.join(UPLOAD_FOLDER, job_id)
    output_dir = os.path.join(OUTPUT_FOLDER, job_id)
    os.makedirs(upload_dir, exist_ok=True)
    os.makedirs(output_dir, exist_ok=True)
    
    # Save uploaded file
    filename = file.filename
    zip_path = os.path.join(upload_dir, filename)
    
    logger.info(f"Starting async upload for job {job_id}: {filename}")
    
    # Initialize job in database
    initial_job_data = {
        'job_id': job_id,
        'status': 'uploading', 
        'progress': 0,
        'started_at': datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    }
    save_job_to_db(job_id, initial_job_data)
    
    # Save the uploaded file
    with open(zip_path, "wb") as buffer:
        shutil.copyfileobj(file.file, buffer)
    
    file_size = os.path.getsize(zip_path)
    logger.info(f"Upload completed for job {job_id}: {file_size} bytes")
    
    # Extract the zip file
    extract_path = os.path.join(upload_dir, 'extracted')
    os.makedirs(extract_path, exist_ok=True)
    
    # Update job status to extracting
    update_job_progress(job_id, 'extracting', 25)
    logger.info(f"Extracting ZIP file for job {job_id}")
    
    with zipfile.ZipFile(zip_path, 'r') as zip_ref:
        zip_ref.extractall(extract_path)
    
    # Detect input type
    input_type, quantms_config = detect_input_type(extract_path)
    logger.info(f"Detected input type: {input_type}")
    
    if input_type == 'unknown':
        update_job_progress(job_id, 'failed', error='Could not detect input type', 
                          finished_at=datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        raise HTTPException(status_code=400, detail="Could not detect input type. Please ensure your ZIP file contains valid data files.")
    
    # Start async processing
    update_job_progress(job_id, 'queued', 50, input_type=input_type)
    
    # Start background task for processing
    import threading
    
    thread = threading.Thread(target=process_job_async, args=(job_id, extract_path, output_dir, input_type, quantms_config))
    thread.daemon = True
    thread.start()
    
    return {
        'job_id': job_id,
        'status': 'queued',
        'message': 'File uploaded successfully. Processing started.',
        'status_url': f'/job-status/{job_id}',
        'file_size': file_size
    }

@app.get("/job-status/{job_id}")
async def job_status_api(job_id: str):
    """
    Get job status and progress
    
    Retrieve the current status and progress of a job by its ID.
    Returns detailed information including console output and error messages.
    """
    try:
        # Validate job_id format
        try:
            uuid.UUID(job_id)
        except ValueError:
            raise HTTPException(status_code=400, detail="Invalid job ID")
        
        logger.info(f"Job status request for {job_id}")
        
        # Try to get job from database first
        job_data = get_job_from_db(job_id)
        
        if job_data:
            logger.info(f"Found job {job_id} in database with status: {job_data.get('status')}")
            
            # Ensure all required fields have default values
            job_data_with_defaults = {
                'job_id': job_data.get('job_id', job_id),
                'status': job_data.get('status', 'unknown'),
                'progress': job_data.get('progress', 0),
                'input_type': job_data.get('input_type'),
                'accession': job_data.get('accession'),
                'started_at': job_data.get('started_at'),
                'finished_at': job_data.get('finished_at'),
                'error': job_data.get('error'),
                'files_processed': job_data.get('files_processed'),
                'total_files': job_data.get('total_files'),
                'files_downloaded': job_data.get('files_downloaded'),
                'processing_stage': job_data.get('processing_stage'),
                'console_output': job_data.get('console_output', []),
                'console_errors': job_data.get('console_errors', []),
                'html_report_urls': job_data.get('html_report_urls'),
                'download_url': job_data.get('download_url'),
                'download_details': job_data.get('download_details', {})
            }
            
            return job_data_with_defaults
        else:
            logger.warning(f"Job {job_id} NOT found in database")
        
        # Check if job is completed (for backward compatibility with old jobs)
        output_dir = os.path.join(OUTPUT_FOLDER, job_id)
        zip_report_path = os.path.join(output_dir, f'pmultiqc_report_{job_id}.zip')
        
        if os.path.exists(zip_report_path):
            logger.info(f"Job {job_id} not in database but zip file exists")
            
            # Try to recover job data from output directory
            console_output = []
            console_errors = []
            input_type = None
            
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
                    logger.info(f"Recovered job info for job {job_id}: input_type={input_type}")
                except Exception as e:
                    logger.error(f"Error reading job info file: {e}")
            
            # Fallback response with recovered data
            base_url = BASE_URL.rstrip('/')
            response_data = {
                'job_id': job_id,
                'status': 'completed',
                'progress': 100,
                'input_type': input_type,
                'accession': None,
                'started_at': None,
                'finished_at': None,
                'error': None,
                'files_processed': None,
                'total_files': None,
                'files_downloaded': None,
                'processing_stage': 'Completed',
                'console_output': console_output,
                'console_errors': console_errors,
                'html_report_urls': None,
                'download_url': f'{base_url}/download-report/{job_id}',
                'download_details': {}
            }
            
            # Save recovered job to database for future requests
            save_job_to_db(job_id, response_data)
            
            logger.info("Returning fallback response for completed job " + job_id)
            return response_data
        else:
            raise HTTPException(status_code=404, detail=f"Job {job_id} not found")
            
    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Error in job_status_api for job {job_id}: {e}")
        raise HTTPException(status_code=500, detail=f"Internal server error: {str(e)}")

@app.post("/generate-report")
async def generate_report(file: UploadFile = File(...)):
    """
    Generate a PMultiQC report from uploaded data.
    
    Expected form data:
    - file: ZIP file containing the data to analyze
    
    Returns:
    - JSON response with job ID and status
    """
    try:
        # Check file extension
        if not file.filename or not file.filename.lower().endswith('.zip'):
            raise HTTPException(status_code=400, detail="Only ZIP files are allowed")
        
        # Generate unique job ID
        job_id = str(uuid.uuid4())
        
        # Create job directories
        upload_dir = os.path.join(UPLOAD_FOLDER, job_id)
        output_dir = os.path.join(OUTPUT_FOLDER, job_id)
        os.makedirs(upload_dir, exist_ok=True)
        os.makedirs(output_dir, exist_ok=True)
        
        # Save uploaded file
        filename = file.filename
        zip_path = os.path.join(upload_dir, filename)
        
        with open(zip_path, "wb") as buffer:
            shutil.copyfileobj(file.file, buffer)
        
        # Extract the zip file
        extract_path = os.path.join(upload_dir, 'extracted')
        os.makedirs(extract_path, exist_ok=True)
        
        with zipfile.ZipFile(zip_path, 'r') as zip_ref:
            zip_ref.extractall(extract_path)
        
        # Detect input type
        input_type, quantms_config = detect_input_type(extract_path)
        logger.info(f"Detected input type: {input_type}")
        
        if input_type == 'unknown':
            raise HTTPException(status_code=400, detail="Could not detect input type. Please ensure your ZIP file contains valid data files.")
        
        # Run PMultiQC
        result = run_pmultiqc_with_progress(extract_path, output_dir, input_type, quantms_config, job_id)
        
        if not result['success']:
            raise HTTPException(status_code=500, detail=f"Failed to generate report: {result['message']}")
        
        # Create zip report
        zip_report_path = os.path.join(output_dir, f'pmultiqc_report_{job_id}.zip')
        if not create_zip_report(output_dir, zip_report_path):
            logger.error(f"Failed to create zip report for job {job_id}")
            raise HTTPException(status_code=500, detail="Failed to create zip report. Please try again.")
        
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
            'finished_at_seconds': int(time.time()),
            'finished_at': datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        }
        
        if html_urls:
            job_update['html_report_urls'] = html_urls

        save_job_to_db(job_id, job_update)

        base_url = BASE_URL.rstrip('/')
        response_data = {
            'job_id': job_id,
            'status': 'completed',
            'input_type': input_type,
            'download_url': f'{base_url}/download-report/{job_id}',
            'message': 'Report generated successfully',
            'finished_at': datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        }
        
        if html_urls:
            response_data['html_report_urls'] = html_urls

        return response_data
        
    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Error in generate_report: {e}")
        raise HTTPException(status_code=500, detail="Internal server error")


@app.post("/submit-pride")
async def api_submit_pride(request: Request):
    """API endpoint for direct PRIDE submission."""
    try:
        logger.info(f"API submit-pride called")
        
        data = await request.json()
        if not data or 'accession' not in data:
            raise HTTPException(status_code=400, detail="Missing accession parameter")
        
        accession = data['accession'].strip().upper()
        
        # Validate accession format
        if not accession.startswith('PXD'):
            raise HTTPException(status_code=400, detail="Invalid PRIDE accession format. Should start with PXD")
        
        # Generate unique job ID
        job_id = str(uuid.uuid4())
        
        # Create output directory
        output_dir = os.path.join(OUTPUT_FOLDER, job_id)
        os.makedirs(output_dir, exist_ok=True)
        
        # Initialize job in database
        initial_job_data = {
            'job_id': job_id,
            'status': 'queued',
            'accession': accession,
            'progress': 0,
            'started_at': datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
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
        save_job_to_db(job_id, initial_job_data)
        
        # Start async processing
        thread = threading.Thread(
            target=process_pride_job_async,
            args=(job_id, accession, output_dir)
        )
        thread.daemon = True
        thread.start()
        
        logger.info(f"Started PRIDE processing for accession {accession} with job ID {job_id}")
        
        return {
            'job_id': job_id,
            'accession': accession,
            'status': 'queued',
            'message': 'PRIDE processing started successfully',
            'status_url': f'/job-status/{job_id}',
            'redirect_url': f'/results/{job_id}'
        }
        
    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Error starting PRIDE processing: {e}")
        raise HTTPException(status_code=500, detail=f"Failed to start PRIDE processing: {str(e)}")

@app.get("/download-report/{job_id}")
async def download_report_api(job_id: str):
    """API version of download-report endpoint."""
    return await download_report(job_id)

async def download_report(job_id: str):
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
            raise HTTPException(status_code=400, detail="Invalid job ID")
        
        output_dir = os.path.join(OUTPUT_FOLDER, job_id)
        zip_report_path = os.path.join(output_dir, f'pmultiqc_report_{job_id}.zip')
        
        if not os.path.exists(zip_report_path):
            raise HTTPException(status_code=404, detail="Report not found")
        
        return FileResponse(
            path=zip_report_path,
            filename=f'pmultiqc_report_{job_id}.zip',
            media_type='application/zip'
        )
        
    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Error in download_report: {e}")
        raise HTTPException(status_code=500, detail="Internal server error")

@app.get("/swagger.json")
async def swagger_json():
    """Return OpenAPI specification as JSON."""
    return JSONResponse(content=app.openapi())

@app.get("/openapi.json")
async def openapi_json():
    """Return OpenAPI specification as JSON."""
    return JSONResponse(content=app.openapi())

# Custom docs endpoint to ensure proper OpenAPI URL resolution
@app.get("/docs", include_in_schema=False)
async def custom_docs():
    """Custom docs endpoint that redirects to Swagger UI with correct configuration."""
    from fastapi.openapi.docs import get_swagger_ui_html
    
    return get_swagger_ui_html(
        openapi_url=f"{BASE_URL.rstrip('/')}/openapi.json",
        title=f"{app.title} - Swagger UI",
        swagger_js_url="https://cdn.jsdelivr.net/npm/swagger-ui-dist@5.9.0/swagger-ui-bundle.js",
        swagger_css_url="https://cdn.jsdelivr.net/npm/swagger-ui-dist@5.9.0/swagger-ui.css",
    )

if __name__ == "__main__":
    uvicorn.run(app, host="0.0.0.0", port=5000) 