"""
pmultiqc Service API - FastAPI Version
A FastAPI-based web service for generating pmultiqc reports from uploaded data files.
"""

import hashlib
import json
import logging
import os
import re
import shutil
import subprocess
import sys
import threading
import time
import traceback
import uuid
import zipfile
from datetime import datetime
from pathlib import Path
from typing import Dict, Any, List, Optional

import redis
import requests
import uvicorn
from fastapi import FastAPI, File, UploadFile, HTTPException, Request
from fastapi.middleware.cors import CORSMiddleware
from fastapi.openapi.docs import get_swagger_ui_html
from fastapi.openapi.utils import get_openapi
from fastapi.responses import JSONResponse, FileResponse, HTMLResponse
from fastapi.staticfiles import StaticFiles
from fastapi.templating import Jinja2Templates
from tuspyserver import create_tus_router

# Configuration
# Use environment variables with fallback to current working directory subdirectories
UPLOAD_FOLDER = os.environ.get("UPLOAD_FOLDER", os.path.join(os.getcwd(), "pmultiqc_uploads"))
OUTPUT_FOLDER = os.environ.get("OUTPUT_FOLDER", os.path.join(os.getcwd(), "pmultiqc_outputs"))
HTML_REPORTS_FOLDER = os.environ.get(
    "HTML_REPORTS_FOLDER", os.path.join(os.getcwd(), "pmultiqc_html_reports")
)
BASE_URL = os.environ.get("BASE_URL", "http://localhost:5000")

# Security: File upload limits
MAX_FILE_SIZE = int(os.environ.get("MAX_FILE_SIZE", str(10 * 1024 * 1024 * 1024)))  # 10GB default
MAX_UPLOAD_FILES = int(os.environ.get("MAX_UPLOAD_FILES", "100"))  # Max files in a zip


# PRIDE button visibility configuration
# Priority: Environment variable > Default (False)
def get_pride_button_visible():
    """Get PRIDE button visibility setting from environment variable."""
    env_val = os.environ.get('PRIDE_BUTTON_VISIBLE')
    if env_val is not None:
        return env_val.lower() in ('true', '1', 'yes', 'on')

    # Default to False - users must explicitly enable PRIDE functionality
    return False


PRIDE_BUTTON_VISIBLE = get_pride_button_visible()

# Ensure directories exist
os.makedirs(UPLOAD_FOLDER, exist_ok=True)
os.makedirs(OUTPUT_FOLDER, exist_ok=True)
os.makedirs(HTML_REPORTS_FOLDER, exist_ok=True)

# Initialize Jinja2 templates
templates = Jinja2Templates(directory=str(Path(__file__).parent / "templates"))

# Allowed file extensions
ALLOWED_EXTENSIONS = {"zip"}

# Redis configuration
REDIS_URL = os.environ.get("REDIS_URL", "redis://localhost:6379")
REDIS_USERNAME = os.environ.get("REDIS_USERNAME", None)
REDIS_PASSWORD = os.environ.get("REDIS_PASSWORD", None)
REDIS_DB = int(os.environ.get("REDIS_DB", "0"))
REDIS_KEY_PREFIX = "pmultiqc:job:"


# Configure logging
def setup_logging():
    """Configure logging for Kubernetes environment."""
    log_level = os.environ.get("LOG_LEVEL", "INFO").upper()

    formatter = logging.Formatter(
        "%(asctime)s [%(levelname)s] %(name)s: %(message)s", datefmt="%Y-%m-%d %H:%M:%S"
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

    # Apply the filter to the uvicorn access logger
    logging.getLogger("uvicorn.access").addFilter(NoHealthAccessLogFilter())

    return logger


class NoHealthAccessLogFilter(logging.Filter):
    def filter(self, record):
        return "/health HTTP/" not in record.getMessage()


# Initialize logging
logger = setup_logging()

# Log PRIDE button configuration
logger.info(f"PRIDE button visibility: {PRIDE_BUTTON_VISIBLE}")

def init_redis():
    """Initialize Redis connection for job storage."""
    try:
        # Create Redis connection with authentication if provided
        if REDIS_PASSWORD:
            if REDIS_USERNAME:
                # Username + password authentication
                redis_client = redis.from_url(
                    REDIS_URL,
                    db=REDIS_DB,
                    decode_responses=True,
                    username=REDIS_USERNAME,
                    password=REDIS_PASSWORD,
                )
                logger.info(
                    f"Redis connected with username/password authentication to {REDIS_URL}"
                )
            else:
                # Password-only authentication
                redis_client = redis.from_url(
                    REDIS_URL, db=REDIS_DB, decode_responses=True, password=REDIS_PASSWORD
                )
                logger.info(f"Redis connected with password authentication to {REDIS_URL}")
        else:
            # No authentication
            redis_client = redis.from_url(REDIS_URL, db=REDIS_DB, decode_responses=True)
            logger.info(f"Redis connected without authentication to {REDIS_URL}")

        # Test connection
        redis_client.ping()
        logger.info(f"Redis connection test successful")

        return redis_client

    except Exception as e:
        logger.error(f"Error connecting to Redis: {e}")
        raise


def get_redis_client():
    """Get Redis client instance."""
    try:
        logger.info(f"Connecting to Redis: {REDIS_URL}")
        logger.info(f"Redis password provided: {REDIS_PASSWORD is not None}")
        logger.info(f"Redis username provided: {REDIS_USERNAME is not None}")
        logger.info(f"Redis database: {REDIS_DB}")
        # Security: Never log actual credentials

        if REDIS_PASSWORD:
            if REDIS_USERNAME:
                # Username + password authentication
                redis_client = redis.from_url(
                    REDIS_URL,
                    db=REDIS_DB,
                    decode_responses=True,
                    username=REDIS_USERNAME,
                    password=REDIS_PASSWORD,
                )
            else:
                # Password-only authentication
                redis_client = redis.from_url(
                    REDIS_URL, db=REDIS_DB, decode_responses=True, password=REDIS_PASSWORD
                )
        else:
            # No authentication
            redis_client = redis.from_url(REDIS_URL, db=REDIS_DB, decode_responses=True)

        # Test connection
        redis_client.ping()
        logger.info("Redis connection successful")
        return redis_client

    except Exception as e:
        logger.error(f"Failed to connect to Redis: {e}")
        logger.error(f"Redis connection traceback: {traceback.format_exc()}")
        return None


def save_job_to_db(job_id: str, job_data: dict):
    """Save or update job data in Redis."""
    try:
        redis_client = get_redis_client()

        # Add timestamp if not present
        if "updated_at" not in job_data:
            job_data["updated_at"] = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

        # Convert job data to JSON string
        job_data_json = json.dumps(job_data)

        # Store in Redis with key prefix
        key = f"{REDIS_KEY_PREFIX}{job_id}"
        redis_client.set(key, job_data_json)

        # Set expiration (24 hours)
        redis_client.expire(key, 24 * 3600)

        logger.info(f"Job {job_id} saved to Redis")

    except Exception as e:
        logger.error(f"Error saving job {job_id} to Redis: {e}")


def get_job_from_db(job_id: str) -> Optional[dict]:
    """Retrieve job data from Redis."""
    try:
        redis_client = get_redis_client()

        # Get job data from Redis
        key = f"{REDIS_KEY_PREFIX}{job_id}"
        job_data_json = redis_client.get(key)

        if job_data_json:
            # Parse JSON data
            job_data = json.loads(job_data_json)
            logger.info(f"Retrieved job {job_id} from Redis")
            return job_data
        else:
            logger.info(f"Job {job_id} not found in Redis")
            return None

    except Exception as e:
        logger.error(f"Error retrieving job {job_id} from Redis: {e}")
        return None


def update_job_progress(job_id: str, status: str, progress: int = None, **kwargs):
    """Update job progress and status in Redis."""
    try:
        redis_client = get_redis_client()

        if redis_client is None:
            logger.error(f"Cannot update job {job_id} progress: Redis client is None")
            return

        # Get existing job data
        key = f"{REDIS_KEY_PREFIX}{job_id}"
        existing_data_json = redis_client.get(key)

        if existing_data_json:
            job_data = json.loads(existing_data_json)
        else:
            job_data = {}

        # Update fields
        job_data["status"] = status
        job_data["updated_at"] = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

        if progress is not None:
            job_data["progress"] = progress

        # Update additional fields
        for key_name, value in kwargs.items():
            job_data[key_name] = value

        # Save updated data back to Redis
        job_data_json = json.dumps(job_data)
        redis_client.set(key, job_data_json)

        # Set expiration (24 hours)
        redis_client.expire(key, 24 * 3600)

        logger.info(f"Updated job {job_id} in Redis: status={status}, progress={progress}")
        logger.info(f"Job data keys: {list(job_data.keys())}")
        if "html_report_urls" in job_data:
            logger.info(f"HTML report URLs in job data: {job_data['html_report_urls']}")

    except Exception as e:
        logger.error(f"Error updating job {job_id} in Redis: {e}")


def cleanup_old_jobs(days_to_keep: int = 30):
    """Clean up old completed/failed jobs from Redis."""
    try:
        redis_client = get_redis_client()

        # Get all job keys
        pattern = f"{REDIS_KEY_PREFIX}*"
        job_keys = redis_client.keys(pattern)

        deleted_count = 0
        cutoff_time = datetime.now().timestamp() - (days_to_keep * 24 * 3600)

        for key in job_keys:
            try:
                job_data_json = redis_client.get(key)
                if job_data_json:
                    job_data = json.loads(job_data_json)

                    # Check if job is old and completed/failed
                    if (
                        job_data.get("status") in ["completed", "failed"]
                        and "updated_at" in job_data
                    ):
                        try:
                            job_time = datetime.strptime(
                                job_data["updated_at"], "%Y-%m-%d %H:%M:%S"
                            ).timestamp()
                            if job_time < cutoff_time:
                                redis_client.delete(key)
                                deleted_count += 1
                        except ValueError:
                            # If timestamp parsing fails, skip
                            continue

            except Exception as e:
                logger.warning(f"Error processing job key {key}: {e}")
                continue

        if deleted_count > 0:
            logger.info(f"Cleaned up {deleted_count} old jobs from Redis")

    except Exception as e:
        logger.error(f"Error cleaning up old jobs: {e}")


# Initialize Redis on startup
init_redis()

# Clean up old jobs on startup
cleanup_old_jobs(days_to_keep=30)

# Create FastAPI app
app = FastAPI(
    title="pmultiqc Service API",
    description="A FastAPI-based web service for generating pmultiqc reports from uploaded data files and PRIDE datasets.",
    version="1.0.0",
    contact={
        "name": "pmultiqc Team",
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
    openapi_tags=None,
)

# Add CORS middleware
# Security: Configure CORS based on environment
# Get allowed origins from environment variable or use secure defaults
ALLOWED_ORIGINS = os.environ.get("ALLOWED_ORIGINS", "").split(",")
if not ALLOWED_ORIGINS or ALLOWED_ORIGINS == [""]:
    # Default to localhost for development
    ALLOWED_ORIGINS = ["http://localhost", "http://localhost:5000", "http://127.0.0.1", "http://127.0.0.1:5000"]
    logger.warning("Using default CORS origins for development. Set ALLOWED_ORIGINS environment variable for production.")

app.add_middleware(
    CORSMiddleware,
    allow_origins=ALLOWED_ORIGINS,
    allow_credentials=True,
    allow_methods=["GET", "POST", "PUT", "PATCH", "DELETE", "OPTIONS", "HEAD"],  # TUS needs PATCH and HEAD
    allow_headers=["*"],  # TUS needs custom headers (Upload-*, Tus-*)
    expose_headers=["*"],  # TUS needs to expose custom headers
)


# TUS Upload Protocol Implementation (using tuspyserver)
# Callback function for when TUS upload completes
def handle_upload_complete(file_path: str, metadata: dict):
    """
    Called by tuspyserver when upload completes.
    
    Args:
        file_path: Path to the uploaded file
        metadata: Upload metadata (filename, filetype, etc.)
    """
    try:
        # Extract metadata
        filename = metadata.get("filename", "upload.zip")
        filetype = metadata.get("filetype", "application/zip")
        
        logger.info(f"TUS upload complete: file={file_path}, filename={filename}")
        
        # Validate filename
        if not filename.lower().endswith(".zip"):
            logger.error(f"Invalid file type: {filename}")
            # Clean up uploaded file
            if os.path.exists(file_path):
                os.remove(file_path)
            raise ValueError(f"Only ZIP files are allowed. Received: {filename}")
        
        # Get file size
        file_size = os.path.getsize(file_path)
        
        # Generate job ID
        job_id = str(uuid.uuid4())
        
        # Create job directories
        job_upload_dir = os.path.join(UPLOAD_FOLDER, job_id)
        output_dir = os.path.join(OUTPUT_FOLDER, job_id)
        os.makedirs(job_upload_dir, exist_ok=True)
        os.makedirs(output_dir, exist_ok=True)
        
        # Move uploaded file to job directory
        final_zip_path = os.path.join(job_upload_dir, filename)
        shutil.move(file_path, final_zip_path)
        
        logger.info(f"Moved upload to {final_zip_path}, job_id={job_id}")
        
        # Initialize job in database
        initial_job_data = {
            "job_id": job_id,
            "status": "extracting",
            "progress": 25,
            "started_at": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "filename": filename,
            "file_size": file_size,
        }
        save_job_to_db(job_id, initial_job_data)
        
        # Extract and process
        extract_path = os.path.join(job_upload_dir, "extracted")
        os.makedirs(extract_path, exist_ok=True)
        
        validate_and_extract_zip(final_zip_path, extract_path, file_size)
        
        # Detect input type
        input_type, quantms_config = detect_input_type(extract_path)
        logger.info(f"Detected input type: {input_type}")
        
        if input_type == "unknown":
            update_job_progress(
                job_id,
                "failed",
                error="Could not detect input type",
                finished_at=datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            )
        else:
            # Start async processing
            update_job_progress(job_id, "queued", 50, input_type=input_type)
            thread = threading.Thread(
                target=process_job_async,
                args=(job_id, extract_path, output_dir, input_type, quantms_config),
            )
            thread.daemon = True
            thread.start()
        
        logger.info(f"Job {job_id} processing started")
        
    except Exception as e:
        logger.error(f"Error handling upload completion: {e}")
        logger.error(traceback.format_exc())
        raise


# Middleware to fix Location headers for TUS when behind ingress
@app.middleware("http")
async def fix_tus_location_header(request: Request, call_next):
    response = await call_next(request)
    
    # Fix Location header for TUS responses
    if "Location" in response.headers and "/files/" in response.headers["Location"]:
        location = response.headers["Location"]
        # If it's a relative path, prepend BASE_URL
        if location.startswith("/files/"):
            # Use BASE_URL from config
            base = BASE_URL.rstrip("/")
            response.headers["Location"] = f"{base}{location}"
            logger.debug(f"Rewrote Location header: {location} -> {response.headers['Location']}")
    
    return response


# Mount TUS upload router
# This handles resumable file uploads via TUS protocol
# Note: In production, ingress rewrites /pride/services/pmultiqc/files/* to /files/*
tus_router = create_tus_router(
    prefix="/files",  # Endpoint: /files
    files_dir=UPLOAD_FOLDER,  # Where to store uploads
    max_size=MAX_FILE_SIZE,  # 10GB default
    on_upload_complete=handle_upload_complete,  # Callback function
    days_to_keep=7,  # Auto-cleanup after 7 days
)
app.include_router(tus_router)
logger.info(f"TUS upload router mounted at /files with max size {MAX_FILE_SIZE / (1024**3):.1f} GB")

# Mount static files
app.mount("/static", StaticFiles(directory=str(Path(__file__).parent / "templates")), name="static")


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
    base_url = BASE_URL.rstrip("/")
    openapi_schema["servers"] = [{"url": base_url, "description": "Production server"}]

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
        all_files = []
        page = 0
        page_size = 100

        logger.info(f"Fetching files from PRIDE for accession: {accession}")

        while True:
            url = f"https://www.ebi.ac.uk/pride/ws/archive/v3/projects/{accession}/files"
            params = {"page": page, "pageSize": page_size}

            response = requests.get(url, params=params, timeout=30)
            response.raise_for_status()

            files = response.json()

            # If we get an empty array, we've reached the end
            if not files:
                logger.info(f"No more files found on page {page}")
                break

            all_files.extend(files)
            logger.info(f"Retrieved {len(files)} files from page {page}")

            # If we got fewer files than page_size, we've reached the end
            if len(files) < page_size:
                break

            page += 1

        logger.info(f"Retrieved {len(all_files)} total files from PRIDE for {accession}")
        return all_files

    except requests.RequestException as e:
        logger.error(f"Error fetching PRIDE files for {accession}: {e}")
        raise Exception(f"Failed to fetch files from PRIDE: {e}")


def filter_search_files(files: List[Dict]) -> tuple[List[Dict], bool]:
    """
    Filter files to include search engine output files and peak list files for COMPLETE submissions.

    Args:
        files: List of file dictionaries from PRIDE API

    Returns:
        tuple: (filtered_files, is_complete) where filtered_files is the list of files to download
               and is_complete indicates if this is a COMPLETE submission (mzIdentML + peak lists)
    """
    search_files = []
    peak_list_files = []
    is_complete = False

    # Log all file categories for debugging
    file_categories = set()
    for file_info in files:
        category = file_info.get("fileCategory", {}).get("value", "UNKNOWN")
        file_categories.add(category)
    logger.info(f"Found file categories in PRIDE submission: {file_categories}")

    # First pass: identify search engine output files and check for mzIdentML
    has_mzidentml = False
    for file_info in files:
        original_filename = file_info.get("fileName", "")
        filename_lower = original_filename.lower()
        file_category = file_info.get("fileCategory", {}).get("value", "")

        # Check if it's a search engine output file
        if (
            file_category in ["SEARCH", "RESULT"]
            or "report" in filename_lower
            or "evidence" in filename_lower
            or "peptides" in filename_lower
            or "proteingroups" in filename_lower
            or "msms" in filename_lower
            or filename_lower.endswith(".mzid")
            or filename_lower.endswith(".mzid.gz")
            or filename_lower.endswith(".mzid.zip")
        ):

            search_files.append(file_info)
            logger.info(
                f"Found search engine output file: {original_filename} (category: {file_category})"
            )

            # Check if this is an mzIdentML file (including compressed)
            if (
                filename_lower.endswith(".mzid")
                or filename_lower.endswith(".mzid.gz")
                or filename_lower.endswith(".mzid.zip")
            ):
                has_mzidentml = True
                logger.info(f"Found mzIdentML file: {original_filename}")

    # If we have mzIdentML files, this might be a COMPLETE submission
    if has_mzidentml:
        logger.info("Detected mzIdentML files - checking for associated peak list files")

        # Second pass: look for peak list files (mgf, mzml) that might be associated
        for file_info in files:
            original_filename = file_info.get("fileName", "")
            filename_lower = original_filename.lower()

            # Check for peak list files (case-insensitive, including compressed)
            if (
                filename_lower.endswith(".mgf")
                or filename_lower.endswith(".mgf.gz")
                or filename_lower.endswith(".mgf.zip")
                or filename_lower.endswith(".mzml")
                or filename_lower.endswith(".mzml.gz")
                or filename_lower.endswith(".mzml.zip")
                or "peak" in filename_lower
                or "spectrum" in filename_lower
            ):

                # Only include if it's not already in search_files
                if file_info not in search_files:
                    peak_list_files.append(file_info)
                    logger.info(f"Found peak list file: {original_filename}")

        # If we found both mzIdentML and peak list files, mark as COMPLETE
        if peak_list_files:
            is_complete = True
            logger.info(
                f"Detected COMPLETE submission with {len(search_files)} search files and {len(peak_list_files)} peak list files"
            )

    # Combine all files to download
    all_files = search_files + peak_list_files

    logger.info(
        f"Found {len(search_files)} search engine output files and {len(peak_list_files)} peak list files"
    )
    logger.info(f"Total files to download: {len(all_files)}")

    # Log some example filenames for debugging
    if all_files:
        example_files = [f.get("fileName", "") for f in all_files[:5]]
        logger.info(f"Example files to download: {example_files}")
    else:
        logger.warning("No files found to download - this might not be a proteomics dataset")
        # Log some example filenames from the original list for debugging
        example_all_files = [f.get("fileName", "") for f in files[:10]]
        logger.info(f"Example files in submission: {example_all_files}")

    return all_files, is_complete


def download_pride_file(file_info: Dict, download_dir: str, job_id: str = None) -> str:
    """
    Download a file from PRIDE with detailed progress tracking and handle compression.

    Args:
        file_info: File information from PRIDE API
        download_dir: Directory to save the file
        job_id: Job ID for progress updates

    Returns:
        Path to the downloaded file (uncompressed if it was compressed)
    """
    try:
        # Get FTP URL (preferred over Aspera for reliability)
        ftp_url = None
        for location in file_info.get("publicFileLocations", []):
            if location.get("accession") == "PRIDE:0000469":  # FTP Protocol
                ftp_url = location.get("value")
                break

        if not ftp_url:
            raise Exception(f"No FTP URL found for file {file_info.get('fileName')}")

        # Convert FTP URL to HTTP URL for easier download
        http_url = ftp_url.replace("ftp://", "https://")

        filename = file_info.get("fileName")
        file_path = os.path.join(download_dir, filename)

        logger.info(f"Downloading {filename} from {http_url}")

        # Download file with detailed progress tracking
        response = requests.get(http_url, stream=True, timeout=300)
        response.raise_for_status()

        # Get total file size for progress calculation
        total_size = int(response.headers.get("content-length", 0))
        downloaded_size = 0
        start_time = time.time()

        with open(file_path, "wb") as f:
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
                                    "downloading_files",
                                    progress=None,  # Keep existing progress
                                    processing_stage=f"Downloading {filename}...",
                                    download_details={
                                        "current_file": filename,
                                        "downloaded_bytes": downloaded_size,
                                        "total_bytes": total_size,
                                        "speed_bytes_per_sec": speed,
                                        "eta_seconds": eta_seconds,
                                        "elapsed_seconds": elapsed_time,
                                    },
                                )

                        start_time = time.time()  # Reset timer

        file_size = os.path.getsize(file_path)
        logger.info(f"Successfully downloaded {filename} ({file_size} bytes)")

        # Handle compression if the file is compressed
        final_file_path = file_path
        if filename.lower().endswith(".gz"):
            logger.info(f"Decompressing gzipped file: {filename}")
            import gzip

            decompressed_path = file_path[:-3]  # Remove .gz extension
            with gzip.open(file_path, "rb") as f_in:
                with open(decompressed_path, "wb") as f_out:
                    shutil.copyfileobj(f_in, f_out)
            # Remove the compressed file
            os.remove(file_path)
            final_file_path = decompressed_path
            logger.info(f"Decompressed {filename} to {os.path.basename(decompressed_path)}")

        elif filename.lower().endswith(".zip"):
            logger.info(f"Extracting zipped file: {filename}")
            with zipfile.ZipFile(file_path, "r") as zip_ref:
                # Extract to the same directory
                zip_ref.extractall(download_dir)
            # Remove the zip file
            os.remove(file_path)
            # Return the directory path since zip files contain multiple files
            final_file_path = download_dir
            logger.info(f"Extracted {filename} to directory")

        return final_file_path

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
        console_output = job_data.get("console_output", [])
        if console_output:
            console_output_file = os.path.join(output_dir, "console_output.txt")
            with open(console_output_file, "w", encoding="utf-8") as f:
                for line in console_output:
                    f.write(line + "\n")
            logger.info(f"Saved {len(console_output)} console output lines for job {job_id}")

        # Save console errors
        console_errors = job_data.get("console_errors", [])
        if console_errors:
            console_errors_file = os.path.join(output_dir, "console_errors.txt")
            with open(console_errors_file, "w", encoding="utf-8") as f:
                for line in console_errors:
                    f.write(line + "\n")
            logger.info(f"Saved {len(console_errors)} console error lines for job {job_id}")

        # Save job info
        job_info = {
            "job_id": job_id,
            "status": job_data.get("status"),
            "input_type": job_data.get("input_type"),
            "accession": job_data.get("accession"),
            "results": job_data.get("results", []),
            "finished_at": job_data.get("finished_at"),
            "files_processed": job_data.get("files_processed"),
            "total_files": job_data.get("total_files"),
            "files_downloaded": job_data.get("files_downloaded"),
            "html_report_urls": job_data.get("html_report_urls"),
        }

        job_info_file = os.path.join(output_dir, "job_info.json")
        with open(job_info_file, "w", encoding="utf-8") as f:
            json.dump(job_info, f, indent=2)
        logger.info(f"Saved job info for job {job_id} to {job_info_file}")

    except Exception as e:
        logger.error(f"Error saving job data to disk for job {job_id}: {e}")


def copy_html_report_for_online_viewing(
    output_path: str, job_id: str, accession: str = None
) -> Optional[Dict[str, str]]:
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
        logger.info(
            f"Starting copy_html_report_for_online_viewing: output_path={output_path}, job_id={job_id}, accession={accession}"
        )

        # Log the contents of the output directory first
        logger.info(f"Contents of output directory {output_path}:")
        logger.info(f"Output path exists: {os.path.exists(output_path)}")
        logger.info(
            f"Output path is directory: {os.path.isdir(output_path) if os.path.exists(output_path) else 'N/A'}"
        )
        if os.path.exists(output_path):
            for root, _, files in os.walk(output_path):
                level = root.replace(output_path, "").count(os.sep)
                indent = " " * 2 * level
                logger.info(f"{indent}{os.path.basename(root)}/")
                subindent = " " * 2 * (level + 1)
                for file in files:
                    logger.info(f"{subindent}{file}")
                    if file == "multiqc_report.html":
                        logger.info(f"{subindent}*** FOUND HTML REPORT: {file} ***")
        else:
            logger.error(f"Output directory does not exist: {output_path}")
            return None

        # Look for all HTML report files
        html_reports = []
        logger.info(f"Searching for HTML reports in: {output_path}")
        for root, _, files in os.walk(output_path):
            logger.info(f"Checking directory: {root}")
            logger.info(f"Files in {root}: {files}")
            for file in files:
                if file == "multiqc_report.html":
                    html_report_path = os.path.join(root, file)
                    # Get the relative path from output_path to determine report name
                    rel_path = os.path.relpath(root, output_path)
                    # Use "." as root if the report is directly in output_path
                    if rel_path == ".":
                        rel_path = "root"
                    html_reports.append((rel_path, html_report_path))
                    logger.info(
                        f"Found HTML report: {html_report_path} (relative path: {rel_path})"
                    )

        if not html_reports:
            logger.warning(f"No HTML reports found in {output_path}")
            logger.info("Searching for any .html files in the output directory...")
            for root, _, files in os.walk(output_path):
                for file in files:
                    if file.endswith(".html"):
                        logger.info(f"Found HTML file: {os.path.join(root, file)}")
            return None

        # Create a unique hash for the report URL
        hash_input = f"{job_id}_{accession}" if accession else job_id
        report_hash = hashlib.sha256(hash_input.encode()).hexdigest()[:12]
        logger.info(f"Generated report hash: {report_hash} (from input: {hash_input})")

        # Create public directory for this report
        public_dir = os.path.join(HTML_REPORTS_FOLDER, report_hash)
        try:
            os.makedirs(public_dir, exist_ok=True)
            logger.info(f"Created public directory: {public_dir}")
            logger.info(f"HTML_REPORTS_FOLDER: {HTML_REPORTS_FOLDER}")
            logger.info(f"Report hash: {report_hash}")

            # Test write permissions
            test_file = os.path.join(public_dir, "test_write.tmp")
            with open(test_file, "w") as f:
                f.write("test")
            os.remove(test_file)
            logger.info(f"Write permissions OK for {public_dir}")

        except Exception as e:
            logger.error(f"Failed to create or write to public directory {public_dir}: {e}")
            logger.error(f"Traceback: {traceback.format_exc()}")
            return None

        # Create a job file for cleanup tracking
        job_file_path = os.path.join(public_dir, ".job_id")
        with open(job_file_path, "w") as f:
            f.write(job_id)

        # Copy all files from output directory to public directory
        files_copied = 0
        logger.info(f"Starting to copy files from {output_path} to {public_dir}")
        for root, _, files in os.walk(output_path):
            for file in files:
                src_path = os.path.join(root, file)
                rel_path = os.path.relpath(src_path, output_path)
                dst_path = os.path.join(public_dir, rel_path)

                # Ensure destination directory exists
                os.makedirs(os.path.dirname(dst_path), exist_ok=True)

                # Copy file
                try:
                    shutil.copy2(src_path, dst_path)
                    files_copied += 1
                    if file == "multiqc_report.html":
                        logger.info(f"Copied HTML report: {src_path} to {dst_path}")
                    else:
                        logger.info(f"Copied {src_path} to {dst_path}")
                except Exception as e:
                    logger.error(f"Failed to copy {src_path} to {dst_path}: {e}")
                    if file == "multiqc_report.html":
                        logger.error(f"CRITICAL: Failed to copy HTML report file!")
                        return None

        logger.info(f"Total files copied: {files_copied}")

        # Create URLs for each HTML report
        html_urls = {}
        for report_name, html_report_path in html_reports:
            # Find the relative path to the HTML report for the URL
            html_rel_path = os.path.relpath(html_report_path, output_path)
            # Use base URL to generate absolute URLs
            base_url = BASE_URL.rstrip("/")
            html_url = f"{base_url}/reports/{report_hash}/{html_rel_path}"
            html_urls[report_name] = html_url
            logger.info(f"HTML report for '{report_name}' available at: {html_url}")
            logger.info(f"  - report_name: {report_name}")
            logger.info(f"  - html_report_path: {html_report_path}")
            logger.info(f"  - html_rel_path: {html_rel_path}")
            logger.info(f"  - base_url: {base_url}")

        logger.info(f"Generated HTML URLs: {html_urls}")
        return html_urls

    except Exception as e:
        logger.error(f"Error copying HTML reports for PRIDE: {e}")
        logger.error(f"Traceback: {traceback.format_exc()}")
        return None


def process_pride_job_async(job_id: str, accession: str, output_dir: str):
    """Process a PRIDE job asynchronously in a separate thread."""
    logger.info(f"Starting process_pride_job_async for job {job_id}")
    try:
        # Initialize job in database
        update_job_progress(
            job_id,
            "fetching_files",
            10,
            accession=accession,
            processing_stage="Fetching file list from PRIDE database...",
        )

        # Fetch files from PRIDE
        logger.info(f"Starting PRIDE processing for job {job_id}, accession {accession}")
        files = get_pride_files(accession)

        update_job_progress(
            job_id,
            "filtering_files",
            20,
            processing_stage="Filtering search engine output files and peak lists...",
        )

        # Filter search engine output files
        search_files, is_complete = filter_search_files(files)

        if not search_files:
            update_job_progress(
                job_id,
                "failed",
                error="No search engine output files found in PRIDE submission",
                finished_at=datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            )
            return

        # Update progress with submission type information
        submission_type_msg = (
            "COMPLETE submission (mzIdentML + peak lists)"
            if is_complete
            else "Standard submission"
        )
        update_job_progress(
            job_id,
            "downloading_files",
            30,
            processing_stage=f"Downloading {len(search_files)} files from PRIDE ({submission_type_msg})...",
            total_files=len(search_files),
            files_downloaded=0,
        )

        # Create download directory
        download_dir = os.path.join(UPLOAD_FOLDER, job_id, "pride_downloads")
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
                update_job_progress(
                    job_id,
                    "downloading_files",
                    progress,
                    files_downloaded=i + 1,
                    total_files=total_files,
                    processing_stage=f"Downloaded {i+1}/{total_files} files...",
                )

            except Exception as e:
                logger.warning(f"Failed to download {file_info.get('fileName')}: {e}")
                # Continue with other files

        if not downloaded_files:
            update_job_progress(
                job_id,
                "failed",
                error="Failed to download any files from PRIDE",
                finished_at=datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            )
            return

        update_job_progress(
            job_id,
            "extracting_files",
            70,
            processing_stage="Extracting and analyzing downloaded files...",
        )

        # Process downloaded files based on submission type
        all_results = []
        total_processed = 0

        if is_complete:
            # For COMPLETE submissions, process all files together
            logger.info(
                "Processing COMPLETE submission - extracting all files to single directory"
            )

            # Create a single extraction directory for all files
            combined_extract_dir = os.path.join(download_dir, "extracted_combined")
            combined_output_dir = os.path.join(output_dir, "report_combined")

            os.makedirs(combined_extract_dir, exist_ok=True)
            os.makedirs(combined_output_dir, exist_ok=True)

            # Copy all downloaded files to the combined directory
            for downloaded_file in downloaded_files:
                if os.path.isdir(downloaded_file):
                    # If it's a directory (from zip extraction), copy all contents
                    logger.info(f"Copying directory contents from {downloaded_file}")
                    for root, _, files in os.walk(downloaded_file):
                        for file in files:
                            src_path = os.path.join(root, file)
                            rel_path = os.path.relpath(src_path, downloaded_file)
                            dst_path = os.path.join(combined_extract_dir, rel_path)
                            os.makedirs(os.path.dirname(dst_path), exist_ok=True)
                            shutil.copy2(src_path, dst_path)
                else:
                    # If it's a file, copy it directly
                    filename = os.path.basename(downloaded_file)
                    dst_path = os.path.join(combined_extract_dir, filename)
                    shutil.copy2(downloaded_file, dst_path)
                    logger.info(f"Copied file {filename} to combined directory")

            # Detect input type for the combined files
            input_type, quantms_config = detect_input_type(combined_extract_dir)
            logger.info(f"Detected input type for COMPLETE submission: {input_type}")

            if input_type == "unknown":
                update_job_progress(
                    job_id,
                    "failed",
                    error="Could not detect input type for COMPLETE submission",
                    finished_at=datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                )
                return

            # Run pmultiqc on the combined files
            update_job_progress(
                job_id,
                "processing",
                80,
                processing_stage="Processing COMPLETE submission with mzIdentML and peak lists...",
            )

            logger.info(f"Starting run_pmultiqc_with_progress for job {job_id}")
            result = run_pmultiqc_with_progress(
                combined_extract_dir, combined_output_dir, input_type, quantms_config, job_id
            )
            logger.info(
                f"run_pmultiqc_with_progress completed for job {job_id}: success={result.get('success')}"
            )

            if result["success"]:
                # Create zip report for the combined analysis
                logger.info(f"Creating zip report for COMPLETE submission for job {job_id}")
                zip_report_path = os.path.join(
                    combined_output_dir, f"pmultiqc_report_cpmplete_{job_id}.zip"
                )
                if create_zip_report(combined_output_dir, zip_report_path):
                    all_results.append(
                        {
                            "file_name": "cpmplete_combined",
                            "input_type": input_type,
                            "report_path": zip_report_path,
                            "output": result.get("output", []),
                            "errors": result.get("errors", []),
                        }
                    )
                    total_processed = 1
                    logger.info(
                        f"Successfully created zip report for COMPLETE submission for job {job_id}"
                    )
                else:
                    logger.warning("Failed to create zip report for COMPLETE submission")
            else:
                logger.warning(f"pmultiqc failed for COMPLETE submission: {result.get('message')}")
        else:
            # For regular submissions, process each file separately
            total_files_to_process = len(downloaded_files)

            for i, downloaded_file in enumerate(downloaded_files):
                try:
                    logger.info(
                        f"Processing downloaded file {i+1}/{len(downloaded_files)}: {downloaded_file}"
                    )
                    logger.info(f"File exists: {os.path.exists(downloaded_file)}")
                    logger.info(
                        f"Is directory: {os.path.isdir(downloaded_file) if os.path.exists(downloaded_file) else 'N/A'}"
                    )

                    # Determine file name and type
                    if os.path.isdir(downloaded_file):
                        # If it's a directory (from zip extraction), use directory name
                        file_name = os.path.basename(downloaded_file)
                        file_extract_dir = downloaded_file
                        logger.info(f"Using directory as-is: {file_name} -> {file_extract_dir}")
                    else:
                        # If it's a file, extract the name without extension
                        file_name = os.path.splitext(os.path.basename(downloaded_file))[0]
                        # Create extraction directory and copy file
                        file_extract_dir = os.path.join(download_dir, f"extracted_{file_name}")
                        os.makedirs(file_extract_dir, exist_ok=True)
                        shutil.copy2(
                            downloaded_file,
                            os.path.join(file_extract_dir, os.path.basename(downloaded_file)),
                        )
                        logger.info(f"Extracted file: {file_name} -> {file_extract_dir}")

                    logger.info(f"Processing file {i+1}/{total_files_to_process}: {file_name}")

                    # Update progress for processing (70-90%)
                    progress = 70 + int((i + 1) / total_files_to_process * 20)
                    update_job_progress(
                        job_id,
                        "processing",
                        progress,
                        files_processed=total_processed,
                        total_files=total_files_to_process,
                        processing_stage=f"Processing {file_name} ({i+1}/{total_files_to_process})...",
                    )

                    # Create output directory for this file
                    file_output_dir = os.path.join(output_dir, f"report_{file_name}")
                    os.makedirs(file_output_dir, exist_ok=True)
                    logger.info(f"Created output directory: {file_output_dir}")

                    # Detect input type for this file
                    input_type, quantms_config = detect_input_type(file_extract_dir)
                    logger.info(f"Detected input type for {file_name}: {input_type}")

                    # Log files found in the directory for debugging
                    try:
                        files_in_dir = os.listdir(file_extract_dir)
                        logger.info(f"Files in {file_name} directory: {files_in_dir}")
                    except Exception as e:
                        logger.warning(f"Could not list files in {file_extract_dir}: {e}")

                    if input_type == "unknown":
                        logger.warning(f"Could not detect input type for {file_name}")
                        continue

                    # Run pmultiqc on this file
                    logger.info(
                        f"Starting run_pmultiqc_with_progress for job {job_id}, file {file_name}"
                    )
                    result = run_pmultiqc_with_progress(
                        file_extract_dir, file_output_dir, input_type, quantms_config, job_id
                    )
                    logger.info(
                        f"run_pmultiqc_with_progress completed for job {job_id}, file {file_name}: success={result.get('success')}"
                    )

                    if result["success"]:
                        # Create zip report for this file
                        zip_report_path = os.path.join(
                            file_output_dir, f"pmultiqc_report_{job_id}.zip"
                        )
                        if create_zip_report(file_output_dir, zip_report_path):
                            all_results.append(
                                {
                                    "file_name": file_name,
                                    "input_type": input_type,
                                    "report_path": zip_report_path,
                                    "output": result.get("output", []),
                                    "errors": result.get("errors", []),
                                }
                            )
                            total_processed += 1
                            logger.info(
                                f"Successfully processed file {file_name}, total_processed now: {total_processed}"
                            )
                    else:
                        logger.warning(f"pmultiqc failed for {file_name}: {result.get('message')}")

                except Exception as e:
                    logger.error(f"Error processing {downloaded_file}: {e}")
                    logger.error(f"Traceback: {traceback.format_exc()}")
                    continue

        logger.info(
            f"File processing loop completed. total_processed: {total_processed}, all_results count: {len(all_results)}"
        )

        if not all_results:
            error_msg = f"No files were successfully processed for job {job_id}"
            logger.error(error_msg)
            logger.error(f"Total files downloaded: {len(downloaded_files)}")
            logger.error(f"Files processed: {total_processed}")
            logger.error(f"All results count: {len(all_results)}")
            update_job_progress(
                job_id,
                "failed",
                error=error_msg,
                finished_at=datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            )
            return

        # Create combined zip report
        combined_zip_path = os.path.join(output_dir, f"pmultiqc_report_{job_id}.zip")

        if not create_zip_report(output_dir, combined_zip_path):
            logger.error(f"Failed to create combined zip report for job {job_id}")
            update_job_progress(
                job_id,
                "failed",
                error="Failed to create zip report",
                finished_at=datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            )
            return

        # Update job status to completed immediately after pmultiqc succeeds
        logger.info(f"Updating PRIDE job {job_id} status to completed after pmultiqc success")
        try:
            update_job_progress(
                job_id,
                "completed",
                100,
                processing_stage="Analysis completed successfully!",
                finished_at=datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            )
            logger.info(f"PRIDE job {job_id} status updated to completed successfully")
        except Exception as e:
            logger.error(f"Failed to update job {job_id} status to completed: {e}")
            logger.error(f"Traceback: {traceback.format_exc()}")
            # Continue with HTML report copying even if status update fails

        logger.info(f"Starting HTML report processing for PRIDE job {job_id}")

        # Copy HTML report to public directory for online viewing
        try:
            logger.info(
                f"Starting HTML report copying for PRIDE job {job_id}, output_dir: {output_dir}"
            )
            # Copy HTML reports without signal-based timeout (signals don't work in background threads)
            html_urls = copy_html_report_for_online_viewing(output_dir, job_id, accession)
            logger.info(f"Generated html_urls for PRIDE job {job_id}: {html_urls}")
            if html_urls is None:
                logger.warning(f"HTML report copying returned None for PRIDE job {job_id}")

        except Exception as e:
            logger.error(f"Failed to copy HTML reports for PRIDE job {job_id}: {e}")
            logger.error(f"Traceback: {traceback.format_exc()}")
            html_urls = None

        logger.info(
            f"Completed HTML report processing for PRIDE job {job_id}, proceeding to final steps"
        )

        # Save job data to disk for persistence
        job_data = {
            "job_id": job_id,
            "status": "completed",
            "progress": 100,
            "accession": accession,
            "submission_type": "complete" if is_complete else "standard",
            "files_processed": total_processed,
            "total_files": len(downloaded_files),
            "results": all_results,
            "finished_at": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "download_url": f'{BASE_URL.rstrip("/")}/download-report/{job_id}',
        }

        # Add HTML report URLs if available
        if html_urls:
            job_data["html_report_urls"] = html_urls
            # Also update the job in Redis with HTML URLs
            try:
                update_job_progress(job_id, "completed", 100, html_report_urls=html_urls)
            except Exception as e:
                logger.error(f"Failed to update job {job_id} with HTML URLs in Redis: {e}")

        save_job_data_to_disk(job_id, job_data, output_dir)

        # Save to database with additional fallback mechanisms
        try:
            save_job_to_db(job_id, job_data)
            logger.info(f"Successfully saved PRIDE job {job_id} to database")
        except Exception as e:
            logger.error(f"Failed to save job {job_id} to database: {e}")
            # Create a fallback completed job file
            fallback_file = os.path.join(output_dir, "job_completed.json")
            try:
                with open(fallback_file, "w") as f:
                    json.dump(job_data, f, indent=2)
                logger.info(f"Created fallback job completion file: {fallback_file}")

                # Also try to update Redis directly as final fallback
                try:
                    update_job_progress(
                        job_id,
                        "completed",
                        100,
                        processing_stage="Analysis completed successfully!",
                        finished_at=datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                    )
                    logger.info(f"Final fallback: Updated PRIDE job {job_id} status in Redis")
                except Exception as redis_e:
                    logger.error(f"Final fallback Redis update failed for job {job_id}: {redis_e}")

            except Exception as fallback_e:
                logger.error(f"Failed to create fallback job file: {fallback_e}")

        logger.info(
            f"Successfully completed PRIDE job {job_id}: {total_processed}/{len(downloaded_files)} files processed"
        )
        logger.info(f"process_pride_job_async finished successfully for job {job_id}")

        # Final status verification and logging
        try:
            final_job_status = get_job_from_db(job_id)
            if final_job_status:
                logger.info(
                    f"Final job status verification for {job_id}: {final_job_status.get('status', 'unknown')}"
                )
            else:
                logger.warning(
                    f"Could not verify final status for job {job_id} - job not found in database"
                )
        except Exception as verify_e:
            logger.warning(f"Could not verify final status for job {job_id}: {verify_e}")

    except Exception as e:
        logger.error(f"Error in PRIDE async processing for job {job_id}: {e}")
        logger.error(f"process_pride_job_async failed for job {job_id}: {traceback.format_exc()}")
        update_job_progress(
            job_id,
            "failed",
            error=str(e),
            finished_at=datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        )


# Files with completion time exceeding 24 hours will be cleaned up.
EXPIRATION_SECONDS = 24 * 3600
CLEANUP_JOB_SECONDS = 3600


def allowed_file(filename):
    """Check if the uploaded file has an allowed extension."""
    return "." in filename and filename.rsplit(".", 1)[1].lower() in ALLOWED_EXTENSIONS


async def validate_and_save_upload(file: UploadFile, upload_dir: str, job_id: str) -> tuple[str, int]:
    """
    Validate and save an uploaded file with security checks.

    Args:
        file: The uploaded file
        upload_dir: Directory to save the file
        job_id: Job ID for error reporting

    Returns:
        tuple: (file_path, file_size)

    Raises:
        HTTPException: If validation fails
    """
    # Check file extension
    if not file.filename or not file.filename.lower().endswith(".zip"):
        raise HTTPException(status_code=400, detail="Only ZIP files are allowed")

    # Security: Validate filename to prevent path traversal
    if ".." in file.filename or "/" in file.filename or "\\" in file.filename:
        raise HTTPException(status_code=400, detail="Invalid filename")

    filename = file.filename
    zip_path = os.path.join(upload_dir, filename)

    # Security: Check file size while saving
    file_size = 0
    with open(zip_path, "wb") as buffer:
        chunk_size = 8192
        while chunk := await file.read(chunk_size):
            file_size += len(chunk)
            # Security: Enforce maximum file size
            if file_size > MAX_FILE_SIZE:
                # Clean up partial file
                buffer.close()
                os.remove(zip_path)
                shutil.rmtree(upload_dir, ignore_errors=True)
                raise HTTPException(
                    status_code=413,
                    detail=f"File too large. Maximum size is {MAX_FILE_SIZE / (1024**3):.1f} GB"
                )
            buffer.write(chunk)

    return zip_path, file_size


def validate_and_extract_zip(zip_path: str, extract_path: str, file_size: int):
    """
    Validate and extract a ZIP file with security checks.

    Args:
        zip_path: Path to the ZIP file
        extract_path: Directory to extract to
        file_size: Size of the ZIP file

    Raises:
        HTTPException: If validation fails
    """
    try:
        with zipfile.ZipFile(zip_path, "r") as zip_ref:
            # Security: Check number of files in zip
            file_list = zip_ref.namelist()
            if len(file_list) > MAX_UPLOAD_FILES:
                raise HTTPException(
                    status_code=400,
                    detail=f"Too many files in ZIP archive. Maximum is {MAX_UPLOAD_FILES} files"
                )

            # Security: Check for zip bombs and path traversal in zip entries
            total_uncompressed_size = 0
            for file_info in zip_ref.infolist():
                # Check for path traversal in zip entry names
                if ".." in file_info.filename or file_info.filename.startswith("/"):
                    raise HTTPException(
                        status_code=400,
                        detail="ZIP file contains invalid file paths"
                    )
                total_uncompressed_size += file_info.file_size

            # Check for zip bombs (compression ratio > 100:1)
            compression_ratio = total_uncompressed_size / file_size if file_size > 0 else 0
            if compression_ratio > 100:
                raise HTTPException(
                    status_code=400,
                    detail="Suspicious ZIP file detected (possible zip bomb)"
                )

            zip_ref.extractall(extract_path)
    except zipfile.BadZipFile:
        raise HTTPException(status_code=400, detail="Invalid ZIP file")
    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Error extracting ZIP file: {e}")
        raise HTTPException(status_code=400, detail=f"Error extracting ZIP file: {str(e)}")


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
        if "multiqc_config.yml" in files and any(
            any(substr in f for substr in quantms_files) for f in files
        ):
            return "quantms", file_paths["multiqc_config.yml"]

        # Check for MaxQuant files
        maxquant_files = ["evidence.txt", "msms.txt", "peptides.txt", "proteinGroups.txt"]
        if any(f in files for f in maxquant_files):
            return "maxquant", None

        # Check for DIANN files - more comprehensive detection
        diann_files = ["report.tsv", "report.parquet", "diann_report.tsv", "diann_report.parquet"]
        diann_patterns = ["*report.tsv", "*report.parquet", "*diann_report*"]
        diann_indicators = ["diann", "dia", "dda"]

        # Check for exact DIANN files
        has_diann_files = any(f in files for f in diann_files)

        # Check for DIANN files with patterns
        has_diann_patterns = any(
            any(f.endswith(pattern.replace("*", "")) for pattern in diann_patterns) for f in files
        )

        # Check for DIA indicators in filenames
        has_diann_indicators = any(
            any(indicator in f.lower() for indicator in diann_indicators) for f in files
        )

        if has_diann_files or has_diann_patterns or has_diann_indicators:
            logger.info(
                f"Detected DIANN files: exact={has_diann_files}, patterns={has_diann_patterns}, indicators={has_diann_indicators}"
            )
            logger.info(
                f"Files found: {[f for f in files if any(indicator in f.lower() for indicator in diann_indicators)]}"
            )
            return "diann", None

        # Check for mzIdentML files (CPMPLETE submissions)
        has_mzid = any(f.endswith(".mzid") for f in files)
        has_peak_lists = any(
            f.lower().endswith((".mgf", ".mzml", ".mgf.gz", ".mzml.gz", ".mgf.zip", ".mzml.zip"))
            for f in files
        )

        if has_mzid:
            if has_peak_lists:
                logger.info("Detected COMPLETE submission with mzIdentML and peak list files")
            else:
                logger.info("Detected mzIdentML files without peak lists")
            return "mzidentml", None

        return "unknown", None

    except Exception as e:
        logger.error(f"Error detecting input type: {e}")
        return "unknown", None


def run_pmultiqc_with_progress(
    input_path: str, output_path: str, input_type: str, pmultiqc_config: str, job_id: str
) -> Dict[str, Any]:
    """
    Run pmultiqc with real-time progress updates stored in job_status_dict.
    """
    try:
        # Security: Validate input_type to prevent command injection
        allowed_input_types = ["maxquant", "quantms", "diann", "mzidentml"]
        if input_type not in allowed_input_types:
            logger.error(f"Invalid input_type: {input_type}")
            return {
                "success": False,
                "message": f"Invalid input type: {input_type}",
                "output": [],
                "errors": []
            }

        # Security: Validate paths to ensure they don't contain injection attempts
        if not os.path.exists(input_path):
            logger.error(f"Input path does not exist: {input_path}")
            return {
                "success": False,
                "message": "Input path does not exist",
                "output": [],
                "errors": []
            }

        # Set up MultiQC arguments based on input type
        args = ["multiqc", input_path, "-o", output_path, "--force"]

        # Add type-specific arguments
        if input_type == "maxquant":
            args.extend(["--maxquant-plugin", "--ignore", "summary.txt"])
        elif input_type == "quantms":
            if pmultiqc_config:
                # Security: Validate config file path
                if not os.path.exists(pmultiqc_config):
                    logger.error(f"Config file does not exist: {pmultiqc_config}")
                    return {
                        "success": False,
                        "message": "Config file does not exist",
                        "output": [],
                        "errors": []
                    }
                args.extend(["--quantms-plugin", "--config", pmultiqc_config])
            else:
                logger.error("The function 'run_pmultiqc_with_progress' is missing the parameter: pmultiqc_config")
        elif input_type == "diann":
            # DIANN files are handled automatically by pmultiqc
            # Add memory optimization arguments for large files
            args.extend(["--diann-plugin", "--no-megaqc-upload", "--verbose"])
        elif input_type == "mzidentml":
            args.extend(["--mzid-plugin"])

        # Run MultiQC with pmultiqc plugin
        logger.info(f"Running pmultiqc with args: {args}")
        logger.info(f"Detected input type: {input_type}")

        # Test if pmultiqc plugin is available
        try:
            import pmultiqc

            logger.info(f"pmultiqc plugin version: {pmultiqc.__version__}")
        except ImportError as e:
            logger.error(f"pmultiqc plugin not available: {e}")
        except Exception as e:
            logger.error(f"Error checking pmultiqc plugin: {e}")

        # Add timeout for large datasets (30 minutes)

        # Check if we have a large file and set timeout accordingly
        large_file_detected = False
        timeout_seconds = 30 * 60  # 30 minutes default

        if input_type == "diann":
            # Check if we have a large file
            input_dir = args[1]
            try:
                for root, _, files in os.walk(input_dir):
                    for file in files:
                        if file.endswith(".tsv"):
                            file_path = os.path.join(root, file)
                            file_size = os.path.getsize(file_path)
                            if file_size > 500 * 1024 * 1024:  # > 500MB
                                large_file_detected = True
                                break
                    if large_file_detected:
                        break
            except (OSError, IOError) as e:
                logger.warning(f"Could not check file sizes in {input_dir}: {e}")

            if large_file_detected:
                logger.info("Large file detected, setting 30-minute timeout")
                timeout_seconds = 30 * 60  # 30 minutes for large files

        # Set environment variables to fix matplotlib and other issues
        import tempfile

        matplotlib_dir = os.environ.get(
            "MPLCONFIGDIR", os.path.join(tempfile.gettempdir(), "matplotlib")
        )
        os.environ["MPLCONFIGDIR"] = matplotlib_dir
        os.makedirs(matplotlib_dir, exist_ok=True)

        # Set output directory environment variable to help pmultiqc plugin
        output_dir = (
            args[args.index("-o") + 1]
            if "-o" in args
            else os.path.join(tempfile.gettempdir(), "multiqc_output")
        )
        os.environ["MULTIQC_OUTPUT_DIR"] = output_dir
        os.makedirs(output_dir, exist_ok=True)

        # Set MultiQC config to use the output directory
        os.environ["MULTIQC_CONFIG"] = f"output_dir: {output_dir}"

        # Preprocess DIANN data to handle NaN values in is_contaminant column
        # Run preprocessing for both 'diann' and 'quantms' input types since they both use DIANN reports
        # Also run as fallback if input type detection failed but report.tsv files exist
        should_preprocess = input_type in ["diann", "quantms"]

        # Fallback: check if report.tsv files exist even if input type detection failed
        if not should_preprocess:
            try:
                import glob

                input_dir = args[1]  # Second argument is the input directory (first is 'multiqc')
                if os.path.exists(input_dir):
                    report_files = glob.glob(f"{input_dir}/**/report.tsv", recursive=True)
                    if report_files:
                        logger.info(
                            f"Found report.tsv files despite input_type={input_type}, running preprocessing as fallback"
                        )
                        should_preprocess = True
            except (OSError, IOError, ImportError) as e:
                logger.debug(f"Could not check for report files in fallback detection: {e}")

        if should_preprocess:
            logger.info(f"Preprocessing {input_type} data to handle NaN values...")
            try:
                import pandas as pd
                import glob

                # Find all report.tsv files in the input directory
                input_dir = args[1]  # Second argument is the input directory (first is 'multiqc')
                logger.info(f"Looking for report files in: {input_dir}")

                # Check if input directory exists
                if not os.path.exists(input_dir):
                    logger.warning(f"Input directory does not exist: {input_dir}")
                    # Continue with MultiQC even if preprocessing fails

                # List all files in the input directory
                all_files = []
                for root, _, files in os.walk(input_dir):
                    for file in files:
                        all_files.append(os.path.join(root, file))
                logger.info(f"All files in input directory: {all_files}")

                report_files = glob.glob(f"{input_dir}/**/report.tsv", recursive=True)
                logger.info(f"Found report files: {report_files}")

                for report_file in report_files:
                    logger.info(f"Processing report file: {report_file}")
                    try:
                        # Check file size first
                        file_size = os.path.getsize(report_file)
                        file_size_mb = file_size / (1024 * 1024)
                        logger.info(f"Report file size: {file_size_mb:.2f} MB")

                        # Warn if file is very large
                        if file_size_mb > 1000:  # > 1GB
                            logger.warning(
                                f"Large dataset detected: {file_size_mb:.2f} MB. This may take a long time to process."
                            )

                        # Read the report file with chunking for large files
                        if file_size_mb > 500:  # > 500MB, use chunking
                            logger.info("Large file detected, using chunked reading...")
                            chunk_size = 100000  # 100k rows per chunk
                            chunks_processed = 0
                            total_rows = 0

                            # First, count total rows
                            logger.info("Counting total rows in file...")
                            with open(report_file, "r") as f:
                                total_rows = sum(1 for line in f) - 1  # Subtract header
                            logger.info(f"Total rows in file: {total_rows:,}")

                            # Process in chunks
                            for chunk_df in pd.read_csv(
                                report_file, sep="\t", chunksize=chunk_size
                            ):
                                chunks_processed += 1
                                rows_in_chunk = len(chunk_df)
                                total_processed = chunks_processed * chunk_size
                                progress_pct = min(100, (total_processed / total_rows) * 100)

                                logger.info(
                                    f"Processing chunk {chunks_processed}: {rows_in_chunk:,} rows (Progress: {progress_pct:.1f}%)"
                                )

                                # Check if is_contaminant column exists and has NaN values
                                if "is_contaminant" in chunk_df.columns:
                                    nan_count = chunk_df["is_contaminant"].isna().sum()
                                    if nan_count > 0:
                                        logger.info(
                                            f"Found {nan_count:,} NaN values in chunk {chunks_processed}"
                                        )
                                        chunk_df["is_contaminant"] = chunk_df[
                                            "is_contaminant"
                                        ].fillna(False)

                                # Write chunk to temporary file
                                temp_file = f"{report_file}.chunk_{chunks_processed}"
                                chunk_df.to_csv(temp_file, sep="\t", index=False)

                                if chunks_processed % 10 == 0:  # Log every 10 chunks
                                    logger.info(
                                        f"Processed {chunks_processed} chunks, {total_processed:,} rows"
                                    )

                            # Combine chunks back into original file
                            logger.info("Combining processed chunks back into original file...")
                            first_chunk = True
                            with open(report_file, "w") as outfile:
                                for i in range(1, chunks_processed + 1):
                                    temp_file = f"{report_file}.chunk_{i}"
                                    with open(temp_file, "r") as infile:
                                        if first_chunk:
                                            outfile.write(
                                                infile.read()
                                            )  # Write header and first chunk
                                            first_chunk = False
                                        else:
                                            next(infile)  # Skip header
                                            outfile.write(infile.read())  # Write data
                                    os.remove(temp_file)  # Clean up temp file

                            logger.info(
                                f"Successfully processed {total_rows:,} rows in {chunks_processed} chunks"
                            )
                        else:
                            # Regular processing for smaller files
                            logger.info("Reading file normally...")
                            df = pd.read_csv(report_file, sep="\t")
                            logger.info(
                                f"Loaded report file with {len(df):,} rows and columns: {list(df.columns)}"
                            )

                            # Check if is_contaminant column exists and has NaN values
                            if "is_contaminant" in df.columns:
                                nan_count = df["is_contaminant"].isna().sum()
                                logger.info(
                                    f"Found {nan_count:,} NaN values in is_contaminant column"
                                )
                                if nan_count > 0:
                                    logger.info(
                                        f"Filling {nan_count:,} NaN values in is_contaminant column with False"
                                    )
                                    df["is_contaminant"] = df["is_contaminant"].fillna(False)

                                    # Write the cleaned data back
                                    df.to_csv(report_file, sep="\t", index=False)
                                    logger.info(f"Cleaned report file saved: {report_file}")
                            else:
                                logger.info("No is_contaminant column found in report file")

                    except Exception as e:
                        logger.warning(f"Could not preprocess report file {report_file}: {e}")

            except Exception as e:
                logger.warning(f"Error during DIANN data preprocessing: {e}")

        # Initialize console output in job status
        update_job_progress(
            job_id,
            "running_pmultiqc",
            0,
            console_output=[],
            console_errors=[],
            processing_stage="Starting pmultiqc analysis...",
        )

        # Run multiqc as a subprocess with real-time output
        try:
            # Handle timeout for large files
            if input_type == "diann" and large_file_detected:
                logger.info("Starting MultiQC with timeout handling for large file")

            # Run the multiqc command with real-time output
            env = os.environ.copy()
            process = subprocess.Popen(
                args,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                bufsize=1,
                universal_newlines=True,
                env=env,
            )

            # Stream output in real-time
            output_lines = []
            error_lines = []

            # Set a timeout for the process
            start_time = time.time()
            logger.info(f"Starting MultiQC with {timeout_seconds} second timeout")

            while True:
                # Check for timeout
                if time.time() - start_time > timeout_seconds:
                    logger.error(
                        f"MultiQC process timed out after {timeout_seconds} seconds for job {job_id}"
                    )
                    process.terminate()
                    return {
                        "success": False,
                        "message": f"MultiQC process timed out after {timeout_seconds} seconds",
                        "output": output_lines,
                        "errors": error_lines
                        + [f"Process timed out after {timeout_seconds} seconds"],
                    }

                output = process.stdout.readline()
                error = process.stderr.readline()

                if output:
                    output_line = output.strip()
                    output_lines.append(output_line)
                    update_job_progress(
                        job_id,
                        "running_pmultiqc",
                        0,
                        console_output=output_lines,
                        console_errors=error_lines,
                        processing_stage="Running pmultiqc...",
                    )
                    logger.info(f"MultiQC: {output_line}")

                if error:
                    error_line = error.strip()
                    error_lines.append(error_line)
                    update_job_progress(
                        job_id,
                        "running_pmultiqc",
                        0,
                        console_output=output_lines,
                        console_errors=error_lines,
                        processing_stage="Running pmultiqc...",
                    )
                    logger.warning(f"MultiQC: {error_line}")

                # Check if process has finished
                if process.poll() is not None:
                    logger.info(f"MultiQC process finished with return code: {process.returncode}")
                    break

            # Get any remaining output
            remaining_output, remaining_error = process.communicate()
            if remaining_output:
                for line in remaining_output.strip().split("\n"):
                    if line.strip():
                        output_lines.append(line.strip())
                        update_job_progress(
                            job_id,
                            "running_pmultiqc",
                            0,
                            console_output=output_lines,
                            console_errors=error_lines,
                            processing_stage="Running pmultiqc...",
                        )
                        logger.info(f"MultiQC remaining output: {line.strip()}")
            if remaining_error:
                for line in remaining_error.strip().split("\n"):
                    if line.strip():
                        error_lines.append(line.strip())
                        update_job_progress(
                            job_id,
                            "running_pmultiqc",
                            0,
                            console_output=output_lines,
                            console_errors=error_lines,
                            processing_stage="Running pmultiqc...",
                        )
                        logger.warning(f"MultiQC remaining error: {line.strip()}")

            if process.returncode == 0:
                logger.info(f"MultiQC completed successfully for job {job_id}")
                return {
                    "success": True,
                    "message": "Report generated successfully",
                    "output": output_lines,
                    "errors": error_lines,
                }
            else:
                logger.error(
                    f"MultiQC failed with return code {process.returncode} for job {job_id}"
                )
                return {
                    "success": False,
                    "message": f"MultiQC failed with return code {process.returncode}",
                    "output": output_lines,
                    "errors": error_lines,
                }

        except subprocess.CalledProcessError as e:
            logger.error(f"MultiQC failed with return code {e.returncode} for job {job_id}")
            logger.error(f"MultiQC stdout: {e.stdout}")
            logger.error(f"MultiQC stderr: {e.stderr}")
            return {
                "success": False,
                "message": f"MultiQC failed: {e.stderr}",
                "output": e.stdout.split("\n") if e.stdout else [],
                "errors": e.stderr.split("\n") if e.stderr else [],
            }

    except Exception as e:
        logger.error(f"Error running pmultiqc: {e}")
        return {"success": False, "message": str(e), "output": [], "errors": []}


def create_zip_report(output_path: str, zip_path: str) -> bool:
    """
    Create a zip file containing the pmultiqc report.

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
        for root, _, files in os.walk(output_path):
            for file in files:
                file_path = os.path.join(root, file)
                # Skip the zip file we're creating to avoid infinite loop
                if file_path != zip_path:
                    files_found.append(file_path)

        logger.info(f"Found {len(files_found)} files to zip")

        if not files_found:
            logger.error("No files found to zip")
            return False

        # Create zip file with better error handling
        with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as zipf:
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
            logger.info(
                f"Successfully created zip file: {zip_path} ({os.path.getsize(zip_path)} bytes)"
            )
            return True
        else:
            logger.error(f"Zip file was not created or is empty: {zip_path}")
            return False

    except Exception as e:
        logger.error(f"Error creating zip report: {e}")
        return False


def process_job_async(
    job_id: str, extract_path: str, output_dir: str, input_type: str, quantms_config: str
):
    """Process a job asynchronously."""
    try:
        # Update job status to running
        update_job_progress(
            job_id,
            "running_pmultiqc",
            0,
            console_output=[],
            console_errors=[],
            processing_stage="Starting pmultiqc analysis...",
        )

        # Run pmultiqc with real-time progress
        logger.info(f"Starting async processing for job {job_id}")
        result = run_pmultiqc_with_progress(
            extract_path, output_dir, input_type, quantms_config, job_id
        )

        if not result["success"]:
            update_job_progress(
                job_id,
                "failed",
                error=result["message"],
                finished_at=datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            )
            return

        # Update job status to creating report
        update_job_progress(
            job_id, "creating_report", 75, processing_stage="Creating zip report..."
        )

        # Create zip report
        zip_report_path = os.path.join(output_dir, f"pmultiqc_report_{job_id}.zip")
        if not create_zip_report(output_dir, zip_report_path):
            update_job_progress(
                job_id,
                "failed",
                error="Failed to create zip report",
                finished_at=datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            )
            return

        # Update job status to completed immediately after pmultiqc succeeds
        logger.info(f"Updating async job {job_id} status to completed after pmultiqc success")
        try:
            update_job_progress(
                job_id,
                "completed",
                100,
                processing_stage="Analysis completed successfully!",
                finished_at=datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            )
            logger.info(f"Async job {job_id} status updated to completed successfully")
        except Exception as e:
            logger.error(f"Failed to update async job {job_id} status to completed: {e}")
            logger.error(f"Traceback: {traceback.format_exc()}")
            # Continue with HTML report copying even if status update fails

        # Copy HTML report to public directory for online viewing
        try:
            logger.info(
                f"Starting HTML report copying for async job {job_id}, output_dir: {output_dir}"
            )
            html_urls = copy_html_report_for_online_viewing(output_dir, job_id)
            logger.info(f"Generated html_urls for async job {job_id}: {html_urls}")
            if html_urls is None:
                logger.warning(f"HTML report copying returned None for async job {job_id}")
        except Exception as e:
            logger.error(f"Failed to copy HTML reports for async job {job_id}: {e}")
            logger.error(f"Traceback: {traceback.format_exc()}")
            html_urls = None

        # Update job to completed
        final_job_data = {
            "status": "completed",
            "progress": 100,
            "input_type": input_type,
            "processing_stage": "Analysis completed successfully!",
            "finished_at": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "download_url": f'{BASE_URL.rstrip("/")}/download-report/{job_id}',
        }

        # Add HTML report URLs if available
        if html_urls:
            final_job_data["html_report_urls"] = html_urls
            # Also update the job in Redis with HTML URLs
            try:
                update_job_progress(job_id, "completed", 100, html_report_urls=html_urls)
            except Exception as e:
                logger.error(f"Failed to update async job {job_id} with HTML URLs in Redis: {e}")

        # Save final job data to database
        try:
            save_job_to_db(job_id, final_job_data)
        except Exception as e:
            logger.error(f"Failed to save async job {job_id} to database: {e}")
            # Create a fallback completed job file
            fallback_file = os.path.join(output_dir, "job_completed.json")
            try:
                with open(fallback_file, "w") as f:
                    json.dump(final_job_data, f, indent=2)
                logger.info(f"Created fallback async job completion file: {fallback_file}")
            except Exception as fallback_e:
                logger.error(f"Failed to create fallback async job file: {fallback_e}")

        logger.info(f"Successfully completed async job {job_id}")

    except Exception as e:
        logger.error(f"Error in async processing for job {job_id}: {e}")
        update_job_progress(
            job_id,
            "failed",
            error=str(e),
            finished_at=datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        )


# API Routes
@app.get("/", response_class=HTMLResponse)
async def index(request: Request):
    """Main page of index.html, with upload interface and PRIDE submission page

    Parameters:
        request: Request object

    Returns:
        HTMLResponse: index.html template with upload interface and PRIDE submission page
    """
    return templates.TemplateResponse(
        "index.html", {"request": request, "config": {"BASE_URL": BASE_URL, "PRIDE_BUTTON_VISIBLE": PRIDE_BUTTON_VISIBLE}}
    )


@app.get("/submit")
async def submit_pride(request: Request):
    """Direct submission page for PRIDE datasets via URL parameter.

    Parameters:
        request: Request object

    Returns:
        HTMLResponse: submit.html template with PRIDE submission page
    """
    accession = request.query_params.get("accession", "").strip().upper()

    if not accession:
        # No accession provided, show the main page
        return templates.TemplateResponse(
            "index.html", {"request": request, "config": {"BASE_URL": BASE_URL, "PRIDE_BUTTON_VISIBLE": PRIDE_BUTTON_VISIBLE}}
        )

    # Security: Validate accession format (PXD followed by 6 digits)
    if not re.match(r'^PXD\d{6}$', accession):
        return templates.TemplateResponse(
            "index.html",
            {
                "request": request,
                "config": {"BASE_URL": BASE_URL, "PRIDE_BUTTON_VISIBLE": PRIDE_BUTTON_VISIBLE},
                "error": f"Invalid PRIDE accession format: {accession}. Should be PXD followed by 6 digits (e.g., PXD012345).",
            },
        )

    # Render the submission page with the accession pre-filled
    return templates.TemplateResponse(
        "submit.html",
        {"request": request, "accession": accession, "config": {"BASE_URL": BASE_URL, "PRIDE_BUTTON_VISIBLE": PRIDE_BUTTON_VISIBLE}},
    )


@app.get("/results")
async def view_results(request: Request):
    """View results page for a specific job.

    Parameters:
        request: Request object

    Returns:
        HTMLResponse: results.html template with results page
    """
    try:
        # Get job_id from query parameter
        job_id = request.query_params.get("job", "").strip()

        if not job_id:
            return templates.TemplateResponse(
                "index.html",
                {
                    "request": request,
                    "config": {"BASE_URL": BASE_URL, "PRIDE_BUTTON_VISIBLE": PRIDE_BUTTON_VISIBLE},
                    "error": f"No job ID provided. Please provide a job ID using ?job={job_id}",
                },
            )

        # Validate job_id format
        try:
            uuid.UUID(job_id)
        except ValueError:
            return templates.TemplateResponse(
                "index.html",
                {
                    "request": request,
                    "config": {"BASE_URL": BASE_URL, "PRIDE_BUTTON_VISIBLE": PRIDE_BUTTON_VISIBLE},
                    "error": f"Invalid job ID format. {job_id}",
                },
            )

        logger.info(f"Results page request for {job_id}")

        # Try to get job from database first
        job_data = get_job_from_db(job_id)

        if job_data:
            logger.info(f"Found job {job_id} in database with status: {job_data.get('status')}")
        else:
            logger.warning(f"Job {job_id} NOT found in database")

            # Check if job is completed (for backward compatibility with old jobs)
            output_dir = os.path.join(OUTPUT_FOLDER, job_id)
            zip_report_path = os.path.join(output_dir, f"pmultiqc_report_{job_id}.zip")

            if os.path.exists(zip_report_path):
                logger.info(f"Job {job_id} not in database but zip file exists")

                # Try to recover job data from output directory
                console_output = []
                console_errors = []
                input_type = None

                # Look for console output files
                console_output_file = os.path.join(output_dir, "console_output.txt")
                console_errors_file = os.path.join(output_dir, "console_errors.txt")
                job_info_file = os.path.join(output_dir, "job_info.json")

                if os.path.exists(console_output_file):
                    try:
                        with open(console_output_file, "r", encoding="utf-8") as f:
                            console_output = [
                                line.strip() for line in f.readlines() if line.strip()
                            ]
                        logger.info(
                            f"Recovered {len(console_output)} console output lines for job {job_id}"
                        )
                    except Exception as e:
                        logger.error(f"Error reading console output file: {e}")

                if os.path.exists(console_errors_file):
                    try:
                        with open(console_errors_file, "r", encoding="utf-8") as f:
                            console_errors = [
                                line.strip() for line in f.readlines() if line.strip()
                            ]
                        logger.info(
                            f"Recovered {len(console_errors)} console error lines for job {job_id}"
                        )
                    except Exception as e:
                        logger.error(f"Error reading console errors file: {e}")

                if os.path.exists(job_info_file):
                    try:
                        with open(job_info_file, "r", encoding="utf-8") as f:
                            job_info = json.load(f)
                            input_type = job_info.get("input_type")
                        logger.info(
                            f"Recovered job info for job {job_id}: input_type={input_type}"
                        )
                    except Exception as e:
                        logger.error(f"Error reading job info file: {e}")

                # Create job_data from recovered information
                job_data = {
                    "job_id": job_id,
                    "status": "completed",
                    "progress": 100,
                    "input_type": input_type,
                    "accession": None,
                    "started_at": None,
                    "finished_at": None,
                    "error": None,
                    "files_processed": None,
                    "total_files": None,
                    "console_output": console_output,
                    "console_errors": console_errors,
                    "html_report_urls": {},  # Initialize as empty dict
                    "download_url": f'{BASE_URL.rstrip("/")}/download-report/{job_id}',
                }

                # Save recovered job to database for future requests
                save_job_to_db(job_id, job_data)

                logger.info("Returning fallback response for completed job " + job_id)
            else:
                # Job not found - show results page with "not found" status
                job_data = {
                    "job_id": job_id,
                    "status": "not_found",
                    "progress": 0,
                    "input_type": None,
                    "accession": None,
                    "started_at": None,
                    "finished_at": None,
                    "error": f"Job {job_id} not found. This job may have expired, been deleted, or never existed.",
                    "files_processed": None,
                    "total_files": None,
                    "console_output": [],
                    "console_errors": [],
                    "html_report_urls": {},
                    "download_url": None,
                }

        # Ensure download_url is present in job_data
        if "download_url" not in job_data:
            job_data["download_url"] = f'{BASE_URL.rstrip("/")}/download-report/{job_id}'

        return templates.TemplateResponse(
            "results.html",
            {
                "request": request,
                "job_id": job_id,
                "job_data": job_data,
                "config": {"BASE_URL": BASE_URL, "PRIDE_BUTTON_VISIBLE": PRIDE_BUTTON_VISIBLE},
            },
        )

    except Exception as e:
        logger.error(f"Error in view_results for job {job_id}: {e}")
        # Show results page with error status instead of redirecting to index
        job_data = {
            "job_id": job_id if "job_id" in locals() else "unknown",
            "status": "error",
            "progress": 0,
            "input_type": None,
            "accession": None,
            "started_at": None,
            "finished_at": None,
            "error": f"Error loading job: {str(e)}",
            "files_processed": None,
            "total_files": None,
            "console_output": [],
            "console_errors": [],
            "html_report_urls": {},
            "download_url": None,
        }
        return templates.TemplateResponse(
            "results.html",
            {
                "request": request,
                "job_id": job_data["job_id"],
                "job_data": job_data,
                "config": {"BASE_URL": BASE_URL, "PRIDE_BUTTON_VISIBLE": PRIDE_BUTTON_VISIBLE},
            },
        )


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
        if (
            not report_hash
            or len(report_hash) != 12
            or not all(c in "0123456789abcdef" for c in report_hash)
        ):
            logger.error(f"Invalid report hash format: {report_hash}")
            raise HTTPException(status_code=400, detail="Invalid report hash")

        # Security: Validate filename to prevent path traversal
        # Reject any filename containing path traversal patterns
        if ".." in filename or filename.startswith("/") or "\\" in filename:
            logger.error(f"Potential path traversal attempt detected: {filename}")
            raise HTTPException(status_code=400, detail="Invalid filename")

        # Construct file path
        file_path = os.path.join(HTML_REPORTS_FOLDER, report_hash, filename)
        logger.info(f"Looking for file at: {file_path}")

        # List contents of the report directory for debugging
        report_dir = os.path.join(HTML_REPORTS_FOLDER, report_hash)
        if os.path.exists(report_dir):
            logger.info(f"Report directory exists. Contents:")
            for root, _, files in os.walk(report_dir):
                level = root.replace(report_dir, "").count(os.sep)
                indent = " " * 2 * level
                logger.info(f"{indent}{os.path.basename(root)}/")
                subindent = " " * 2 * (level + 1)
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
        mime_type = "text/html"  # Default for HTML files
        if filename.endswith(".css"):
            mime_type = "text/css"
        elif filename.endswith(".js"):
            mime_type = "application/javascript"
        elif filename.endswith(".png"):
            mime_type = "image/png"
        elif filename.endswith(".jpg") or filename.endswith(".jpeg"):
            mime_type = "image/jpeg"
        elif filename.endswith(".svg"):
            mime_type = "image/svg+xml"
        elif filename.endswith(".ico"):
            mime_type = "image/x-icon"
        elif filename.endswith(".json"):
            mime_type = "application/json"
        elif filename.endswith(".txt"):
            mime_type = "text/plain"
        elif filename.endswith(".parquet"):
            mime_type = "application/octet-stream"
        elif filename.endswith(".log"):
            mime_type = "text/plain"

        # Get file size for logging
        file_size = os.path.getsize(file_path)
        logger.info(f"Serving file: {file_path} ({file_size} bytes) with MIME type: {mime_type}")

        # Log timing for large files
        if file_size > 1024 * 1024:  # Files larger than 1MB
            logger.warning(
                f"Large file being served: {filename} ({file_size / (1024*1024):.1f} MB)"
            )

        # Add cache headers for static assets
        if filename.endswith((".css", ".js", ".png", ".jpg", ".jpeg", ".svg", ".ico")):
            response = FileResponse(
                path=file_path,
                media_type=mime_type,
                headers={"Cache-Control": "max-age=3600"},  # Cache for 1 hour
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


@app.get("/health-check")
async def health_check():
    """Health check endpoint."""
    return {
        "schemaVersion": 1,
        "label": "pmultiqc-service",
        "message": "alive",
        "color": "brightgreen",
    }


@app.get("/health")
async def health():
    """Health check endpoint."""
    return {
        "status": "healthy",
        "timestamp": datetime.now().isoformat(),
        "service": "pmultiqc-service",
        "base_url": BASE_URL,
        "endpoints": [
            "/upload-async",
            "/generate-report",
            "/process-pride",
            "/submit-pride",
            "/job-status/{job_id}",
            "/download-report/{job_id}",
            "/reports/{hash}/{filename}",
            "/docs",
            "/swagger.json",
        ],
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
    accession = data.get("accession", "").strip().upper()

    # Security: Validate accession format (PXD followed by 6 digits)
    if not re.match(r'^PXD\d{6}$', accession):
        raise HTTPException(
            status_code=400, detail="Invalid PRIDE accession format. Should be PXD followed by 6 digits (e.g., PXD012345)"
        )

    # Generate unique job ID
    job_id = str(uuid.uuid4())

    # Create output directory
    output_dir = os.path.join(OUTPUT_FOLDER, job_id)
    os.makedirs(output_dir, exist_ok=True)

    # Initialize job status
    update_job_progress(job_id, "queued", 0, accession=accession)

    logger.info(f"Started PRIDE processing for accession {accession} with job ID {job_id}")

    return {
        "job_id": job_id,
        "accession": accession,
        "status": "queued",
        "message": "PRIDE processing started successfully",
        "status_url": f"/job-status/{job_id}",
        "redirect_url": f"/results?job={job_id}",
    }


@app.post("/upload-async")
async def upload_async(file: UploadFile = File(..., alias="files")):
    """
    Upload and process files asynchronously

    Upload a ZIP file containing data files for pmultiqc analysis.
    The job will be processed asynchronously and you can check the status
    using the returned job_id.
    """
    # Generate unique job ID
    job_id = str(uuid.uuid4())

    # Create job directories
    upload_dir = os.path.join(UPLOAD_FOLDER, job_id)
    output_dir = os.path.join(OUTPUT_FOLDER, job_id)
    os.makedirs(upload_dir, exist_ok=True)
    os.makedirs(output_dir, exist_ok=True)

    logger.info(f"Starting async upload for job {job_id}: {file.filename}")

    # Initialize job in database
    initial_job_data = {
        "job_id": job_id,
        "status": "uploading",
        "progress": 0,
        "started_at": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
    }
    save_job_to_db(job_id, initial_job_data)

    # Validate and save uploaded file
    zip_path, file_size = await validate_and_save_upload(file, upload_dir, job_id)
    logger.info(f"Upload completed for job {job_id}: {file_size} bytes")

    # Extract the zip file
    extract_path = os.path.join(upload_dir, "extracted")
    os.makedirs(extract_path, exist_ok=True)

    # Update job status to extracting
    update_job_progress(job_id, "extracting", 25)
    logger.info(f"Extracting ZIP file for job {job_id}")

    validate_and_extract_zip(zip_path, extract_path, file_size)

    # Detect input type
    input_type, quantms_config = detect_input_type(extract_path)
    logger.info(f"Detected input type: {input_type}")

    if input_type == "unknown":
        update_job_progress(
            job_id,
            "failed",
            error="Could not detect input type",
            finished_at=datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        )
        raise HTTPException(
            status_code=400,
            detail="Could not detect input type. Please ensure your ZIP file contains valid data files.",
        )

    # Start async processing
    update_job_progress(job_id, "queued", 50, input_type=input_type)

    # Start background task for processing
    thread = threading.Thread(
        target=process_job_async,
        args=(job_id, extract_path, output_dir, input_type, quantms_config),
    )
    thread.daemon = True
    thread.start()

    return {
        "job_id": job_id,
        "status": "queued",
        "message": "File uploaded successfully. Processing started.",
        "status_url": f"/job-status/{job_id}",
        "redirect_url": f"/results?job={job_id}",
        "file_size": file_size,
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

        # Try to get job from Redis first (for active jobs)
        redis_client = get_redis_client()
        redis_key = f"{REDIS_KEY_PREFIX}{job_id}"
        redis_data_json = redis_client.get(redis_key)

        if redis_data_json:
            try:
                job_data = json.loads(redis_data_json)
                logger.info(f"Found job {job_id} in Redis with status: {job_data.get('status')}")

                # Ensure all required fields have default values
                job_data_with_defaults = {
                    "job_id": job_data.get("job_id", job_id),
                    "status": job_data.get("status", "unknown"),
                    "progress": job_data.get("progress", 0),
                    "input_type": job_data.get("input_type"),
                    "accession": job_data.get("accession"),
                    "started_at": job_data.get("started_at"),
                    "finished_at": job_data.get("finished_at"),
                    "error": job_data.get("error"),
                    "files_processed": job_data.get("files_processed"),
                    "total_files": job_data.get("total_files"),
                    "files_downloaded": job_data.get("files_downloaded"),
                    "processing_stage": job_data.get("processing_stage"),
                    "console_output": job_data.get("console_output", []),
                    "console_errors": job_data.get("console_errors", []),
                    "html_report_urls": job_data.get("html_report_urls", {}),
                    "download_url": job_data.get("download_url"),
                    "download_details": job_data.get("download_details", {}),
                }

                return job_data_with_defaults
            except Exception as e:
                logger.error(f"Error parsing Redis data for job {job_id}: {e}")

        # Try to get job from database (for completed jobs)
        job_data = get_job_from_db(job_id)

        if job_data:
            logger.info(f"Found job {job_id} in database with status: {job_data.get('status')}")

            # Ensure all required fields have default values
            job_data_with_defaults = {
                "job_id": job_data.get("job_id", job_id),
                "status": job_data.get("status", "unknown"),
                "progress": job_data.get("progress", 0),
                "input_type": job_data.get("input_type"),
                "accession": job_data.get("accession"),
                "started_at": job_data.get("started_at"),
                "finished_at": job_data.get("finished_at"),
                "error": job_data.get("error"),
                "files_processed": job_data.get("files_processed"),
                "total_files": job_data.get("total_files"),
                "files_downloaded": job_data.get("files_downloaded"),
                "processing_stage": job_data.get("processing_stage"),
                "console_output": job_data.get("console_output", []),
                "console_errors": job_data.get("console_errors", []),
                "html_report_urls": job_data.get("html_report_urls", {}),
                "download_url": job_data.get("download_url"),
                "download_details": job_data.get("download_details", {}),
            }

            return job_data_with_defaults
        else:
            logger.warning(f"Job {job_id} NOT found in database")

        # Check if job is completed (for backward compatibility with old jobs)
        output_dir = os.path.join(OUTPUT_FOLDER, job_id)
        zip_report_path = os.path.join(output_dir, f"pmultiqc_report_{job_id}.zip")

        if os.path.exists(zip_report_path):
            logger.info(f"Job {job_id} not in database but zip file exists")

            # Try to recover job data from output directory
            console_output = []
            console_errors = []
            input_type = None

            # Look for console output files
            console_output_file = os.path.join(output_dir, "console_output.txt")
            console_errors_file = os.path.join(output_dir, "console_errors.txt")
            job_info_file = os.path.join(output_dir, "job_info.json")

            if os.path.exists(console_output_file):
                try:
                    with open(console_output_file, "r", encoding="utf-8") as f:
                        console_output = [line.strip() for line in f.readlines() if line.strip()]
                    logger.info(
                        f"Recovered {len(console_output)} console output lines for job {job_id}"
                    )
                except Exception as e:
                    logger.error(f"Error reading console output file: {e}")

            if os.path.exists(console_errors_file):
                try:
                    with open(console_errors_file, "r", encoding="utf-8") as f:
                        console_errors = [line.strip() for line in f.readlines() if line.strip()]
                    logger.info(
                        f"Recovered {len(console_errors)} console error lines for job {job_id}"
                    )
                except Exception as e:
                    logger.error(f"Error reading console errors file: {e}")

            if os.path.exists(job_info_file):
                try:
                    with open(job_info_file, "r", encoding="utf-8") as f:
                        job_info = json.load(f)
                        input_type = job_info.get("input_type")
                    logger.info(f"Recovered job info for job {job_id}: input_type={input_type}")
                except Exception as e:
                    logger.error(f"Error reading job info file: {e}")

            # Fallback response with recovered data
            base_url = BASE_URL.rstrip("/")
            response_data = {
                "job_id": job_id,
                "status": "completed",
                "progress": 100,
                "input_type": input_type,
                "accession": None,
                "started_at": None,
                "finished_at": None,
                "error": None,
                "files_processed": None,
                "total_files": None,
                "files_downloaded": None,
                "processing_stage": "Completed",
                "console_output": console_output,
                "console_errors": console_errors,
                "html_report_urls": None,
                "download_url": f"{base_url}/download-report/{job_id}",
                "download_details": {},
            }

            # Save recovered job to database for future requests
            save_job_to_db(job_id, response_data)

            logger.info("Returning fallback response for completed job " + job_id)
            return response_data
        else:
            # Check if job is completed by looking for output files (final fallback)
            output_dir = os.path.join(OUTPUT_FOLDER, job_id)
            zip_report_path = os.path.join(output_dir, f"pmultiqc_report_{job_id}.zip")
            fallback_file = os.path.join(output_dir, "job_completed.json")

            # Check for fallback completion file first
            if os.path.exists(fallback_file):
                try:
                    with open(fallback_file, "r") as f:
                        job_data = json.load(f)
                    logger.info(f"Found fallback completion file for job {job_id}")
                    return job_data
                except Exception as e:
                    logger.error(f"Error reading fallback file for job {job_id}: {e}")

            if os.path.exists(zip_report_path):
                logger.info(
                    f"Job {job_id} not found in Redis or database, but output files exist - job likely completed"
                )

                # Try to copy HTML reports if they don't exist
                try:
                    html_urls = copy_html_report_for_online_viewing(output_dir, job_id)
                    logger.info(f"Generated html_urls for completed job {job_id}: {html_urls}")
                except Exception as e:
                    logger.error(f"Failed to copy HTML reports for completed job {job_id}: {e}")
                    html_urls = None

                # Create a completed job response
                base_url = BASE_URL.rstrip("/")
                response_data = {
                    "job_id": job_id,
                    "status": "completed",
                    "progress": 100,
                    "input_type": None,
                    "accession": None,
                    "started_at": None,
                    "finished_at": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                    "error": None,
                    "files_processed": None,
                    "total_files": None,
                    "files_downloaded": None,
                    "processing_stage": "Analysis completed successfully!",
                    "console_output": [],
                    "console_errors": [],
                    "html_report_urls": html_urls or {},
                    "download_url": f"{base_url}/download-report/{job_id}",
                    "download_details": {},
                }

                # Save to database for future requests
                save_job_to_db(job_id, response_data)

                logger.info(f"Returning completed job response for {job_id}")
                return response_data

            # Return a proper JSON response with not_found status instead of 404
            return {
                "job_id": job_id,
                "status": "not_found",
                "progress": 0,
                "input_type": None,
                "accession": None,
                "started_at": None,
                "finished_at": None,
                "error": f"Job {job_id} not found. This job may have expired, been deleted, or never existed.",
                "files_processed": None,
                "total_files": None,
                "files_downloaded": None,
                "processing_stage": "Job not found",
                "console_output": [],
                "console_errors": [],
                "html_report_urls": None,
                "download_url": None,
                "download_details": {},
            }

    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Error in job_status_api for job {job_id}: {e}")
        raise HTTPException(status_code=500, detail=f"Internal server error: {str(e)}")


@app.post("/generate-report")
async def generate_report(file: UploadFile = File(...)):
    """
    Generate a pmultiqc report from uploaded data.

    Expected form data:
    - file: ZIP file containing the data to analyze

    Returns:
    - JSON response with job ID and status
    """
    try:
        # Generate unique job ID
        job_id = str(uuid.uuid4())

        # Create job directories
        upload_dir = os.path.join(UPLOAD_FOLDER, job_id)
        output_dir = os.path.join(OUTPUT_FOLDER, job_id)
        os.makedirs(upload_dir, exist_ok=True)
        os.makedirs(output_dir, exist_ok=True)

        # Validate and save uploaded file
        zip_path, file_size = await validate_and_save_upload(file, upload_dir, job_id)

        # Extract the zip file
        extract_path = os.path.join(upload_dir, "extracted")
        os.makedirs(extract_path, exist_ok=True)

        validate_and_extract_zip(zip_path, extract_path, file_size)

        # Detect input type
        input_type, quantms_config = detect_input_type(extract_path)
        logger.info(f"Detected input type: {input_type}")

        if input_type == "unknown":
            raise HTTPException(
                status_code=400,
                detail="Could not detect input type. Please ensure your ZIP file contains valid data files.",
            )

        # Run pmultiqc
        result = run_pmultiqc_with_progress(
            extract_path, output_dir, input_type, quantms_config, job_id
        )

        if not result["success"]:
            raise HTTPException(
                status_code=500, detail=f"Failed to generate report: {result['message']}"
            )

        # Create zip report
        zip_report_path = os.path.join(output_dir, f"pmultiqc_report_{job_id}.zip")
        if not create_zip_report(output_dir, zip_report_path):
            logger.error(f"Failed to create zip report for job {job_id}")
            raise HTTPException(
                status_code=500, detail="Failed to create zip report. Please try again."
            )

        # Copy HTML report to public directory for online viewing
        try:
            html_urls = copy_html_report_for_online_viewing(output_dir, job_id)
            logger.info(f"Generated html_urls for sync job {job_id}: {html_urls}")
        except Exception as e:
            logger.error(f"Failed to copy HTML reports for sync job {job_id}: {e}")
            html_urls = None

        job_update = {
            "job_id": job_id,
            "status": "completed",
            "input_type": input_type,
            "finished_at_seconds": int(time.time()),
            "finished_at": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        }

        if html_urls:
            job_update["html_report_urls"] = html_urls

        save_job_to_db(job_id, job_update)

        base_url = BASE_URL.rstrip("/")
        response_data = {
            "job_id": job_id,
            "status": "completed",
            "input_type": input_type,
            "download_url": f"{base_url}/download-report/{job_id}",
            "message": "Report generated successfully",
            "finished_at": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        }

        if html_urls:
            response_data["html_report_urls"] = html_urls

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
        if not data or "accession" not in data:
            raise HTTPException(status_code=400, detail="Missing accession parameter")

        accession = data["accession"].strip().upper()

        # Security: Validate accession format (PXD followed by 6 digits)
        if not re.match(r'^PXD\d{6}$', accession):
            raise HTTPException(
                status_code=400, detail="Invalid PRIDE accession format. Should be PXD followed by 6 digits (e.g., PXD012345)"
            )

        # Generate unique job ID
        job_id = str(uuid.uuid4())

        # Create output directory
        output_dir = os.path.join(OUTPUT_FOLDER, job_id)
        os.makedirs(output_dir, exist_ok=True)

        # Initialize job in database
        initial_job_data = {
            "job_id": job_id,
            "status": "queued",
            "accession": accession,
            "progress": 0,
            "started_at": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "finished_at": None,
            "error": None,
            "input_type": None,
            "files_processed": 0,
            "total_files": 0,
            "files_downloaded": 0,
            "console_output": [],
            "console_errors": [],
            "results": [],
            "html_report_urls": None,
        }
        save_job_to_db(job_id, initial_job_data)

        # Start async processing
        thread = threading.Thread(
            target=process_pride_job_async, args=(job_id, accession, output_dir)
        )
        thread.daemon = True
        thread.start()

        logger.info(f"Started PRIDE processing for accession {accession} with job ID {job_id}")

        return {
            "job_id": job_id,
            "accession": accession,
            "status": "queued",
            "message": "PRIDE processing started successfully",
            "status_url": f"/job-status/{job_id}",
            "redirect_url": f"/results?job={job_id}",
        }

    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Error starting PRIDE processing: {e}")
        raise HTTPException(status_code=500, detail=f"Failed to start PRIDE processing: {str(e)}")


@app.post("/force-complete-job/{job_id}")
async def force_complete_job(job_id: str):
    """
    Force complete a job that appears to be stuck in processing.
    This endpoint checks if MultiQC has completed and manually updates the job status.
    """
    try:
        # Validate job_id format
        try:
            uuid.UUID(job_id)
        except ValueError:
            raise HTTPException(status_code=400, detail="Invalid job ID")

        logger.info(f"Force completion request for job {job_id}")

        # Check if job exists and is in a running state
        job_data = get_job_from_db(job_id)
        if not job_data:
            # Try to get from Redis
            redis_client = get_redis_client()
            redis_key = f"{REDIS_KEY_PREFIX}{job_id}"
            redis_data_json = redis_client.get(redis_key)
            if redis_data_json:
                job_data = json.loads(redis_data_json)

        if not job_data:
            raise HTTPException(status_code=404, detail=f"Job {job_id} not found")

        current_status = job_data.get("status", "unknown")
        if current_status in ["completed", "failed"]:
            return {
                "job_id": job_id,
                "message": f"Job is already {current_status}",
                "current_status": current_status,
            }

        # Check if output files exist (indicating completion)
        output_dir = os.path.join(OUTPUT_FOLDER, job_id)
        zip_report_path = os.path.join(output_dir, f"pmultiqc_report_{job_id}.zip")

        # More comprehensive completion detection
        def detect_multiqc_completion():
            """Detect if MultiQC has completed by checking for actual output files."""
            completion_indicators = {
                "console_multiqc_complete": False,
                "zip_report_exists": False,
                "html_reports_found": [],
                "html_reports_in_public": [],
                "total_html_reports": 0,
            }

            # Check console output for MultiQC complete message
            console_errors = job_data.get("console_errors", [])
            completion_indicators["console_multiqc_complete"] = any(
                "MultiQC complete" in error for error in console_errors
            )

            # Check if main zip report exists
            completion_indicators["zip_report_exists"] = os.path.exists(zip_report_path)

            # Look for MultiQC HTML reports in output directory
            if os.path.exists(output_dir):
                for root, _, files in os.walk(output_dir):
                    for file in files:
                        if file == "multiqc_report.html":
                            html_report_path = os.path.join(root, file)
                            # Get relative path from output_dir for identification
                            rel_path = os.path.relpath(root, output_dir)
                            completion_indicators["html_reports_found"].append(
                                {
                                    "path": html_report_path,
                                    "relative_dir": rel_path,
                                    "exists": os.path.exists(html_report_path),
                                    "size": (
                                        os.path.getsize(html_report_path)
                                        if os.path.exists(html_report_path)
                                        else 0
                                    ),
                                }
                            )

            completion_indicators["total_html_reports"] = len(
                completion_indicators["html_reports_found"]
            )

            # Check if HTML reports have been copied to public directory
            # Look for existing HTML report hashes
            if os.path.exists(HTML_REPORTS_FOLDER):
                for item in os.listdir(HTML_REPORTS_FOLDER):
                    hash_dir = os.path.join(HTML_REPORTS_FOLDER, item)
                    if os.path.isdir(hash_dir):
                        job_file = os.path.join(hash_dir, ".job_id")
                        if os.path.exists(job_file):
                            try:
                                with open(job_file, "r") as f:
                                    stored_job_id = f.read().strip()
                                if stored_job_id == job_id:
                                    # Check for HTML reports in this directory
                                    for root, _, files in os.walk(hash_dir):
                                        for file in files:
                                            if file == "multiqc_report.html":
                                                completion_indicators[
                                                    "html_reports_in_public"
                                                ].append(
                                                    {
                                                        "hash": item,
                                                        "path": os.path.join(root, file),
                                                        "relative_dir": os.path.relpath(
                                                            root, hash_dir
                                                        ),
                                                    }
                                                )
                            except Exception as e:
                                logger.warning(f"Could not read job file {job_file}: {e}")

            return completion_indicators

        completion_info = detect_multiqc_completion()
        logger.info(f"Completion detection for job {job_id}: {completion_info}")

        # Determine if job should be considered completed
        job_appears_completed = (
            completion_info["console_multiqc_complete"]
            or completion_info["zip_report_exists"]
            or completion_info["total_html_reports"] > 0
        )

        if job_appears_completed:
            logger.info(f"Force completing job {job_id} based on file detection")
            logger.info(
                f"  - Console MultiQC complete: {completion_info['console_multiqc_complete']}"
            )
            logger.info(f"  - Zip report exists: {completion_info['zip_report_exists']}")
            logger.info(f"  - HTML reports found: {completion_info['total_html_reports']}")
            logger.info(
                f"  - HTML reports in public: {len(completion_info['html_reports_in_public'])}"
            )

            # Force update to completed
            update_job_progress(
                job_id,
                "completed",
                100,
                processing_stage="Analysis completed successfully! (Force completed)",
                finished_at=datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            )

            # Try to copy HTML reports if they're not already in public directory
            html_urls = None
            if (
                completion_info["total_html_reports"] > 0
                and len(completion_info["html_reports_in_public"]) == 0
            ):
                logger.info(
                    f"Found {completion_info['total_html_reports']} HTML reports not yet in public directory, copying..."
                )
                try:
                    html_urls = copy_html_report_for_online_viewing(
                        output_dir, job_id, job_data.get("accession")
                    )
                    if html_urls:
                        update_job_progress(job_id, "completed", 100, html_report_urls=html_urls)
                        logger.info(
                            f"Successfully copied HTML reports to public directory: {html_urls}"
                        )
                    else:
                        logger.warning("HTML report copying returned None")
                except Exception as e:
                    logger.warning(f"Failed to copy HTML reports during force completion: {e}")
            elif len(completion_info["html_reports_in_public"]) > 0:
                logger.info("HTML reports already exist in public directory")
                # Reconstruct HTML URLs from existing public reports
                try:
                    base_url = BASE_URL.rstrip("/")
                    html_urls = {}
                    for report_info in completion_info["html_reports_in_public"]:
                        report_name = (
                            report_info["relative_dir"]
                            if report_info["relative_dir"] != "."
                            else "root"
                        )
                        html_url = f"{base_url}/reports/{report_info['hash']}/{report_info['relative_dir']}/multiqc_report.html"
                        if report_info["relative_dir"] == ".":
                            html_url = (
                                f"{base_url}/reports/{report_info['hash']}/multiqc_report.html"
                            )
                        html_urls[report_name] = html_url

                    if html_urls:
                        update_job_progress(job_id, "completed", 100, html_report_urls=html_urls)
                        logger.info(
                            f"Reconstructed HTML URLs from existing public reports: {html_urls}"
                        )
                except Exception as e:
                    logger.warning(f"Failed to reconstruct HTML URLs: {e}")

            return {
                "job_id": job_id,
                "message": "Job force completed successfully",
                "previous_status": current_status,
                "new_status": "completed",
                "completion_details": completion_info,
                "html_urls_generated": html_urls is not None,
                "html_urls": html_urls or {},
            }
        else:
            raise HTTPException(
                status_code=400,
                detail="Job does not appear to be completed - no MultiQC completion detected",
            )

    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Error in force_complete_job for job {job_id}: {e}")
        raise HTTPException(status_code=500, detail=f"Internal server error: {str(e)}")


@app.get("/job-diagnostic/{job_id}")
async def job_diagnostic(job_id: str):
    """
    Comprehensive diagnostic information for a job to help debug completion issues.
    """
    try:
        # Validate job_id format
        try:
            uuid.UUID(job_id)
        except ValueError:
            raise HTTPException(status_code=400, detail="Invalid job ID")

        logger.info(f"Diagnostic request for job {job_id}")

        # Get job data from both Redis and database
        redis_data = None
        db_data = None

        try:
            redis_client = get_redis_client()
            redis_key = f"{REDIS_KEY_PREFIX}{job_id}"
            redis_data_json = redis_client.get(redis_key)
            if redis_data_json:
                redis_data = json.loads(redis_data_json)
        except Exception as e:
            logger.warning(f"Could not get Redis data for job {job_id}: {e}")

        try:
            db_data = get_job_from_db(job_id)
        except Exception as e:
            logger.warning(f"Could not get database data for job {job_id}: {e}")

        # Check file system
        output_dir = os.path.join(OUTPUT_FOLDER, job_id)
        upload_dir = os.path.join(UPLOAD_FOLDER, job_id)

        file_system_info = {
            "output_dir_exists": os.path.exists(output_dir),
            "upload_dir_exists": os.path.exists(upload_dir),
            "output_files": [],
            "html_reports": [],
            "zip_reports": [],
            "public_html_reports": [],
        }

        # Scan output directory
        if os.path.exists(output_dir):
            for root, _, files in os.walk(output_dir):
                for file in files:
                    file_path = os.path.join(root, file)
                    rel_path = os.path.relpath(file_path, output_dir)
                    file_info = {
                        "path": rel_path,
                        "full_path": file_path,
                        "size": os.path.getsize(file_path),
                        "modified": os.path.getmtime(file_path),
                    }

                    file_system_info["output_files"].append(file_info)

                    if file == "multiqc_report.html":
                        file_system_info["html_reports"].append(file_info)
                    elif file.endswith(".zip"):
                        file_system_info["zip_reports"].append(file_info)

        # Check public HTML reports
        if os.path.exists(HTML_REPORTS_FOLDER):
            for item in os.listdir(HTML_REPORTS_FOLDER):
                hash_dir = os.path.join(HTML_REPORTS_FOLDER, item)
                if os.path.isdir(hash_dir):
                    job_file = os.path.join(hash_dir, ".job_id")
                    if os.path.exists(job_file):
                        try:
                            with open(job_file, "r") as f:
                                stored_job_id = f.read().strip()
                            if stored_job_id == job_id:
                                # Found public reports for this job
                                for root, _, files in os.walk(hash_dir):
                                    for file in files:
                                        if file == "multiqc_report.html":
                                            rel_path = os.path.relpath(root, hash_dir)
                                            base_url = BASE_URL.rstrip("/")
                                            if rel_path == ".":
                                                url = f"{base_url}/reports/{item}/multiqc_report.html"
                                                report_name = "root"
                                            else:
                                                url = f"{base_url}/reports/{item}/{rel_path}/multiqc_report.html"
                                                report_name = rel_path

                                            file_system_info["public_html_reports"].append(
                                                {
                                                    "hash": item,
                                                    "report_name": report_name,
                                                    "relative_path": rel_path,
                                                    "url": url,
                                                    "file_size": os.path.getsize(
                                                        os.path.join(root, file)
                                                    ),
                                                }
                                            )
                        except Exception as e:
                            logger.warning(f"Error checking public report {hash_dir}: {e}")

        # Determine completion status
        appears_completed = (
            len(file_system_info["html_reports"]) > 0
            or len(file_system_info["zip_reports"]) > 0
            or (redis_data and "MultiQC complete" in str(redis_data.get("console_errors", [])))
            or (db_data and "MultiQC complete" in str(db_data.get("console_errors", [])))
        )

        return {
            "job_id": job_id,
            "redis_data": redis_data,
            "database_data": db_data,
            "file_system": file_system_info,
            "appears_completed": appears_completed,
            "completion_indicators": {
                "html_reports_generated": len(file_system_info["html_reports"]) > 0,
                "zip_reports_generated": len(file_system_info["zip_reports"]) > 0,
                "public_html_available": len(file_system_info["public_html_reports"]) > 0,
                "total_output_files": len(file_system_info["output_files"]),
            },
            "recommended_action": (
                "force_complete"
                if appears_completed
                and (
                    (redis_data and redis_data.get("status") not in ["completed", "failed"])
                    or (db_data and db_data.get("status") not in ["completed", "failed"])
                    or (not redis_data and not db_data)
                )
                else "none"
            ),
        }

    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Error in job_diagnostic for job {job_id}: {e}")
        raise HTTPException(status_code=500, detail=f"Internal server error: {str(e)}")


@app.get("/download-report/{job_id}")
async def download_report_api(job_id: str):
    """API version of download-report endpoint."""
    return await download_report(job_id)


async def download_report(job_id: str):
    """
    Download the generated pmultiqc report.

    Args:
        job_id: The job ID returned from generate-report

    Returns:
        ZIP file containing the pmultiqc report
    """
    try:
        # Validate job_id format
        try:
            uuid.UUID(job_id)
        except ValueError:
            raise HTTPException(status_code=400, detail="Invalid job ID")

        output_dir = os.path.join(OUTPUT_FOLDER, job_id)
        zip_report_path = os.path.join(output_dir, f"pmultiqc_report_{job_id}.zip")

        if not os.path.exists(zip_report_path):
            raise HTTPException(status_code=404, detail="Report not found")

        return FileResponse(
            path=zip_report_path,
            filename=f"pmultiqc_report_{job_id}.zip",
            media_type="application/zip",
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
    return get_swagger_ui_html(
        openapi_url=f"{BASE_URL.rstrip('/')}/openapi.json",
        title=f"{app.title} - Swagger UI",
        swagger_js_url="https://cdn.jsdelivr.net/npm/swagger-ui-dist@5.9.0/swagger-ui-bundle.js",
        swagger_css_url="https://cdn.jsdelivr.net/npm/swagger-ui-dist@5.9.0/swagger-ui.css",
    )


if __name__ == "__main__":
    uvicorn.run(app, host="0.0.0.0", port=5000)