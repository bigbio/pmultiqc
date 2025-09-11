FROM python:3.10-slim

# Set working directory
WORKDIR /app

# Install system dependencies
RUN apt-get update && apt-get install -y \
    gcc \
    g++ \
    git \
    libxml2-dev \
    libxslt-dev \
    zlib1g-dev \
    libhdf5-dev \
    libnetcdf-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libffi-dev \
    libglib2.0-0 \
    && rm -rf /var/lib/apt/lists/*

# Upgrade pip and install wheel
RUN pip install --upgrade pip setuptools wheel

# Copy requirements and install Python dependencies
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Install pmultiqc from GitHub main branch (overrides PyPI version)
RUN pip install --no-cache-dir --force-reinstall git+https://github.com/bigbio/pmultiqc.git@main

# Copy service files
COPY app.py .
COPY templates/ ./templates/

# Create necessary directories (will be created by the application using environment variables)
RUN mkdir -p /app/pmultiqc_uploads /app/pmultiqc_outputs /app/pmultiqc_html_reports

# Expose port
EXPOSE 5000

# Set environment variables
ENV FLASK_APP=app.py
ENV UPLOAD_FOLDER=/app/pmultiqc_uploads
ENV OUTPUT_FOLDER=/app/pmultiqc_outputs
ENV HTML_REPORTS_FOLDER=/app/pmultiqc_html_reports

# Run the application
CMD ["uvicorn", "app:app", "--host", "0.0.0.0", "--port", "5000", "--workers", "2"]