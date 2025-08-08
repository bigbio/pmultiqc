FROM python:3.10-slim

# Set working directory
WORKDIR /app

# Install system dependencies
RUN apt-get update && apt-get install -y \
    gcc \
    g++ \
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

# Copy service files
COPY app.py .
COPY templates/ ./templates/

# Create necessary directories
RUN mkdir -p /tmp/pmultiqc_uploads /tmp/pmultiqc_outputs

# Expose port
EXPOSE 5000

# Set environment variables
ENV FLASK_APP=app.py
ENV FLASK_ENV=production

# Run the application
CMD ["uvicorn", "app:app", "--host", "0.0.0.0", "--port", "5000", "--workers", "2"]