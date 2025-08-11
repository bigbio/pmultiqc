# PMultiQC Service Configuration

The PMultiQC service supports multiple configuration methods with the following priority order:

1. **Environment Variables** (highest priority)
2. **YAML Configuration File**
3. **Default Values** (lowest priority)

## Configuration Methods

### 1. Environment Variables

Set environment variables to override configuration:

```bash
export MAX_CONTENT_LENGTH=4294967296  # 4GB
export UPLOAD_FOLDER=/data/uploads
export OUTPUT_FOLDER=/data/outputs
export HTML_REPORTS_FOLDER=/data/reports
export CONFIG_FILE=/path/to/config.yaml
```

### 2. YAML Configuration File

Create a `config.yaml` file in one of these locations (in order of preference):

- Path specified by `CONFIG_FILE` environment variable
- `config.yaml` in current directory
- `config.yml` in current directory
- `/etc/pmultiqc/config.yaml` (system-wide)
- `config.yaml` next to `app.py`

Example `config.yaml`:

```yaml
# Maximum file upload size in bytes (default: 2GB)
max_content_length: 2147483648

# Directory paths for file storage
upload_folder: /tmp/pmultiqc_uploads
output_folder: /tmp/pmultiqc_outputs
html_reports_folder: /tmp/pmultiqc_html_reports

# Optional: Additional configuration options
logging:
  level: INFO

cleanup:
  expiration_seconds: 86400  # 24 hours
  cleanup_job_seconds: 3600  # 1 hour
```

### 3. Default Values

If no environment variables or YAML file is found, the service uses these defaults:

- `MAX_CONTENT_LENGTH`: 2GB (2147483648 bytes)
- `UPLOAD_FOLDER`: `/tmp/pmultiqc_uploads`
- `OUTPUT_FOLDER`: `/tmp/pmultiqc_outputs`
- `HTML_REPORTS_FOLDER`: `/tmp/pmultiqc_html_reports`

## Deployment Examples

### Docker Compose

```yaml
version: '3.8'
services:
  pmultiqc-service:
    build:
      context: ../
      dockerfile: Dockerfile
    image: pmultiqc-service
    container_name: pmultiqc-service
    ports:
      - "5000:5000"
    environment:
      - FLASK_ENV=production
      - MAX_CONTENT_LENGTH=4294967296  # 4GB upload limit
      - UPLOAD_FOLDER=/app/uploads
      - OUTPUT_FOLDER=/app/outputs
      - HTML_REPORTS_FOLDER=/app/reports
      - BASE_URL=http://localhost:5000
      - LOG_LEVEL=INFO
    volumes:
      - ./uploads:/app/uploads
      - ./outputs:/app/outputs
      - ./reports:/app/reports
      - ./config.yaml:/app/config.yaml:ro
    restart: unless-stopped
    deploy:
      resources:
        limits:
          memory: 8Gi
          cpus: '4.0'
        reservations:
          memory: 4Gi
          cpus: '1.0'
    healthcheck:
      test: ["CMD", "curl", "-f", "http://localhost:5000/health"]
      interval: 30s
      timeout: 10s
      retries: 3
      start_period: 40s
```

### Kubernetes

```yaml
apiVersion: v1
kind: ConfigMap
metadata:
  name: pmultiqc-config
  namespace: pmultiqc
data:
  FLASK_ENV: "production"
  MAX_CONTENT_LENGTH: "4294967296"  # 4GB in bytes
  UPLOAD_FOLDER: "/tmp/pmultiqc_uploads"
  OUTPUT_FOLDER: "/tmp/pmultiqc_outputs"
  HTML_REPORTS_FOLDER: "/tmp/pmultiqc_html_reports"
  BASE_URL: "https://www.ebi.ac.uk/pride/services/pmultiqc"
  LOG_LEVEL: "INFO"
---
apiVersion: apps/v1
kind: Deployment
metadata:
  name: pmultiqc-service
  namespace: pmultiqc
spec:
  replicas: 2
  selector:
    matchLabels:
      app: pmultiqc-service
  template:
    metadata:
      labels:
        app: pmultiqc-service
    spec:
      containers:
      - name: pmultiqc-service
        image: pmultiqc-service:latest
        ports:
        - containerPort: 5000
        env:
        - name: FLASK_ENV
          valueFrom:
            configMapKeyRef:
              name: pmultiqc-config
              key: FLASK_ENV
        - name: MAX_CONTENT_LENGTH
          valueFrom:
            configMapKeyRef:
              name: pmultiqc-config
              key: MAX_CONTENT_LENGTH
        - name: UPLOAD_FOLDER
          valueFrom:
            configMapKeyRef:
              name: pmultiqc-config
              key: UPLOAD_FOLDER
        - name: OUTPUT_FOLDER
          valueFrom:
            configMapKeyRef:
              name: pmultiqc-config
              key: OUTPUT_FOLDER
        - name: HTML_REPORTS_FOLDER
          valueFrom:
            configMapKeyRef:
              name: pmultiqc-config
              key: HTML_REPORTS_FOLDER
        - name: BASE_URL
          valueFrom:
            configMapKeyRef:
              name: pmultiqc-config
              key: BASE_URL
        resources:
          limits:
            cpu: '4'
            memory: 8Gi
          requests:
            cpu: 1000m
            memory: 8Gi
        livenessProbe:
          httpGet:
            path: /health
            port: 5000
          initialDelaySeconds: 120
          periodSeconds: 60
        readinessProbe:
          httpGet:
            path: /health
            port: 5000
          initialDelaySeconds: 30
          periodSeconds: 15
```

### Docker Run

```bash
docker run -d \
  --memory=8g \
  --cpus=4.0 \
  -e FLASK_ENV=production \
  -e MAX_CONTENT_LENGTH=4294967296 \
  -e UPLOAD_FOLDER=/data/uploads \
  -e OUTPUT_FOLDER=/data/outputs \
  -e HTML_REPORTS_FOLDER=/data/reports \
  -e BASE_URL=https://www.ebi.ac.uk/pride/services/pmultiqc \
  -e LOG_LEVEL=INFO \
  -v /path/to/config.yaml:/app/config.yaml:ro \
  -v /path/to/data:/data \
  -p 5000:5000 \
  pmultiqc-service
```

## System Requirements

### Memory Requirements
- **Minimum Memory**: 8GB RAM
- **Recommended Memory**: 16GB RAM for production workloads
- **Memory Usage**: The application processes large multiQC files and requires significant memory for optimal performance

### CPU Requirements
- **Minimum CPU**: 1 CPU core
- **Recommended CPU**: 4 CPU cores for production workloads

## Configuration Options

| Option | Environment Variable | YAML Key | Default | Description |
|--------|---------------------|----------|---------|-------------|
| Max Upload Size | `MAX_CONTENT_LENGTH` | `max_content_length` | 2GB | Maximum file upload size in bytes |
| Upload Directory | `UPLOAD_FOLDER` | `upload_folder` | `/tmp/pmultiqc_uploads` | Directory for uploaded files |
| Output Directory | `OUTPUT_FOLDER` | `output_folder` | `/tmp/pmultiqc_outputs` | Directory for generated reports |
| HTML Reports Directory | `HTML_REPORTS_FOLDER` | `html_reports_folder` | `/tmp/pmultiqc_html_reports` | Directory for HTML reports |
| Base URL | `BASE_URL` | `base_url` | `http://localhost:5000` | Base URL for generating absolute URLs |

## Dependencies

The YAML configuration feature requires the `PyYAML` package. If not installed, the service will:

1. Log a warning that YAML is not available
2. Continue using environment variables and defaults
3. Not fail to start

To install PyYAML:

```bash
pip install PyYAML
```

## Troubleshooting

### Check Configuration

The service logs all configuration values at startup. Look for lines like:

```
Configuration loaded:
  MAX_CONTENT_LENGTH: 2147483648 bytes (2.0 GB)
  UPLOAD_FOLDER: /tmp/pmultiqc_uploads
  OUTPUT_FOLDER: /tmp/pmultiqc_outputs
  HTML_REPORTS_FOLDER: /tmp/pmultiqc_html_reports
  BASE_URL: http://localhost:5000
```

### Debug Configuration Loading

Set log level to DEBUG to see detailed configuration loading:

```bash
export LOG_LEVEL=DEBUG
```

Or in YAML:

```yaml
logging:
  level: DEBUG
``` 