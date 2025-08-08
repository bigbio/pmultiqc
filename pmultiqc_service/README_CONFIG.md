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
  pmultiqc:
    image: pmultiqc-service
    environment:
      - MAX_CONTENT_LENGTH=4294967296
      - UPLOAD_FOLDER=/app/uploads
      - OUTPUT_FOLDER=/app/outputs
      - HTML_REPORTS_FOLDER=/app/reports
    volumes:
      - ./config.yaml:/app/config.yaml
      - ./uploads:/app/uploads
      - ./outputs:/app/outputs
      - ./reports:/app/reports
```

### Kubernetes

```yaml
apiVersion: v1
kind: ConfigMap
metadata:
  name: pmultiqc-config
data:
  config.yaml: |
    max_content_length: 4294967296
    upload_folder: /data/uploads
    output_folder: /data/outputs
    html_reports_folder: /data/reports
    base_url: https://www.ebi.ac.uk/pride/services/pmultiqc
---
apiVersion: apps/v1
kind: Deployment
metadata:
  name: pmultiqc
spec:
  template:
    spec:
      containers:
      - name: pmultiqc
        image: pmultiqc-service
        env:
        - name: CONFIG_FILE
          value: /etc/pmultiqc/config.yaml
        volumeMounts:
        - name: config
          mountPath: /etc/pmultiqc
      volumes:
      - name: config
        configMap:
          name: pmultiqc-config
```

### Docker Run

```bash
docker run -d \
  -e MAX_CONTENT_LENGTH=4294967296 \
  -e UPLOAD_FOLDER=/data/uploads \
  -e OUTPUT_FOLDER=/data/outputs \
  -e HTML_REPORTS_FOLDER=/data/reports \
  -e BASE_URL=https://www.ebi.ac.uk/pride/services/pmultiqc \
  -v /path/to/config.yaml:/app/config.yaml \
  -v /path/to/data:/data \
  pmultiqc-service
```

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