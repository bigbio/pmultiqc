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

This section provides comprehensive Docker Compose configurations for deploying PMultiQC service in different environments.

#### Quick Start

**For Local Development:**

1. Copy the development example environment file:
   ```bash
   cd pmultiqc_service/docker_compose
   cp .env.development.example .env
   ```

2. Start the services:
   ```bash
   docker-compose up -d
   ```

3. Access the service at: http://localhost:8080

**For Production Deployment:**

1. Copy and customize the production example:
   ```bash
   cd pmultiqc_service/docker_compose
   cp .env.production.example .env
   ```

2. **IMPORTANT**: Edit `.env` and set the `BASE_URL` to your actual deployment URL:
   ```bash
   # For example, if deploying at https://your-domain.com/pmultiqc
   BASE_URL=https://your-domain.com/pmultiqc
   ```

3. Start the services using the production configuration:
   ```bash
   docker-compose -f docker-compose.prod.yml up -d
   ```

#### Configuration Files

**Environment Files:**
- **`.env.template`** - Template with all available configuration options
- **`.env.development.example`** - Pre-configured for local development  
- **`.env.production.example`** - Pre-configured for production deployment
- **`.env`** - Your actual configuration (not in version control)

**Docker Compose Files:**
- **`docker-compose.yml`** - Basic configuration with environment variable support
- **`docker-compose.prod.yml`** - Production-optimized configuration with .env file support
- **`docker-compose.dev.yml`** - Development configuration (legacy)

#### Critical Configuration: BASE_URL

The `BASE_URL` environment variable is **critical** for proper operation. It must be set to the exact URL where users will access your PMultiQC service.

**Examples:**
- **Local testing**: `BASE_URL=http://localhost:8080`
- **Domain root**: `BASE_URL=https://pmultiqc.example.com`
- **Subpath deployment**: `BASE_URL=https://example.com/services/pmultiqc`
- **Reverse proxy**: `BASE_URL=https://abi-services.cs.uni-tuebingen.de/pmultiqc`

**Why BASE_URL matters:**
- **Static assets**: CSS, JavaScript files must load from correct paths
- **API endpoints**: Frontend JavaScript makes API calls to correct URLs
- **Downloads**: Report download links must point to accessible URLs
- **Redirects**: Internal redirects must use the correct base URL
- **OpenAPI documentation**: Swagger UI must know the correct server URL

#### Deployment Scenarios

**1. Simple Domain Deployment**

If deploying at a domain root (e.g., `https://pmultiqc.example.com`):

```env
BASE_URL=https://pmultiqc.example.com
HOST_PORT=80
```

**2. Subpath Deployment**

If deploying at a subpath (e.g., `https://example.com/pmultiqc`):

```env
BASE_URL=https://example.com/pmultiqc
HOST_PORT=8080
```

**Note**: Your reverse proxy must strip the subpath when forwarding to the container.

**3. Behind Reverse Proxy**

If using nginx, Apache, or another reverse proxy:

```env
BASE_URL=https://your-domain.com/pmultiqc
HOST_PORT=8080
```

Example nginx configuration:
```nginx
location /pmultiqc/ {
    proxy_pass http://localhost:8080/;
    proxy_set_header Host $host;
    proxy_set_header X-Real-IP $remote_addr;
    proxy_set_header X-Forwarded-For $proxy_add_x_forwarded_for;
    proxy_set_header X-Forwarded-Proto $scheme;
}
```

#### Available Environment Variables

| Variable | Default | Description |
|----------|---------|-------------|
| `BASE_URL` | `http://localhost:8080` | **CRITICAL**: Public URL for your deployment |
| `HOST_PORT` | `8080` | Host port to bind container port 5000 |
| `FLASK_ENV` | `production` | Flask environment (development/production) |
| `MAX_CONTENT_LENGTH` | `4294967296` | Max upload size in bytes (4GB) |
| `LOG_LEVEL` | `INFO` | Logging level (DEBUG/INFO/WARNING/ERROR) |
| `REDIS_URL` | `redis://redis:6379` | Redis connection URL |
| `REDIS_HOST_PORT` | `6379` | Host port to bind Redis container port 6379 |
| `MEMORY_LIMIT` | `16384M` | Container memory limit |
| `CPU_LIMIT` | `4.0` | Container CPU limit |

#### Docker Compose Troubleshooting

**Static Assets Not Loading**
- **Problem**: CSS/JavaScript files return 404 errors
- **Solution**: Check that `BASE_URL` matches exactly how users access your service

**API Calls Failing**
- **Problem**: Frontend can't connect to backend API
- **Solution**: Verify `BASE_URL` is accessible from user browsers

**Download Links Broken**
- **Problem**: Report download links don't work
- **Solution**: Ensure `BASE_URL` points to the correct public URL

**Service Not Accessible**
- **Problem**: Can't access the service after deployment
- **Solution**: 
  1. Check `HOST_PORT` is correctly mapped
  2. Verify firewall/security group settings
  3. Check reverse proxy configuration

#### Health Checks and Logs

The service includes health checks at `/health`. You can verify deployment:

```bash
curl http://localhost:8080/health
```

For production deployments, replace `localhost:8080` with your actual `BASE_URL`.

View service logs:
```bash
docker-compose logs -f pmultiqc-service
```

#### Comparison with Kubernetes

This Docker Compose approach provides similar functionality to the Kubernetes ConfigMap approach:

- **Kubernetes**: Uses ConfigMap to set environment variables
- **Docker Compose**: Uses .env files to set environment variables
- **Both**: Support the same configuration options and BASE_URL customization

The key difference is deployment complexity - Docker Compose is simpler for single-host deployments, while Kubernetes provides better scaling and orchestration for multi-host deployments.

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