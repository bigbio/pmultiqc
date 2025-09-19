# PMultiQC Docker Compose Deployment Guide

This directory contains Docker Compose configurations for deploying PMultiQC service in different environments.

## Quick Start

### For Local Development

1. Copy the development example environment file:
   ```bash
   cp .env.development.example .env
   ```

2. Start the services:
   ```bash
   docker-compose up -d
   ```

3. Access the service at: http://localhost:8080

### For Production Deployment

1. Copy and customize the production example:
   ```bash
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

## Configuration Files

### Environment Files

- **`.env.template`** - Template with all available configuration options
- **`.env.development.example`** - Pre-configured for local development  
- **`.env.production.example`** - Pre-configured for production deployment
- **`.env`** - Your actual configuration (not in version control)

### Docker Compose Files

- **`docker-compose.yml`** - Basic configuration with environment variable support
- **`docker-compose.prod.yml`** - Production-optimized configuration with .env file support
- **`docker-compose.dev.yml`** - Development configuration (legacy)

## Critical Configuration: BASE_URL

The `BASE_URL` environment variable is **critical** for proper operation. It must be set to the exact URL where users will access your PMultiQC service.

### Examples:

- **Local testing**: `BASE_URL=http://localhost:8080`
- **Domain root**: `BASE_URL=https://pmultiqc.example.com`
- **Subpath deployment**: `BASE_URL=https://example.com/services/pmultiqc`
- **Reverse proxy**: `BASE_URL=https://abi-services.cs.uni-tuebingen.de/pmultiqc`

### Why BASE_URL matters:

- **Static assets**: CSS, JavaScript files must load from correct paths
- **API endpoints**: Frontend JavaScript makes API calls to correct URLs
- **Downloads**: Report download links must point to accessible URLs
- **Redirects**: Internal redirects must use the correct base URL
- **OpenAPI documentation**: Swagger UI must know the correct server URL

## Deployment Scenarios

### 1. Simple Domain Deployment

If deploying at a domain root (e.g., `https://pmultiqc.example.com`):

```env
BASE_URL=https://pmultiqc.example.com
HOST_PORT=80
```

### 2. Subpath Deployment

If deploying at a subpath (e.g., `https://example.com/pmultiqc`):

```env
BASE_URL=https://example.com/pmultiqc
HOST_PORT=8080
```

**Note**: Your reverse proxy must strip the subpath when forwarding to the container.

### 3. Behind Reverse Proxy

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

## Available Environment Variables

| Variable | Default | Description |
|----------|---------|-------------|
| `BASE_URL` | `http://localhost:8080` | **CRITICAL**: Public URL for your deployment |
| `HOST_PORT` | `8080` | Host port to bind container port 5000 |
| `FLASK_ENV` | `production` | Flask environment (development/production) |
| `MAX_CONTENT_LENGTH` | `4294967296` | Max upload size in bytes (4GB) |
| `LOG_LEVEL` | `INFO` | Logging level (DEBUG/INFO/WARNING/ERROR) |
| `REDIS_URL` | `redis://redis:6379` | Redis connection URL |
| `MEMORY_LIMIT` | `16384M` | Container memory limit |
| `CPU_LIMIT` | `4.0` | Container CPU limit |

## Troubleshooting

### Static Assets Not Loading

**Problem**: CSS/JavaScript files return 404 errors
**Solution**: Check that `BASE_URL` matches exactly how users access your service

### API Calls Failing

**Problem**: Frontend can't connect to backend API
**Solution**: Verify `BASE_URL` is accessible from user browsers

### Download Links Broken

**Problem**: Report download links don't work
**Solution**: Ensure `BASE_URL` points to the correct public URL

### Service Not Accessible

**Problem**: Can't access the service after deployment
**Solution**: 
1. Check `HOST_PORT` is correctly mapped
2. Verify firewall/security group settings
3. Check reverse proxy configuration

## Health Checks

The service includes health checks at `/health`. You can verify deployment:

```bash
curl http://localhost:8080/health
```

For production deployments, replace `localhost:8080` with your actual `BASE_URL`.

## Logs

View service logs:
```bash
docker-compose logs -f pmultiqc-service
```

## Comparison with Kubernetes

This Docker Compose approach provides similar functionality to the Kubernetes ConfigMap approach:

- **Kubernetes**: Uses ConfigMap to set environment variables
- **Docker Compose**: Uses .env files to set environment variables
- **Both**: Support the same configuration options and BASE_URL customization

The key difference is deployment complexity - Docker Compose is simpler for single-host deployments, while Kubernetes provides better scaling and orchestration for multi-host deployments.