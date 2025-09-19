# Quick Deployment Guide

This guide shows how to deploy PMultiQC with the correct BASE_URL configuration for your domain, specifically for deployments like https://abi-services.cs.uni-tuebingen.de/pmultiqc/.

## Before (Problem)

The old docker-compose.yml had hardcoded BASE_URL:
```yaml
environment:
  - BASE_URL=http://localhost:8080  # Wrong for production!
```

This caused issues:
- ❌ Static assets (CSS/JS) wouldn't load from correct URLs
- ❌ API calls from frontend would fail  
- ❌ Download links would be broken
- ❌ OpenAPI docs would show wrong server URL

## After (Solution)

### Option 1: Use environment variable override
```bash
BASE_URL=https://abi-services.cs.uni-tuebingen.de/pmultiqc docker compose up -d
```

### Option 2: Use production configuration with .env file
```bash
# Copy the configuration template
cp .env.abi-services.example .env

# Deploy with production settings
docker compose -f docker-compose.prod.yml up -d
```

## What this fixes

✅ **Static assets load correctly**: CSS and JavaScript files load from the right BASE_URL  
✅ **API calls work**: Frontend can make API calls to the correct backend URL  
✅ **Download links work**: Report downloads use the correct public URL  
✅ **OpenAPI docs work**: Swagger UI shows the correct server URL  
✅ **Redirects work**: Internal redirects use the correct base URL  

## Example configuration for your deployment

```env
# .env file
BASE_URL=https://abi-services.cs.uni-tuebingen.de/pmultiqc
HOST_PORT=8080
FLASK_ENV=production
LOG_LEVEL=INFO
```

## Equivalent to Kubernetes ConfigMap

This Docker Compose approach provides the same functionality as the Kubernetes ConfigMap:

**Kubernetes ConfigMap:**
```yaml
apiVersion: v1
kind: ConfigMap
metadata:
  name: pmultiqc-config
data:
  BASE_URL: "https://www.ebi.ac.uk/pride/services/pmultiqc"
```

**Docker Compose .env:**
```env
BASE_URL=https://abi-services.cs.uni-tuebingen.de/pmultiqc
```

Both approaches:
- ✅ Allow configuring BASE_URL for different deployments
- ✅ Support all the same environment variables
- ✅ Enable proper subpath deployments
- ✅ Work with reverse proxies