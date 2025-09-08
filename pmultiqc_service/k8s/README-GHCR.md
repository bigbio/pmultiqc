# GitHub Container Registry (GHCR) Setup for PMultiQC

This document explains how to use GitHub Container Registry for PMultiQC Docker images instead of Docker Hub.

## Overview

The PMultiQC service has been configured to build and push Docker images to GitHub Container Registry (ghcr.io) with the following benefits:

- **Cross-platform builds**: Automatically builds for `linux/amd64` (required for Kubernetes)
- **No local Docker required**: Builds happen in GitHub Actions runners
- **Integrated with GitHub**: Uses GitHub's built-in container registry
- **Free for public repositories**: No additional costs for public repos

## Image Location

- **Registry**: `ghcr.io`
- **Repository**: `bigbio/pmultiqc`
- **Full image name**: `ghcr.io/bigbio/pmultiqc:latest`

## GitHub Action Workflow

The workflow file `.github/workflows/build-and-push.yml` automatically:

1. **Triggers on**:
   - Push to `main` or `master` branch
   - Tag pushes (e.g., `v1.0.0`)
   - Pull requests
   - Manual workflow dispatch

2. **Builds for**:
   - `linux/amd64` platform (required for Kubernetes)

3. **Pushes to**:
   - `ghcr.io/bigbio/pmultiqc:latest`
   - `ghcr.io/bigbio/pmultiqc:simple`
   - Branch-specific tags
   - SHA-based tags

**Note**: The workflow only builds and pushes images. Manual deployment to Kubernetes is required.

## Kubernetes Deployment

### Current Configuration

The Kubernetes deployment has been updated to use the GHCR image:

```yaml
containers:
  - name: pmultiqc-service
    image: ghcr.io/bigbio/pmultiqc:latest
```

### Updating the Deployment

To update your Kubernetes cluster to use the new GHCR image:

```bash
# Navigate to the k8s directory
cd pmultiqc_service/k8s

# Run the update script
./update-to-ghcr.sh
```


## Authentication

### For GitHub Actions

The workflow uses the built-in `GITHUB_TOKEN` which has the necessary permissions to push to GHCR.

### For Kubernetes

If your Kubernetes cluster needs to pull from GHCR, you may need to create a secret:

```bash
kubectl create secret docker-registry ghcr-secret \
  --docker-server=ghcr.io \
  --docker-username=<YOUR_GITHUB_USERNAME> \
  --docker-password=<YOUR_GITHUB_TOKEN> \
  --docker-email=<YOUR_EMAIL> \
  -n pmultiqc
```

Then update the deployment to use the secret:

```yaml
spec:
  template:
    spec:
      imagePullSecrets:
        - name: ghcr-secret
```

## Benefits of GHCR vs Docker Hub

| Feature | GHCR | Docker Hub |
|---------|------|------------|
| Cross-platform builds | ✅ Native support | ❌ Requires setup |
| Free for public repos | ✅ Unlimited | ❌ Limited |
| GitHub integration | ✅ Built-in | ❌ External |
| Rate limits | ✅ Generous | ❌ Strict |
| Security scanning | ✅ Included | ❌ Paid feature |

## Troubleshooting

### Image Pull Errors

If you get image pull errors:

1. Check if the image exists: https://github.com/bigbio/pmultiqc/pkgs/container/pmultiqc
2. Verify the image name in the deployment
3. Check if authentication is required

### Build Failures

If the GitHub Action fails:

1. Check the Actions tab in your GitHub repository
2. Verify the Dockerfile is correct
3. Check for any dependency issues

### Deployment Issues

If the Kubernetes deployment fails:

1. Check pod logs: `kubectl logs -n pmultiqc deployment/pmultiqc-service`
2. Verify the image exists and is accessible
3. Check resource limits and requests

## Next Steps

1. **Push your changes** to trigger the GitHub Action
2. **Wait for the build** to complete (check Actions tab)
3. **Update your Kubernetes deployment** using the provided script
4. **Verify the deployment** is working correctly

The GitHub Action will automatically build and push new images whenever you push changes to the repository.
