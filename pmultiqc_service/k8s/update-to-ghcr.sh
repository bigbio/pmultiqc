#!/bin/bash

# Script to update PMultiQC deployment to use GitHub Container Registry
# This script updates the deployment to use the new GHCR image

set -e

# Colors for output
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m' # No Color

echo -e "${GREEN}=== Updating PMultiQC to use GitHub Container Registry ===${NC}"
echo ""

# Check if kubectl is available
if ! command -v kubectl &> /dev/null; then
    echo -e "${RED}❌ kubectl is not installed or not in PATH${NC}"
    exit 1
fi

# Set the kubeconfig
KUBECONFIG_FILE="$HOME/.kube/config_bckp"
if [ ! -f "$KUBECONFIG_FILE" ]; then
    echo -e "${RED}❌ Kubeconfig file not found: $KUBECONFIG_FILE${NC}"
    exit 1
fi

echo -e "${YELLOW}Using kubeconfig: $KUBECONFIG_FILE${NC}"

# Check current deployment
echo -e "${YELLOW}Current deployment status:${NC}"
kubectl --kubeconfig "$KUBECONFIG_FILE" get deployment pmultiqc-service -n pmultiqc -o jsonpath='{.spec.template.spec.containers[0].image}' 2>/dev/null || echo "Deployment not found"
echo ""

# Apply the updated deployment
echo -e "${YELLOW}Applying updated deployment with GHCR image...${NC}"
kubectl --kubeconfig "$KUBECONFIG_FILE" apply -f deployment.yaml -n pmultiqc

if [ $? -eq 0 ]; then
    echo -e "${GREEN}✅ Deployment updated successfully${NC}"
else
    echo -e "${RED}❌ Failed to update deployment${NC}"
    exit 1
fi

echo ""

# Wait for rollout to complete
echo -e "${YELLOW}Waiting for rollout to complete...${NC}"
kubectl --kubeconfig "$KUBECONFIG_FILE" rollout status deployment/pmultiqc-service -n pmultiqc --timeout=300s

if [ $? -eq 0 ]; then
    echo -e "${GREEN}✅ Rollout completed successfully${NC}"
else
    echo -e "${RED}❌ Rollout failed or timed out${NC}"
    exit 1
fi

echo ""

# Show final status
echo -e "${YELLOW}Final deployment status:${NC}"
kubectl --kubeconfig "$KUBECONFIG_FILE" get pods -n pmultiqc
echo ""

echo -e "${YELLOW}New image being used:${NC}"
kubectl --kubeconfig "$KUBECONFIG_FILE" get deployment pmultiqc-service -n pmultiqc -o jsonpath='{.spec.template.spec.containers[0].image}'
echo ""

echo -e "${GREEN}🎉 PMultiQC has been successfully updated to use GitHub Container Registry!${NC}"
echo -e "${YELLOW}Image: ghcr.io/bigbio/pmultiqc:latest${NC}"
