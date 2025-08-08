#!/bin/bash

# Local testing script for PMultiQC service
# This script builds and runs the container locally

set -e

# Colors for output
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m' # No Color

echo -e "${GREEN}=== PMultiQC Service Local Testing (FastAPI) ===${NC}"
echo ""

# Check if Docker is running
if ! docker info > /dev/null 2>&1; then
    echo -e "${RED}❌ Docker is not running. Please start Docker and try again.${NC}"
    exit 1
fi

# Build the Docker image
echo -e "${YELLOW}Building Docker image...${NC}"
docker build -t pmultiqc-service:latest .

if [ $? -eq 0 ]; then
    echo -e "${GREEN}✅ Image built successfully${NC}"
else
    echo -e "${RED}❌ Failed to build image${NC}"
    exit 1
fi

echo ""

# Create local directories for testing
echo -e "${YELLOW}Creating local directories...${NC}"
mkdir -p ./local_uploads ./local_outputs ./local_html_reports
echo -e "${GREEN}✅ Directories created${NC}"

echo ""

# Run the container
echo -e "${YELLOW}Starting container...${NC}"
echo "Container will be available at: http://localhost:5000"
echo "Press Ctrl+C to stop the container"
echo ""

docker run --rm \
    -p 5000:5000 \
    -v $(pwd)/local_uploads:/tmp/pmultiqc_uploads \
    -v $(pwd)/local_outputs:/tmp/pmultiqc_outputs \
    -v $(pwd)/local_html_reports:/tmp/pmultiqc_html_reports \
    -e LOG_LEVEL=DEBUG \
    -e BASE_URL=http://localhost:5000 \
    pmultiqc-service:latest

echo ""
echo -e "${GREEN}Container stopped${NC}" 