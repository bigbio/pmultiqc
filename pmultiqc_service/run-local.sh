#!/bin/bash

# Local testing script for pmultiqc service with Redis
# This script starts Redis and runs the container locally

set -e

# Colors for output
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m' # No Color

echo -e "${GREEN}=== pmultiqc Service Local Testing with Redis ===${NC}"
echo ""

# Check if Docker is running
if ! docker info > /dev/null 2>&1; then
    echo -e "${RED}❌ Docker is not running. Please start Docker and try again.${NC}"
    exit 1
fi

# Check if Redis container exists (running or stopped)
if docker ps -a --format "table {{.Names}}" | grep -q "pmultiqc-redis"; then
    echo -e "${YELLOW}Redis container exists, stopping and removing it...${NC}"
    docker stop pmultiqc-redis > /dev/null 2>&1 || true
    docker rm pmultiqc-redis > /dev/null 2>&1 || true
    echo -e "${GREEN}✅ Existing Redis container cleaned up${NC}"
fi

# Start Redis container
echo -e "${YELLOW}Starting Redis container...${NC}"
REDIS_PASSWORD=${REDIS_PASSWORD:-"pmultiqc123"}
docker run -d \
    --name pmultiqc-redis \
    -p 6379:6379 \
    -e REDIS_PASSWORD=$REDIS_PASSWORD \
    redis:7-alpine redis-server --requirepass $REDIS_PASSWORD

if [ $? -eq 0 ]; then
    echo -e "${GREEN}✅ Redis container started successfully${NC}"
else
    echo -e "${RED}❌ Failed to start Redis container${NC}"
    exit 1
fi

# Wait for Redis to be ready
echo -e "${YELLOW}Waiting for Redis to be ready...${NC}"
sleep 3

# Test Redis connection
if docker exec pmultiqc-redis redis-cli -a $REDIS_PASSWORD ping | grep -q "PONG"; then
    echo -e "${GREEN}✅ Redis is ready${NC}"
else
    echo -e "${RED}❌ Redis is not responding${NC}"
    exit 1
fi

echo ""

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

# Run the container with Redis connection
echo -e "${YELLOW}Starting pmultiqc service container...${NC}"
echo "Service will be available at: http://localhost:5000"
echo "Redis is available at: localhost:6379"
echo "Press Ctrl+C to stop the containers"
echo ""

docker run --rm \
    --name pmultiqc-service \
    --link pmultiqc-redis:redis \
    -p 5000:5000 \
    -v $(pwd)/local_uploads:/tmp/pmultiqc_uploads \
    -v $(pwd)/local_outputs:/tmp/pmultiqc_outputs \
    -v $(pwd)/local_html_reports:/tmp/pmultiqc_html_reports \
    -e LOG_LEVEL=DEBUG \
    -e BASE_URL=http://localhost:5000 \
    -e REDIS_URL=redis://redis:6379 \
    -e REDIS_PASSWORD=$REDIS_PASSWORD \
    pmultiqc-service:latest

echo ""
echo -e "${YELLOW}Stopping Redis container...${NC}"
docker stop pmultiqc-redis > /dev/null 2>&1 || true
docker rm pmultiqc-redis > /dev/null 2>&1 || true
echo -e "${GREEN}✅ All containers stopped${NC}" 