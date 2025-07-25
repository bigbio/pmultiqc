#!/bin/bash

# PMultiQC Service Kubernetes Deployment Script

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Configuration
NAMESPACE="pmultiqc"
IMAGE_NAME="pmultiqc-service"
IMAGE_TAG="latest"

echo -e "${GREEN}PMultiQC Service Kubernetes Deployment${NC}"
echo "=============================================="

# Function to check if kubectl is available
check_kubectl() {
    if ! command -v kubectl &> /dev/null; then
        echo -e "${RED}Error: kubectl is not installed${NC}"
        exit 1
    fi
}

# Function to check if namespace exists
check_namespace() {
    if kubectl get namespace $NAMESPACE &> /dev/null; then
        echo -e "${YELLOW}Namespace $NAMESPACE already exists${NC}"
    else
        echo -e "${GREEN}Creating namespace $NAMESPACE${NC}"
        kubectl apply -f k8s/namespace.yaml
    fi
}

# Function to build and push Docker image
build_image() {
    echo -e "${GREEN}Building Docker image...${NC}"
    
    # Check if Docker is available
    if ! command -v docker &> /dev/null; then
        echo -e "${RED}Error: Docker is not installed${NC}"
        exit 1
    fi
    
    # Build the image
    docker build -t $IMAGE_NAME:$IMAGE_TAG .
    
    echo -e "${GREEN}Image built successfully${NC}"
}

# Function to deploy to Kubernetes
deploy() {
    echo -e "${GREEN}Deploying to Kubernetes...${NC}"
    
    # Apply all resources
    kubectl apply -k k8s/
    
    echo -e "${GREEN}Deployment completed${NC}"
}

# Function to check deployment status
check_status() {
    echo -e "${GREEN}Checking deployment status...${NC}"
    
    # Wait for pods to be ready
    kubectl wait --for=condition=ready pod -l app=pmultiqc-service -n $NAMESPACE --timeout=300s
    
    # Show pod status
    kubectl get pods -n $NAMESPACE
    
    # Show service status
    kubectl get svc -n $NAMESPACE
    
    # Show ingress status
    kubectl get ingress -n $NAMESPACE
}

# Function to show logs
show_logs() {
    echo -e "${GREEN}Showing logs...${NC}"
    kubectl logs -f deployment/pmultiqc-service -n $NAMESPACE
}

# Function to delete deployment
delete() {
    echo -e "${YELLOW}Deleting PMultiQC service...${NC}"
    kubectl delete -k k8s/
    echo -e "${GREEN}Deletion completed${NC}"
}

# Function to show usage
usage() {
    echo "Usage: $0 [COMMAND]"
    echo ""
    echo "Commands:"
    echo "  build     Build Docker image"
    echo "  deploy    Deploy to Kubernetes"
    echo "  status    Check deployment status"
    echo "  logs      Show application logs"
    echo "  delete    Delete deployment"
    echo "  all       Build, deploy, and check status"
    echo "  help      Show this help message"
}

# Main script logic
case "${1:-help}" in
    build)
        check_kubectl
        build_image
        ;;
    deploy)
        check_kubectl
        check_namespace
        deploy
        ;;
    status)
        check_kubectl
        check_status
        ;;
    logs)
        check_kubectl
        show_logs
        ;;
    delete)
        check_kubectl
        delete
        ;;
    all)
        check_kubectl
        check_namespace
        build_image
        deploy
        check_status
        ;;
    help|*)
        usage
        ;;
esac 