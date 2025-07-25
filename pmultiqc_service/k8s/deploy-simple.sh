#!/bin/bash

# PMultiQC Service Simple Kubernetes Deployment Script

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Configuration
NAMESPACE="pmultiqc"

echo -e "${GREEN}PMultiQC Service Kubernetes Deployment${NC}"
echo "=============================================="

# Function to check kubectl connection
check_kubectl_connection() {
    echo -e "${GREEN}Checking kubectl connection...${NC}"
    
    if ! kubectl cluster-info &> /dev/null; then
        echo -e "${RED}Error: Cannot connect to Kubernetes cluster${NC}"
        echo "Please check:"
        echo "1. kubectl is installed"
        echo "2. kubectl is configured to connect to your cluster"
        echo "3. Your cluster is running"
        echo ""
        echo "Try running: kubectl cluster-info"
        exit 1
    fi
    
    echo -e "${GREEN}✓ Connected to Kubernetes cluster${NC}"
    kubectl cluster-info
}

# Function to deploy resources individually
deploy_resources() {
    echo -e "${GREEN}Deploying PMultiQC resources...${NC}"
    
    # Get the directory where this script is located
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    
    # Create namespace
    echo "Creating namespace..."
    kubectl apply -f "$SCRIPT_DIR/namespace.yaml"
    
    # Create configmap
    echo "Creating configmap..."
    kubectl apply -f "$SCRIPT_DIR/configmap.yaml"
    
    # Create PVCs
    echo "Creating persistent volume claims..."
    kubectl apply -f "$SCRIPT_DIR/persistent-volume-claim.yaml"
    
    # Create deployment
    echo "Creating deployment..."
    kubectl apply -f "$SCRIPT_DIR/deployment.yaml"
    
    # Create service
    echo "Creating service..."
    kubectl apply -f "$SCRIPT_DIR/service.yaml"
    
    # Create HPA
    echo "Creating horizontal pod autoscaler..."
    kubectl apply -f "$SCRIPT_DIR/hpa.yaml"
    
    # Create ingress (optional - comment out if you don't have ingress controller)
    echo "Creating ingress..."
    kubectl apply -f "$SCRIPT_DIR/ingress.yaml"
    
    echo -e "${GREEN}✓ All resources deployed${NC}"
}

# Function to check deployment status
check_status() {
    echo -e "${GREEN}Checking deployment status...${NC}"
    
    echo "Namespace:"
    kubectl get namespace $NAMESPACE
    
    echo ""
    echo "Pods:"
    kubectl get pods -n $NAMESPACE
    
    echo ""
    echo "Services:"
    kubectl get svc -n $NAMESPACE
    
    echo ""
    echo "PVCs:"
    kubectl get pvc -n $NAMESPACE
    
    echo ""
    echo "HPA:"
    kubectl get hpa -n $NAMESPACE
    
    echo ""
    echo "Ingress:"
    kubectl get ingress -n $NAMESPACE
}

# Function to show logs
show_logs() {
    echo -e "${GREEN}Showing logs...${NC}"
    kubectl logs -f deployment/pmultiqc-service -n $NAMESPACE
}

# Function to delete deployment
delete_deployment() {
    echo -e "${YELLOW}Deleting PMultiQC service...${NC}"
    
    # Get the directory where this script is located
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    
    kubectl delete -f "$SCRIPT_DIR/ingress.yaml" --ignore-not-found=true
    kubectl delete -f "$SCRIPT_DIR/hpa.yaml" --ignore-not-found=true
    kubectl delete -f "$SCRIPT_DIR/service.yaml" --ignore-not-found=true
    kubectl delete -f "$SCRIPT_DIR/deployment.yaml" --ignore-not-found=true
    kubectl delete -f "$SCRIPT_DIR/persistent-volume-claim.yaml" --ignore-not-found=true
    kubectl delete -f "$SCRIPT_DIR/configmap.yaml" --ignore-not-found=true
    kubectl delete -f "$SCRIPT_DIR/namespace.yaml" --ignore-not-found=true
    
    echo -e "${GREEN}✓ Deletion completed${NC}"
}

# Function to port forward for local access
port_forward() {
    echo -e "${GREEN}Setting up port forward...${NC}"
    echo "Access the service at: http://localhost:8080"
    echo "Press Ctrl+C to stop"
    kubectl port-forward svc/pmultiqc-service 8080:80 -n $NAMESPACE
}

# Function to show usage
usage() {
    echo "Usage: $0 [COMMAND]"
    echo ""
    echo "Commands:"
    echo "  check     Check kubectl connection"
    echo "  deploy    Deploy all resources"
    echo "  status    Check deployment status"
    echo "  logs      Show application logs"
    echo "  delete    Delete all resources"
    echo "  forward   Port forward for local access"
    echo "  all       Check connection, deploy, and show status"
    echo "  help      Show this help message"
}

# Main script logic
case "${1:-help}" in
    check)
        check_kubectl_connection
        ;;
    deploy)
        check_kubectl_connection
        deploy_resources
        ;;
    status)
        check_status
        ;;
    logs)
        show_logs
        ;;
    delete)
        delete_deployment
        ;;
    forward)
        port_forward
        ;;
    all)
        check_kubectl_connection
        deploy_resources
        echo ""
        echo "Waiting for pods to be ready..."
        sleep 10
        check_status
        ;;
    help|*)
        usage
        ;;
esac 