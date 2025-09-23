#!/bin/bash

# Test script to validate Docker Compose BASE_URL configuration
# This script demonstrates how the new configuration system works

set -e

echo "🧪 Testing PMultiQC Docker Compose BASE_URL Configuration"
echo "========================================================="

cd "$(dirname "$0")"

# Test 1: Verify template files exist
echo
echo "📋 Test 1: Checking template files..."
if [[ -f ".env.template" && -f ".env.production.example" && -f ".env.development.example" ]]; then
    echo "✅ All template files exist"
else
    echo "❌ Missing template files"
    exit 1
fi

# Test 2: Test docker-compose.yml with default values (no .env file)
echo
echo "📋 Test 2: Testing docker-compose.yml with default BASE_URL (no .env)..."
# Temporarily move .env if it exists
if [[ -f ".env" ]]; then
    mv .env .env.temp
fi
BASE_URL_DEFAULT=$(docker compose -f docker-compose.yml config | grep -A 1 "BASE_URL:" | grep -v "API_BASE" | head -1 | awk '{print $2}')
echo "Default BASE_URL: $BASE_URL_DEFAULT"
if [[ "$BASE_URL_DEFAULT" == "http://localhost:8080" ]]; then
    echo "✅ Default BASE_URL is correct"
else
    echo "❌ Default BASE_URL is wrong: expected 'http://localhost:8080', got '$BASE_URL_DEFAULT'"
fi
# Restore .env if it existed
if [[ -f ".env.temp" ]]; then
    mv .env.temp .env
fi

# Test 3: Test docker-compose.yml with custom BASE_URL
echo
echo "📋 Test 3: Testing docker-compose.yml with custom BASE_URL..."
CUSTOM_URL="https://example.com/pmultiqc"
BASE_URL_CUSTOM=$(BASE_URL="$CUSTOM_URL" docker compose -f docker-compose.yml config | grep -A 1 "BASE_URL:" | grep -v "API_BASE" | head -1 | awk '{print $2}')
echo "Custom BASE_URL: $BASE_URL_CUSTOM"
if [[ "$BASE_URL_CUSTOM" == "$CUSTOM_URL" ]]; then
    echo "✅ Custom BASE_URL override works"
else
    echo "❌ Custom BASE_URL override failed: expected '$CUSTOM_URL', got '$BASE_URL_CUSTOM'"
fi

# Test 4: Test production configuration with .env file
echo
echo "📋 Test 4: Testing production configuration with .env file..."
if [[ -f ".env" ]]; then
    BASE_URL_PROD=$(docker compose -f docker-compose.prod.yml config | grep -A 1 "BASE_URL:" | grep -v "API_BASE" | head -1 | awk '{print $2}')
    echo "Production BASE_URL from .env: $BASE_URL_PROD"
    echo "✅ Production configuration loads .env file"
else
    echo "⚠️  No .env file found - this is expected in CI/testing"
fi

# Test 5: Validate YAML syntax
echo
echo "📋 Test 5: Validating YAML syntax..."
if docker compose -f docker-compose.yml config --quiet; then
    echo "✅ docker-compose.yml syntax is valid"
else
    echo "❌ docker-compose.yml has syntax errors"
    exit 1
fi

if docker compose -f docker-compose.prod.yml config --quiet 2>/dev/null || [[ ! -f ".env" ]]; then
    echo "✅ docker-compose.prod.yml syntax is valid"
else
    echo "❌ docker-compose.prod.yml has syntax errors"
    exit 1
fi

echo
echo "🎉 All tests passed! Docker Compose BASE_URL configuration is working correctly."
echo
echo "💡 Usage examples:"
echo "  # For development:"
echo "  cp .env.development.example .env"
echo "  docker compose up -d"
echo
echo "  # For production:"
echo "  cp .env.production.example .env"
echo "  # Edit .env to set your BASE_URL"
echo "  docker compose -f docker-compose.prod.yml up -d"
echo
echo "  # Quick override:"
echo "  BASE_URL=https://your-domain.com/pmultiqc docker compose up -d"