# Security Policy

## Security Improvements

This document outlines the security improvements made to the pmultiqc project.

### 1. Path Traversal Protection

**Issue**: Multiple locations in the codebase were vulnerable to path traversal attacks.

**Fixed in**:
- `pmultiqc_service/app.py`: Enhanced `serve_html_report()` function with path validation
- `pmultiqc/modules/common/file_utils.py`: Added path validation in `extract_zip()` and `extract_tar()`

**Mitigation**: 
- Validate all file paths to prevent ".." sequences and absolute paths
- Use `os.path.normpath()` and verify paths remain within expected directories
- Reject filenames containing path separators in upload endpoints

### 2. Zip Bomb Protection

**Issue**: Malicious zip files with high compression ratios could cause resource exhaustion.

**Fixed in**:
- `pmultiqc_service/app.py`: Added compression ratio checks in upload endpoints
- `pmultiqc/modules/common/file_utils.py`: Added compression ratio validation in `extract_zip()`

**Mitigation**:
- Check compression ratio (reject if > 100:1)
- Limit maximum uncompressed size (10GB default)
- Validate total uncompressed size before extraction

### 3. File Upload Security

**Issue**: No limits on file size or number of files in uploads.

**Fixed in**:
- `pmultiqc_service/app.py`: Added `MAX_FILE_SIZE` and `MAX_UPLOAD_FILES` configuration

**Mitigation**:
- Maximum file size: 10GB (configurable via `MAX_FILE_SIZE` environment variable)
- Maximum files per ZIP: 100 (configurable via `MAX_UPLOAD_FILES` environment variable)
- Stream-based file size checking to prevent memory exhaustion

### 4. Command Injection Prevention

**Issue**: User-controlled input used in subprocess execution without validation.

**Fixed in**:
- `pmultiqc_service/app.py`: Added whitelist validation for `input_type` parameter in `run_pmultiqc_with_progress()`

**Mitigation**:
- Whitelist allowed input types: `maxquant`, `quantms`, `diann`, `mzidentml`
- Validate file paths exist before processing
- Validate config file paths for quantms plugin

### 5. CORS Configuration

**Issue**: Overly permissive CORS policy allowing all origins.

**Fixed in**:
- `pmultiqc_service/app.py`: Configurable CORS origins via `ALLOWED_ORIGINS` environment variable

**Mitigation**:
- Default to localhost origins for development
- Support comma-separated list of allowed origins via environment variable
- Restrict HTTP methods to GET, POST, OPTIONS
- Restrict headers to Content-Type and Authorization

### 6. Sensitive Information Exposure

**Issue**: Redis credentials logged to console/files.

**Fixed in**:
- `pmultiqc_service/app.py`: Added security comment to prevent credential logging

**Mitigation**:
- Log only whether credentials are provided, not the actual values
- Document secure logging practices

### 7. Input Validation

**Issue**: PRIDE accession numbers not properly validated.

**Fixed in**:
- `pmultiqc_service/app.py`: Added regex validation for PRIDE accessions (PXD followed by 6 digits)

**Mitigation**:
- Validate format: `PXD\d{6}` (e.g., PXD012345)
- Applied to all PRIDE-related endpoints

### 8. UUID Validation

**Issue**: Directory cleanup script assumed 36-character names were valid UUIDs.

**Fixed in**:
- `pmultiqc_service/cleanup_disk.py`: Added proper UUID validation

**Mitigation**:
- Use `uuid.UUID()` to validate format before processing
- Skip non-UUID directories to prevent accidental deletions

## Configuration

### Environment Variables for Security

```bash
# CORS Configuration
ALLOWED_ORIGINS="https://example.com,https://app.example.com"

# File Upload Limits
MAX_FILE_SIZE=10737418240  # 10GB in bytes
MAX_UPLOAD_FILES=100       # Maximum files in a ZIP

# Redis Configuration (never log these)
REDIS_PASSWORD=your_secure_password
REDIS_USERNAME=your_username
```

## Best Practices

1. **Always set ALLOWED_ORIGINS in production** - Never use the default wildcard CORS policy in production
2. **Use HTTPS** - Always deploy behind HTTPS in production
3. **Limit file sizes** - Adjust `MAX_FILE_SIZE` based on your infrastructure
4. **Monitor logs** - Regularly review logs for suspicious activity
5. **Keep dependencies updated** - Regularly update Python packages for security patches
6. **Use strong Redis passwords** - If Redis is exposed, use strong authentication

## Reporting Security Issues

If you discover a security vulnerability, please email the maintainers directly rather than opening a public issue:

- Email: [security contact]
- Include: Detailed description, steps to reproduce, potential impact

## Security Testing

Run security scans before deployment:

```bash
# Run CodeQL security analysis
# (Configured in GitHub Actions)

# Check for known vulnerabilities in dependencies
pip install safety
safety check

# Run Bandit for Python security issues
pip install bandit
bandit -r pmultiqc pmultiqc_service
```

## Version History

- **v0.0.40+**: Security improvements implemented
  - Path traversal protection
  - Zip bomb detection
  - File upload limits
  - CORS configuration
  - Input validation
  - Command injection prevention
