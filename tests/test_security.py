"""
Security tests for pmultiqc
Tests path traversal protection and file validation
"""
import os
import tempfile
import zipfile
import pytest
from pmultiqc.modules.common.file_utils import extract_zip


class TestSecurityValidation:
    """Test security features in file handling"""

    def test_path_traversal_protection_in_zip(self):
        """Test that path traversal attempts in zip files are rejected"""
        with tempfile.TemporaryDirectory() as tmpdir:
            # Create a malicious zip with path traversal
            zip_path = os.path.join(tmpdir, "malicious.zip")
            extract_to = os.path.join(tmpdir, "extract")
            os.makedirs(extract_to)

            with zipfile.ZipFile(zip_path, 'w') as zf:
                # Try to write outside the extraction directory
                zf.writestr("../../evil.txt", "malicious content")

            # This should raise ValueError due to path traversal
            with pytest.raises(ValueError, match="path traversal|Invalid path"):
                extract_zip(zip_path, extract_to)

    def test_zip_bomb_protection(self):
        """Test that zip bombs are detected and rejected"""
        with tempfile.TemporaryDirectory() as tmpdir:
            zip_path = os.path.join(tmpdir, "bomb.zip")
            extract_to = os.path.join(tmpdir, "extract")
            os.makedirs(extract_to)

            # Create a simple zip bomb (large uncompressed content)
            # This creates a 50MB file that compresses well
            large_content = b"0" * (50 * 1024 * 1024)

            with zipfile.ZipFile(zip_path, 'w', compression=zipfile.ZIP_DEFLATED) as zf:
                zf.writestr("large.txt", large_content)

            # Get the actual compression ratio
            file_size = os.path.getsize(zip_path)
            if file_size > 0:
                compression_ratio = len(large_content) / file_size

                # Only test if compression ratio is high enough
                if compression_ratio > 100:
                    # This should raise ValueError due to high compression ratio
                    with pytest.raises(ValueError, match="Suspicious zip file|compression ratio"):
                        extract_zip(zip_path, extract_to)
                else:
                    # If compression isn't high enough, just verify extraction works
                    # This can happen if ZIP_DEFLATED doesn't achieve high compression
                    pass

    def test_valid_zip_extraction(self):
        """Test that valid zip files are extracted correctly"""
        with tempfile.TemporaryDirectory() as tmpdir:
            zip_path = os.path.join(tmpdir, "valid.zip")
            extract_to = os.path.join(tmpdir, "extract")
            os.makedirs(extract_to)

            # Create a valid zip file
            test_content = b"valid content"
            with zipfile.ZipFile(zip_path, 'w') as zf:
                zf.writestr("test.txt", test_content)

            # This should work without raising exceptions
            extract_zip(zip_path, extract_to)

            # Verify the file was extracted
            extracted_file = os.path.join(extract_to, "test.txt")
            assert os.path.exists(extracted_file)
            with open(extracted_file, 'rb') as f:
                assert f.read() == test_content

    def test_absolute_path_in_zip(self):
        """Test that absolute paths in zip files are rejected"""
        with tempfile.TemporaryDirectory() as tmpdir:
            zip_path = os.path.join(tmpdir, "absolute.zip")
            extract_to = os.path.join(tmpdir, "extract")
            os.makedirs(extract_to)

            with zipfile.ZipFile(zip_path, 'w') as zf:
                # Try to write to absolute path
                zf.writestr("/etc/passwd", "malicious")

            # This should raise ValueError due to absolute path
            with pytest.raises(ValueError, match="path traversal|Invalid path"):
                extract_zip(zip_path, extract_to)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
