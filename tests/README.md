# pmultiqc Tests

This directory contains tests for the pmultiqc project.

## Test Data

The tests use MaxQuant output files from the PXD003133 dataset. For GitHub storage efficiency, the test data files are individually compressed with gzip in the `tests/resources/maxquant/PXD003133-22min` directory.

The compressed files include:
- proteinGroups.txt.gz (1.2MB)
- evidence.txt.gz (10MB)
- msms.txt.gz (54MB)
- msmsScans.txt.gz (11MB)
- parameters.txt.gz (829B)
- summary.txt.gz (1.0KB)

The test fixtures in `conftest.py` automatically handle reading these compressed files, so you don't need to extract them manually.

Note: Most files are under GitHub's recommended 50MB limit, but msms.txt.gz is 54MB, which is still under GitHub's hard limit of 100MB.

## Running Tests

### Command Line

To run all tests:

```bash
python -m pytest tests/
```

To run a specific test file:

```bash
python -m pytest tests/test_maxquant.py
```

To run tests with verbose output:

```bash
python -m pytest tests/test_maxquant.py -v
```

### Running Tests in PyCharm

To run tests in PyCharm:

1. **Configure pytest as the default test runner**:
   - Go to `File > Settings > Tools > Python Integrated Tools`
   - Set "Default test runner" to "pytest"
   - Click "Apply" and "OK"

2. **Run a specific test file**:
   - Right-click on the test file (e.g., `maxquant_tests.py`) in the Project view
   - Select "Run 'pytest in maxquant_tests.py'"

3. **Run a specific test method**:
   - Open the test file
   - Right-click on the test method name (e.g., `test_proteingroups_file_exists`)
   - Select "Run 'pytest for test_proteingroups_file_exists'"

4. **Run all tests**:
   - Right-click on the `tests` directory in the Project view
   - Select "Run 'pytest in tests'"

5. **Configure a permanent run configuration**:
   - Go to `Run > Edit Configurations`
   - Click the "+" button and select "Python tests > pytest"
   - Name your configuration (e.g., "All MaxQuant Tests")
   - Set the target to "Script path" and select the test file or directory
   - Set the working directory to the project root
   - Click "OK"
   - Now you can select this configuration from the run configurations dropdown

### Running Tests in GitHub CI/CD

The tests are configured to run automatically in GitHub CI/CD. There are three workflow files in the `.github/workflows/` directory:

1. **python-app.yml**: Runs functional tests for different datasets (LFQ, TMT, DIA, mzid_mzML, mzid_mgf, MaxQuant)
2. **python-package.yml**: Runs pytest tests across multiple Python versions (3.9, 3.10, 3.11)
3. **pytest.yml**: Specifically runs pytest tests with coverage reporting

When you push changes to the repository, GitHub Actions will automatically run these workflows to ensure everything is working correctly.

To view the CI/CD results:
1. Go to the GitHub repository
2. Click on the "Actions" tab
3. Select the workflow you want to view
4. Click on a specific workflow run to see details

## Test Files

### maxquant_tests.py

This file contains tests for the MaxQuant module functionality. The tests verify that:

1. MaxQuant output files exist and have the expected structure
2. Basic functionality of the maxquant module can be tested without direct imports
3. Core functions like reading files, filtering contaminants, and analyzing data work as expected

The tests use a mock approach to avoid dependency issues with the full pmultiqc package. This allows the tests to run independently and verify the expected behavior of the functions.

## Adding New Tests

When adding new tests:

1. Follow the existing pattern of creating mock tests that verify expected behavior
2. Use pytest's skip functionality when required data or columns are not present
3. Add appropriate assertions to verify the results
4. Make sure tests can run independently of the full pmultiqc package

## Test Coverage

The current tests cover the following functionality:

- File existence and structure validation
- Reading and filtering MaxQuant output files
- Protein group analysis
- Contaminant analysis
- Intensity distribution analysis
- Charge distribution analysis
- Retention time analysis
- Peptide intensity analysis

## Working with Compressed Test Data

If you need to add or update test data:

1. Compress individual files using gzip:
   ```bash
   gzip -c your_file.txt > tests/resources/maxquant/PXD003133-22min/your_file.txt.gz
   ```

2. Update the fixtures in conftest.py if needed

3. Make sure your tests handle reading compressed files using gzip:
   ```python
   with gzip.open(file_path, 'rt') as f:
       df = pd.read_csv(f, sep="\t", low_memory=False)
   ```

4. Be mindful of GitHub's file size limits:
   - Individual files should be under 100MB (hard limit)
   - Files under 50MB are recommended for better performance
   - Repositories should ideally be under 1GB total