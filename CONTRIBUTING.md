# Contributing to pmultiqc

Thank you for your interest in contributing to pmultiqc! This guide will help you get started.

## Quick Start

1. **Fork and clone** the repository
2. **Install in development mode**: `pip install -e .`
3. **Install dev dependencies**: `pip install pytest black isort`
4. **Create a feature branch**: `git checkout -b feature/your-feature`
5. **Make your changes** and test them
6. **Submit a pull request**

## Reporting Issues

We have issue templates for different types of problems:

- **Bug Reports**: Crashes, incorrect metrics, unexpected behavior
- **Metric Requests**: New proteomics quality control metrics (we encourage these!)
- **Feature Requests**: New visualizations, data format support, functionality
- **Service Issues**: Problems with the PRIDE web service
- **Suggestions**: General improvements and ideas

## Development

### Code Style
- Follow PEP 8 guidelines
- Use type hints where appropriate
- Write docstrings for functions and classes
- Format code with `black` and `isort`

### Testing
```bash
# Run tests
pytest tests/

# Test with sample data
cd tests && multiqc resources/LFQ -o ./test_output
```

### Adding New Metrics

1. **Choose the right module** (quantms, maxquant, proteobench, core)
2. **Implement the calculation** function
3. **Add visualization** function
4. **Write tests** for your metric
5. **Update documentation**

### Adding New Visualizations

- Use consistent styling with matplotlib/seaborn
- Make plots interactive with plotly when possible
- Add plots to the MultiQC report using `add_section()`

## Pull Requests

When submitting a PR:

- [ ] Code follows style guidelines
- [ ] Tests are added/updated
- [ ] Documentation is updated
- [ ] All tests pass
- [ ] PR description is complete

## Community

- **Be respectful** and inclusive
- **Be constructive** in feedback
- **Be patient** with newcomers
- **Help others** when you can

## Getting Help

- **GitHub Discussions**: For questions and ideas
- **Issue Tracker**: For bugs and feature requests
- **Email**: ypriverol@gmail.com for urgent issues

## Good First Issues

Look for issues labeled `good-first-issue` or `help-wanted`:
- Documentation improvements
- Adding examples
- Fixing typos
- Improving error messages
- Adding tests
- Code cleanup

---

Thank you for contributing to pmultiqc! Your contributions help make proteomics analysis more accessible and reliable for the entire community.
