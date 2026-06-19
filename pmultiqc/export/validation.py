from mzqc import MZQCFile
from mzqc.SyntaxCheck import SyntaxCheck
from mzqc.SemanticCheck import SemanticCheck

file_path = "/home/timo/Desktop/BA/data/mq_data/jsons/mq_data_quality.mzQC"

with open(file_path, "r") as f:
    content = f.read()

# Parse
mzqc_obj = MZQCFile.JsonSerialisable.from_json(content)

#  syntax
syntax = SyntaxCheck()
syntax.validate(mzqc_obj)

#  semantics
semantic = SemanticCheck(mzqc_obj)
semantic.validate()

print("Validation completed.\n" + "="*30)

# Check for syntax errors explicitly
# The python-mzqc library stores syntax errors in a dictionary or list attribute (usually .errors)
if hasattr(syntax, 'errors') and syntax.errors:
    print("Syntax issues found:")
    print(syntax.errors)
else:
    # If it has an 'errors' attribute but it's empty, or no errors are flagged:
    print("Syntax check: PASSED (0 issues found)")

print("-" * 30)

# Check semantic issues
# Since this printed a dictionary with empty lists before, we can format it nicely
has_semantic_issues = any(len(errors) > 0 for errors in semantic.values()) if isinstance(semantic, dict) else False

if has_semantic_issues:
    print("Semantic issues found:")
    print(semantic)
else:
    print("Semantic check: PASSED (0 issues found)")