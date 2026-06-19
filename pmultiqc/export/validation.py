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

print("Validation completed.")
print("Syntax issues:", syntax)
print("Semantic issues:", semantic)