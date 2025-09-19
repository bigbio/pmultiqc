#!/usr/bin/env python3
"""
ARM64 compatibility patch for pmultiqc.
This script patches pmultiqc to work on ARM64 systems where pyopenms is not available.
"""

import os
import sys
import pmultiqc


def patch_pmultiqc_for_arm64():
    """Apply ARM64 compatibility patches to pmultiqc."""
    pmultiqc_path = os.path.dirname(pmultiqc.__file__)
    
    # Patch 1: Handle missing pyopenms in all modules
    files_to_patch = [
        'modules/quantms/quantms.py',
        'modules/common/ms_io.py'
    ]
    
    for file_path in files_to_patch:
        full_path = os.path.join(pmultiqc_path, file_path)
        
        if not os.path.exists(full_path):
            print(f"Warning: {full_path} not found")
            continue
        
        with open(full_path, 'r') as f:
            content = f.read()
        
        # Check if already patched
        if 'try:' in content and 'ImportError:' in content and 'pyopenms' in content:
            print(f"{file_path} already patched for ARM64")
            continue
        
        # Apply the patch - handle different import patterns
        patches_applied = False
        
        # Pattern 1: from pyopenms import OpenMSBuildInfo, AASequence
        old_import1 = 'from pyopenms import OpenMSBuildInfo, AASequence'
        new_import1 = '''try:
    from pyopenms import OpenMSBuildInfo, AASequence
except ImportError:
    OpenMSBuildInfo = None
    AASequence = None'''
        
        if old_import1 in content:
            print(f"Patching {file_path} for ARM64 compatibility (pattern 1)")
            content = content.replace(old_import1, new_import1)
            patches_applied = True
        
        # Pattern 2: from pyopenms import IdXMLFile, MzMLFile, MSExperiment
        old_import2 = 'from pyopenms import IdXMLFile, MzMLFile, MSExperiment'
        new_import2 = '''try:
    from pyopenms import IdXMLFile, MzMLFile, MSExperiment
except ImportError:
    IdXMLFile = None
    MzMLFile = None
    MSExperiment = None'''
        
        if old_import2 in content:
            print(f"Patching {file_path} for ARM64 compatibility (pattern 2)")
            content = content.replace(old_import2, new_import2)
            patches_applied = True
        
        # Pattern 3: Handle OpenMSBuildInfo().getOpenMPMaxNumThreads() call
        old_call = 'log.info("pyopenms has: " + str(OpenMSBuildInfo().getOpenMPMaxNumThreads()) + " threads.")'
        new_call = '''try:
    if OpenMSBuildInfo is not None:
        log.info("pyopenms has: " + str(OpenMSBuildInfo().getOpenMPMaxNumThreads()) + " threads.")
    else:
        log.info("pyopenms not available - OpenMSBuildInfo is None")
except Exception as e:
    log.info(f"pyopenms not available: {e}")'''
        
        if old_call in content:
            print(f"Patching {file_path} for ARM64 compatibility (pattern 3)")
            content = content.replace(old_call, new_call)
            patches_applied = True
        
        if patches_applied:
            with open(full_path, 'w') as f:
                f.write(content)
            print(f"{file_path} patched successfully")
        else:
            print(f"No pyopenms import found in {file_path}")


if __name__ == '__main__':
    try:
        patch_pmultiqc_for_arm64()
        print("ARM64 patches applied successfully")
    except Exception as e:
        print(f"Error applying ARM64 patches: {e}")
        sys.exit(1)
