#!/usr/bin/env python3
import sys

args = sys.argv[1:]

# Write output
with open("fileset.txt", "w") as out:
    out.write('fileset = {\n')
    
    i = 0
    while i < len(args):
        input_file = args[i]
        
        # Check if next arg is a number
        if i + 1 < len(args) and args[i + 1].isdigit():
            n_files = int(args[i + 1])
            i += 2
        else:
            n_files = None
            i += 1
        
        # Extract dataset name from filename
        dataset_name = input_file.replace('.txt', '').replace('_', '')
        
        # Read paths from this file
        with open(input_file, "r") as f:
            paths = [line.strip() for line in f if line.strip()]
            if n_files is not None:
                paths = paths[:n_files]
        
        print(f"{dataset_name}: using {len(paths)} files")
        
        out.write(f'    "{dataset_name}": {{\n')
        out.write('        "files": {\n')
        
        for path in paths:
            #out.write(f'            "root://cmsxrootd.fnal.gov/{path}": "Events",\n')
            #out.write(f'            "root://cmsxrootd.fnal.gov/{path}": "KUAnalysis",\n')
            out.write(f'            "root://cms-xrd-global.cern.ch/{path}": "Events",\n')
        
        out.write('        },\n')
        out.write('        "metadata": {"is_mc": True, "is_UL": False},\n')
        out.write('    },\n')
    
    out.write('}\n')

print("Done! Created fileset.txt")
