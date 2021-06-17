#!/usr/bin/env python
# coding: utf-8

# Core
import re, glob, sys, string, argparse
import pandas as pd
from pathlib import Path
import numpy as np
from datetime import datetime
# Bio and tqdm
from tqdm import tqdm
# Bio
from Bio import SeqIO

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def main():
    print("""
   _____ __  __ ______ 
  / ____|  \/  |  ____|
 | |  __| \  / | |__   
 | | |_ | |\/| |  __|  
 | |__| | |  | | |____ 
  \_____|_|  |_|______|
                       
Gene Model Extractor
""")
    parser = argparse.ArgumentParser()
    parser.add_argument("--models", type=str,
                    action="store", dest="models_path",
                    help="Path to file with list of models (one per line)")
    parser.add_argument("--fasta", type=str,
                    action="store", dest="fasta_path",
                    help="Path to fasta")
    parser.add_argument("--output", type=str, default="./",
                    action="store", dest="output_path",
                    help="Output path. Will create directory if it doesn't exist. Default='./output.fasta'")

    # Parse arguments
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args()

    models_path = Path(args.models_path)
    fasta_path = Path(args.fasta_path)
    output_path = Path(args.output_path)

    if not models_path.is_file():
        print("Could not find file {}".format(models_path))
    if not fasta_path.is_file():
        print("Could not find file {}".format(fasta_path))
    if not output_path.parent.exists():
        output_path.parent.mkdir(parents=True, exist_ok=True)
        print("Created output directory {}".format(str(output_path.parent.absolute())))

    print("Model labels path: {}".format(models_path))
    print("Master fasta path: {}".format(fasta_path))
    print("Output path: {}".format(output_path))

    with open(models_path, "r") as f:
        models = [re.sub("\n", "", l) for l in f.readlines()]
    records = list(SeqIO.parse(fasta_path, format="fasta"))
    print("{} contains {} records".format(fasta_path.name, len(records)))

    with open(output_path, "w+") as handle:
        model_records = []
        names = [] # Maintain list of unique names
        for m in tqdm(sorted(set(models))):
            for r in tqdm(records, leave=False):
                # Candidate records have matching description
                if m in r.description:
                    if r.name in names:
                        break # Skip if already found
                    else:
                        names.append(r.name)
                        model_records.append(r)
                        break

        SeqIO.write(sequences=model_records, handle=handle, format="fasta")
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if __name__ == "__main__":
    start = datetime.now()
    main()
    end = datetime.now()
    print("[{}] Time elapsed: {}".format(datetime.now().strftime('%d %b %Y %H:%M:%S'), end-start))
