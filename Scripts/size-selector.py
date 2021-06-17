#!/usr/bin/env python3

import sys, argparse
from pathlib import Path
from Bio import SeqIO

def main():
    print("""Size selector""")
    parser = argparse.ArgumentParser()
    parser.add_argument("--in", type=str,
                    action="store", dest="input",
                    help="Path to sequence file")
    parser.add_argument("--minlen", type=int,
                    action="store", dest="minlen",
                    help="Minimum sequence size (inclusive)")
    parser.add_argument("--out", type=str,
                    action="store", dest="output",
                    help="Output path. Will create directory if it doesn't exist. Default='./<input>.min<minsize>.<filetype>'")
    
    # Parse arguments
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args()
    
    fasta = Path(args.input)
    seqtype = fasta.suffix.split(".")[-1]
    longseqs = []
    i=0
    with open(fasta, "r") as handle:
        for record in SeqIO.parse(handle, seqtype):
            if len(record.seq) >= int(args.minlen):
                longseqs.append(record)
            i+=1

    if args.output is not None:
        output = Path(args.output)
    else:
        output = Path("{}.min{}{}".format(fasta.stem, args.minlen, fasta.suffix))
    
    if not output.parent.exists():
        output.parent.mkdir(parents=True, exist_ok=True)
    
    SeqIO.write(longseqs, output, seqtype)
    print("{} original sequences\n{} less than {}\n{} sequences written to {}".format(i, i-len(longseqs), args.minlen, len(longseqs), output))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if __name__ == "__main__":
    main()
