from __future__ import annotations
from argparse import ArgumentParser, Namespace
import sys
from ..services.clean_service import run_clean_service

def add_clean_parser(parser: ArgumentParser) -> None:
    parser.add_argument("-i","--input_path", required=True, help="Path to input CIF/PDB file.")
    parser.add_argument("-o","--output_dir", required=True, help="Path to a directory for outputting cleaned CIF, PDB, and FASTA files and a JSON report.")
    parser.add_argument("--no_add_H",action="store_false",dest="add_H",help="Disable adding hydrogens using OpenMM (default: enabled).")
    parser.set_defaults(add_H=True)
    parser.add_argument("--pH",type=float,default=7.0,help="pH value for hydrogen addition (default: 7.0).")

    parser.set_defaults(func=run_clean)

def run_clean(args: Namespace) -> None:
    success = run_clean_service(input_path=args.input_path, output_dir=args.output_dir, add_H=args.add_H, pH=args.pH)
    if not success:
        sys.exit(1)


