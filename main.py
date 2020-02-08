import os
import argparse

# Parse the arguments to get the input file.
parser = argparse.ArgumentParser()
parser.add_argument("-i", help="An input file that tells me what to do")
args = parser.parse_args()

if args.i is None:
    raise SystemExit("The input file is a required argument. Specify it via the flag `-i`")
