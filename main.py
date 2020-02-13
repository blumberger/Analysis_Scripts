import os
import argparse

from src.parsing import input_file as inp_file_utils
#from IPython import embed
#c = get_config()
#c.InteractiveShellEmbed.colors = "Linux"

# Parse the arguments to get the input file.
parser = argparse.ArgumentParser()
parser.add_argument("-i", help="An input file that tells me what to do")
args = parser.parse_args()

# Check the input file exists
if args.i is None:
    raise SystemExit("The input file is a required argument. Specify it via the flag `-i`")

input_filepath = os.path.join(os.getcwd(), args.i)

INP_File = inp_file_utils.INP_File(input_filepath)

#embed(config=c)
