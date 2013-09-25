# Copyright (c) 2013 Jose Cruz-Toledo

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

""" Command-line application for getting the ligand neighbourhood of an aptamer
Usage:
  $ python ligandneighbours.py

"""

import os
import sys
import argparse
from Bio.PDB import *
#parser for command-line arguments 
parser = argparse.ArgumentParser(
	description=__doc__,
	formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('-dir','--input_dir', help='a local direcory containing PDB structures for which you wish to find the ligand neighbourhood', required=True)

def main(argv):
	#parse the command-line flags.
	flags = parser.parse_args(argv[1:])
	local_dir = flags.input_dir
	#fetch a list of all the pdb files in the input directory
	filepaths = fetchPdbFilePaths(local_dir)



def fetchPdbFilePaths(alocaldir):
	rm = []
	for fn in os.listdir(alocaldir):
		fileName, fileExtension = os.path.splitext(alocaldir+'/'+fn)
		if fileExtension == '.pdb':
			rm.append(fileName+fileExtension)
	if len(rm) == 0:
		raise Exception("No pdb files found in provided folder!")
		sys.exit()

	return rm

#start it
if __name__ == '__main__':
  main(sys.argv)