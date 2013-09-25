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
import urllib2
import csv
import re
from Bio.PDB import *
from collections import defaultdict
#parser for command-line arguments 
parser = argparse.ArgumentParser(
	description=__doc__,
	formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('-dir','--input_dir', help='a local direcory containing PDB structures for which you wish to find the ligand neighbourhood', required=True)
parser.add_argument('-out', '--output_file', help='the file where the output will be stored as CSV', required=True)
parser.add_argument('--radius', nargs='?', const=5.0, type=float, default=5.0)
ligand_map_url = 'https://docs.google.com/spreadsheet/pub?key=0AnGgKfZdJasrdC00bUxHcVRXaFloSnJYb3VmYkwyVnc&single=true&gid=0&output=csv'
base_uri = 'http://bio2rdf.org'
def main(argv):
	#parse the command-line flags.
	flags = parser.parse_args(argv[1:])
	local_dir = flags.input_dir
	output_dir = flags.output_file
	radius = flags.radius
	#fetch a list of all the pdb files in the input directory
	filepaths = fetchPdbFilePaths(local_dir)
	#fetch the ligand list
	ligands = fetchLigandList(ligand_map_url)
	findNeighbours(filepaths, ligands, radius)

def findNeighbours(someFilePaths, someLigands, aRadius):
	d = defaultdict(list)
	#iterate over the pdb file paths
	for ap in someFilePaths:
		fn, fe = os.path.splitext(ap)
		anId = fn.rsplit('/')[-1]
		match = re.match('^\w{4}$', anId)
		if match:
			p = PDBParser(PERMISSIVE=1, QUIET=1)
			structure = p.get_structure(anId, ap)
			models = structure.get_list()
			#iterate over the models
			for aModel in models:
				chains = aModel.get_list()
				#get all the atoms in this model
				model_atoms = Selection.unfold_entities(aModel, 'A')
				#create a neighbor search
				ns = NeighborSearch(model_atoms)
				#search the chains for any ligands
				for aChain in chains:
					residues = aChain.get_list()
					for aR in residues:
						if findLigandIdInListOfLigands(someLigands, aR.get_resname()) != None:
							#get the atoms of this residue
							atom_list = Selection.unfold_entities(aR, 'A')
							#pick a center
							center = atom_list[0].get_coord()
							neighbors = ns.search(center, aRadius)
							residue_list = Selection.unfold_entities(neighbors,'R')
							for aResidue in residue_list:
								#create rdf statements
								fi= aResidue.get_full_id()
								res_pdbid = fi[0]
								res_model_num = fi[1]
								res_chain = fi[2]
								res_position = fi[3][1]
								res_uri = base_uri+'/pdb_resource:'+res_pdbid+'/chemicalComponent_'+res_chain+str(res_position)
								if res_pdbid in d:
									d[res_pdbid].append(res_uri)
								else:
									d[res_pdbid] = [res_uri]			
		else :
			raise Exception ("Not a valid PDBID found!")

	print d

def fetchLigandList(aUrl):
	rm = []
	r = urllib2.urlopen(ligand_map_url)
	reader = csv.reader(r)
	rownum = 0;
	for row in reader:
		if rownum != 0:
			ligand_check = row[6]
			if ligand_check == 'Y':
				rm.append(row)
		rownum += 1
	if len(rm) == 0:
		raise Exception ("No valid ligands found in list !")
		sys.exit()
	return rm

# search aListOfLigands for a given ligand id 
# aListOfLigands is derived from the google speadsheet
# returns the ligid or None if not found
def findLigandIdInListOfLigands(aListOfLigands, aLigId):
	for al in aListOfLigands:
		anid = al[0]
		if anid == aLigId:
			return anid
	return None


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