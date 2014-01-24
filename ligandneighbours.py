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

""" 
Command-line application for creating NTriples representation of a ligand neighbourhood of a pdb file.
This program takes as input a directory of pdb files and searches those structures for  any residues found in a known dictionary of ligands:
(https://docs.google.com/spreadsheet/pub?key=0AnGgKfZdJasrdC00bUxHcVRXaFloSnJYb3VmYkwyVnc&single=true&gid=0&output=csv). If a ligand is found
in a structure an NTriples file is generated that includes details about the neighbourhood members and 

Usage:
  $ python ligandneighbours.py -dir /path/to/loca/dir/with/pdb/files --radius 4.8 --out /path/to/output 

"""

import os
import sys
import argparse
import urllib2
import csv
import re
import hashlib
import random

from Bio.PDB import *
from collections import defaultdict

#parser for command-line arguments 
parser = argparse.ArgumentParser(
	description=__doc__,
	formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('-dir','--input_dir', help='a local direcory containing PDB structures (in the .pdb format) for which you wish to find the ligand neighbourhood', required=True)
parser.add_argument('-out', '--output_file', help='the file where the output will be stored as CSV', required=True)
parser.add_argument('--radius', nargs='?', const=5.0, type=float, default=5.0)
pdb_to_ligand_list_url = 'https://docs.google.com/spreadsheet/pub?key=0AnGgKfZdJasrdC00bUxHcVRXaFloSnJYb3VmYkwyVnc&single=true&gid=0&output=csv'
base_uri = 'http://bio2rdf.org'
rdfs = 'http://www.w3.org/2000/01/rdf-schema#'
rdf = 'http://www.w3.org/1999/02/22-rdf-syntax-ns#'

def main(argv):
	#parse the command-line flags.
	flags = parser.parse_args(argv[1:])
	local_dir = flags.input_dir
	output_dir = flags.output_file
	radius = flags.radius
	#fetch a list of all the pdb files in the input directory
	filepaths = fetchPdbFilePaths(local_dir)
	#fetch the ligand list 
	ligands = fetchLigandList(pdb_to_ligand_list_url)
	for fp in filepaths:
		#get the file name and extension of the pdb file
		fn, fe = os.path.splitext(fp)
		pdbId = fn.rsplit('/')[-1]
		# dict of ligands (PDB.Residue) to residues (PDB.Residue) 
		ln = findNeighbours(fp, ligands, radius)
		#now we can generate a list of uris for each ligand neighbor
		luri = makeURIHood(ln)
		hoodNTriples = makeHoodNTriplesAnnotation(luri, pdbId, radius)
		#write an N3 file as output
		writeN3Hood(hoodNTriples, pdbId, output_dir)

#Creates a ligand neighborhood
def makeHoodNTriplesAnnotation(ligand_uri_dict, aPdbId, aRadius):
	rm = ''
	#make a hood uri
	hood_uri = base_uri+'/lighood_resource:'+hashlib.sha224(str(aPdbId)+str(aRadius)+str(random.random())).hexdigest()
	#type the hood
	rm += "<"+hood_uri+"> <"+rdf+"type> <"+base_uri+"/lighood_vocabulary:ligand_neighbourhood> .\n"
	rm += "<"+base_uri+"/lighood_vocabulary:ligand_neighbourhood> <http://www.w3.org/1999/02/22-rdf-syntax-ns#type> <http://www.w3.org/2000/01/rdf-schema#Class>  ."
	#link it to the pdb structure
	rm += "<"+base_uri+"/pdb:"+aPdbId+"> <"+base_uri+"/lighood_vocabulary:has_neighborhood> <"+hood_uri+"> .\n"
	#add the radius 
	radius_uri = base_uri+'/lighood_resource:'+hashlib.sha224(str(aRadius)+str(random.random())).hexdigest()
	rm += "<"+hood_uri+"> <"+base_uri+"/lighood_vocabulary:has_attribute> <"+radius_uri+">. \n"
	rm += "<"+radius_uri+"> <"+rdf+"type> <"+base_uri+"/lighood_vocabulary:radius> .\n"
	rm += "<"+radius_uri+"> <"+base_uri+"/lighood_vocabulary:has_value> \""+str(aRadius)+"\". \n"
	for (ligand_uri, res_uri) in ligand_uri_dict.items():
		#add ligand 
		rm += "<"+hood_uri+"> <"+base_uri+"/lighood_vocabulary:has_member> <"+ligand_uri+"> .\n"
		#type the ligand
		rm += "<"+ligand_uri+"> <"+rdf+"type> <"+base_uri+"/lighood_vocabulary:ligand> .\n"
		for aru in res_uri:
			#add parts
			rm += "<"+hood_uri+"> <"+base_uri+"/lighood_vocabulary:has_member> <"+aru+"> .\n"
			#link ligand to neighbors
			rm += "<"+ligand_uri+"> <"+base_uri+"/lighood_vocabulary:has_neighbor> <"+aru+"> .\n"
			#type the neighbors
			rm += "<"+aru+"> <"+rdf+"type> <"+base_uri+"/lighood_vocabulary:neighbor> .\n"
	return rm

#creates an N3 file with aPdbId in the specified anOutputDirectory
# by parsing the ligand_uri_dict
def writeN3Hood(someNTriples, aPdbId, anOutputDirectory):
	if someNTriples:
		f = open(anOutputDirectory+'/'+aPdbId+'-ligand-neighborhood.nt','w')
		f.write(someNTriples)
		f.close()

#returns a defaultdict(list) where the key is the Bio2RDF URI of 
# the ligand residue in this structure and the value is a list of 
# residue URIs that are in the radius of the given ligand
def makeURIHood(aLigandDictList):
	rm = defaultdict(list)
	for (ligand, hood) in aLigandDictList.items():
		ligfi = ligand.get_full_id()
		#build ligand uri
		ligand_uri = base_uri+'/pdb_resource:'+ligfi[0]+'/chemicalComponent_'+ligfi[2]+str(ligfi[3][1])
		residue_uris = []
		for aResidue in hood:
			fi = aResidue.get_full_id()
			res_pdbid = fi[0]
			res_chain = fi[2]
			res_position = fi[3][1]
			res_uri = base_uri+'/pdb_resource:'+res_pdbid+'/chemicalComponent_'+res_chain+str(res_position)
			residue_uris.append(res_uri)
		if not ligand in rm:
			rm[ligand_uri] = residue_uris
	return rm	

# compute the ligand neighbourhood for the given pdb file. 
# someVerifiedLigands is the contents of the spreadsheet with known ligands
# aRadius is the threshold underwhich the comparison will be made
# this method returns a defaultdict(list) where the key is a ligand
# and the value is a list of PDB.Residue that exist within aRadius of the given
# ligand
def findNeighbours(aPdbFilePath, someVerifiedLigands, aRadius):
	rm = defaultdict(list)
	fn, fe = os.path.splitext(aPdbFilePath)
	pdbId = fn.rsplit('/')[-1]
	match = re.match('^\w{4}$', pdbId)
	if match:
		p = PDBParser(PERMISSIVE=1, QUIET=1)
		structure = p.get_structure(pdbId, aPdbFilePath)
		models = structure.get_list()
		#iterate over the models
		for aModel in models:
			chains = aModel.get_list()
			#get all the atoms ('A') in this model
			model_atoms = Selection.unfold_entities(aModel,'A')
			#create a neighbor search
			ns = NeighborSearch(model_atoms)
			#search the chains for any known ligands
			for aChain in chains:
				#get the residues in this chainln
				residues = aChain.get_list()
				#iterate over the residues
				for aR in residues:
					if findLigandIdInListOfLigands(someVerifiedLigands,aR.get_resname()) != None:
						#found a ligand!
						ligand = aR
						neighborset = computeNeighborSet(ns, ligand, aRadius)
						rm[ligand] = neighborset
	return rm

#given a BioPython.NeighborSearch object a ligand residue and a radius
# this method iterates over all the atoms in a ligand and finds the 
# residue neighbors every atom and stores it in a set
# this method then returns a set of all the residues that exist
# within the specified threshold distance (aRadius) of aLigand
def computeNeighborSet(aNeighborSearch, aLigand, aRadius):
	rm = set()
	#get the atoms ('A') of this residue
	ligand_atom_list = Selection.unfold_entities(aLigand, 'A')	
	#iterate over the list of atoms
	for anAtom in ligand_atom_list:
		#set a center
		center = anAtom.get_coord()
		neighbor_atoms = aNeighborSearch.search(center, aRadius)
		if neighbor_atoms:
			#get the residues that those atoms correspond to
			nr = Selection.unfold_entities(neighbor_atoms, 'R')
			#now add these residues to the rm set
			for ar in nr:
				#do not add the residue that may have the ligand residue name in it
				pos = ar.get_resname().find(aLigand.get_resname())
				if pos == -1:
					rm.add(ar)
	return rm

#downloads the latest version of the google doc
def fetchLigandList(aUrl):
	rm = []
	print "fetching latest ligand list ..."
	r = urllib2.urlopen(pdb_to_ligand_list_url)
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

#retrieve a list of files with .pdb as extension from the given local directory (alocaldir)
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

#start the program
if __name__ == '__main__':
  main(sys.argv)