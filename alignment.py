#!/usr/bin/env python
# <Description>
# 7/11/2013 - Jeremy Liu

import os, sys, variants
from subprocess import Popen, PIPE
from variants import list_duplicates_of

# Function that returns true if a stop or frame shift is found, false otherwise
def find_stop_fs(fastxPair):
	item = fastxPair[0]
	original = fastxPair[1]
	for query in item[5:]:
		sequence = query[0][2]
		if '*' in sequence:
			return True
		if '\\' in sequence:
			return True
		if '/' in sequence:
			return True
	return False
	
def compute_query_density(fastxPair):
	item = fastxPair[0]
	queryLength = int(item[2])
	queryStart = int(item[5][0][1])
	queryEnd = int(item[-1][0][3])
	return str(abs(queryEnd - queryStart + 1) / float(queryLength))

def compute_subject_density(fastxPair):
	item = fastxPair[0]
	subjectLength = int(item[4])
	subjectStart = int(item[5][1][1])
	subjectEnd = int(item[-1][1][3])
	return str(abs(subjectEnd - subjectStart + 1) / float(subjectLength))

# Function to read in all 'Query' and 'Subject' lines within a queried alignment.
def getMatchSequences(infile, list, original):
	found = 0
	while True:
		line = infile.readline()
		if not line:
			break
		l = line.strip().split()
		if len(l) > 0:
			if l[0] == '!!':
				return False
			elif '>>>' in l[0] and found == 1:
				infile.seek(-len(line),1)
				return True
			elif l[0][0] == '>' and found == 1:
				return True
			elif l[0] == 'Lambda':
				return True
			elif l[0][0] == '>':
				found = 1
				list.append(l[0][1:])
				original.append(line.strip())
				while True:
					line = infile.readline()
					if not line:
						break
					original.append(line.strip())
					if len(line.strip().split()) > 0 and line.startswith('Length='):
						list.append(line.strip()[line.index('=')+1:])
					elif len(line.strip().split()) > 0 and line.strip().split()[0] == 'Frame':
						break
			elif l[0] == 'Score' and len(l) > 2:
				return True
			elif l[0] == 'Query':
				original.append(line.strip('\n'))
				match = infile.readline()
				original.append(match.strip('\n'))
				subject = infile.readline()
				list.append([l, subject.strip().split()])
				original.append(subject.strip('\n'))
	return True

# Read in FASTX file
# List element: [calculation items, original lines]
# Item: [GeneID, TranscriptID, QueryLen, ProteinID, SubjectLen [[Query, Pos, Sequence, Pos][Sbjct, Pos, Sequence, Pos]]...]
def parse_alignment(fastxFileFd):
	fastx = {}
	while True:
		line = fastxFileFd.readline()
		if not line:
			break
		original = []
		item = []
		l = line.strip().split()
		if len(l) > 0:
			if '>>>' in l[0]:
				query = line.strip()
				gene = query.split('|')[0]
				item.append(gene[(gene.index('>>>') + 3):])
				item.append(query.split('|')[1])
				while True:
					line = fastxFileFd.readline()
					if not line:
						break
					if 'Library:' in line:
						break
					query = query + line.strip()
				item.append(l[-2])
				original.append(query)
				if getMatchSequences(fastxFileFd, item, original) == True:
					fastx[item[1]] = [item, original] 
	return fastx

# Function to write a dictionary of fasta sequences to file, only those specified by TranscriptID
def fasta2file(fastaDict, transcripts, outPath):
	tempFile = open(outPath, 'w')
	for transcript in transcripts:
		if transcript in fastaDict.keys():
			tempFile.write(fastaDict[transcript][0] + '\n')
			tempFile.write(fastaDict[transcript][1] + '\n')
	tempFile.close()

# Function to run fastx alignment given a fasta file path, db path, and an outpath.
def run_fastx_alignment(seqsFastaPath, dbFastaPath, outpath):
	with open(outpath, 'w') as output:
		fastxPipe = Popen(['/home2/jjl83/fasta-36.3.5d/bin/fastx36', '-m', 'B', seqsFastaPath, dbFastaPath], stdout=output)
		fastxPipe.communicate()
	