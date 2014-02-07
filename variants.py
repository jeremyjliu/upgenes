#!/usr/bin/env python
# <Description>
# 7/11/2013 - Jeremy Liu

import os, sys

def list_duplicates_of(seq,item):
    start_at = -1
    locs = []
    while True:
        try:
            loc = seq.index(item,start_at+1)
        except ValueError:
            break
        else:
            locs.append(loc)
            start_at = loc
    return locs

# UNUSED Function to read in variants file in VCF format
# Format: [[CHR#, POS, ID, REF, ALT, QUAL, FILTER, INFO]...]
def parse_variants(inVcfFd):
	records = []
	for line in inVcfFd:
		lineComponents = line.strip().split('\t')
		if line.startswith('#CHROM') or not line.startswith('#'):
			records.append(lineComponents)
	return records

# Function to read in gencode annotation gtf file
# Format: {TranscriptID: [[CHR#, ExonStart, ExonEnd, Strand]...]}
def parse_annotation_gtf(annotationGtf):
	gencodeAnnotation = {}
	for line in annotationGtf:
		lineComponents = line.strip().split('\t')
		key = line.strip().split()[11].rstrip(';').strip('"')
		if key in gencodeAnnotation.keys():
			gencodeAnnotation[key].append([lineComponents[0][3:], lineComponents[3], lineComponents[4], lineComponents[6]])
		else:
			gencodeAnnotation[key] = [[lineComponents[0][3:], lineComponents[3], lineComponents[4], lineComponents[6]]]
	return gencodeAnnotation
	
# Function to read in gencode annotation fasta file
# Format: {TranscriptID: [Header, Sequence]}
def parse_annotation_fa(annotationFasta):
	gencodeAnnotation = {}
	for line in annotationFasta:
		sequence = annotationFasta.next()
		key = line.strip().split('|')[1]
		if key in gencodeAnnotation.keys():
			sys.stderr.write('WARNING: Multiple sequences per unique TranscriptID\n')
		gencodeAnnotation[key] = [line.strip(), sequence.strip()]
	return gencodeAnnotation

# UNUSED Function to determine whether a variants intersects exons in gtf file
# Returns TranscriptID if true, returns False otherwise
def compute_intersect(variantLineComponents, gencodeGtf):
	pos = int(variantLineComponents[1])
	for transcript in gencodeGtf.keys():
		for exon in gencodeGtf[transcript]:
			if exon[0] == variantLineComponents[0] and int(exon[1]) <= int(pos) and int(pos) <= int(exon[2]):
				return transcript
	return False

# Function to translate absolute genomic position into relative 1-indexed based sequence position.
def find_sequence_position(transcript, sequence, genePos, gtf):
	pos = 1
	if transcript in gtf.keys():
		print gtf[transcript]
		for exon in gtf[transcript]:
			start = exon[1]
			end = exon[2]
			if int(end) < int(genePos):
				print 'outside' + str(exon)
				pos += int(end) - int(start) + 1
			elif int(start) <= int(genePos) and int(genePos) <= int(end):
				print 'inside' + str(exon)
				pos += int(genePos) - int(start)			
	else:
		sys.stderr.write('WARNING: Transcript not found in gtf file in sequence substitution\n')
	return pos

# Function to substitute an alternate allele for a reference allele.
def sequence_substitute(transcript, sequence, genePos, ref, alt, relPos):
	#relPos = find_sequence_position(transcript, sequence, genePos, gtf)
	#print(transcript + '\t' + sequence + '\t' + str(genePos) + '\t' + str(relPos) + '\t' + ref)
	print '\t' + sequence + ' : ' + str(relPos)
	if not ref == sequence[relPos - 1]:
		print transcript + '\t' + str(genePos) + '\t' + str(relPos) + '\t' + ref + '\t' + sequence[relPos - 1] + '\t' + alt
	if len(ref) == 1 and len(alt) == 1:
		newSequence = sequence[:relPos-1] + alt + sequence[relPos:]
	else:
		newSequence = sequence[:relPos-1] + alt + sequence[relPos-1+len(ref):]
	return newSequence
	
# Function to substitute an alternate allele for a reference allele, with African style formatting.
def sequence_substitute_african(transcript, sequence, genePos, ref, alt, gtf):
	relPos = find_sequence_position(transcript, sequence, genePos, gtf)
	if not ref == sequence[relPos - 1]:
		print transcript + '\t' + str(genePos) + '\t' + str(relPos) + '\t' + ref + '\t' + sequence[relPos - 1] + '\t' + alt
	if len(ref) == 1 and len(alt) == 1:
		newSequence = sequence[:relPos-1] + alt + sequence[relPos:]
	else:
		newSequence = sequence[:relPos-len(ref)] + alt + sequence[relPos:]
	return newSequence
