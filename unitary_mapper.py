#!/home2/jjl83/python
# Maps variants in vcf file to a gencode interval annotation file 
# (unitary pseudogenes) using VAT genericMapper
# Creates fasta file of selected sequences, run through FASTX protein alignment.
# Generates fasta file of alternative sequences, run through FASTX.
# Two alignments compared for presence of stop and frameshifts.
# 07/11/2013 - Jeremy Liu (jeremy.liu@yale.edu)
# 02/07/2014 Fully incorperate genericMapper output to avoid manual computing
#	relative sequence position and correct parsing of transcript.

import os, sys, gzip
import variants, alignment

# Usage documentation if too few arguments
if len(sys.argv) < 7:
	print "Usage: ./unitary_variants.py <variant.vcf> <annotation_gtf> \
<annotation_interval> <annotation_fa> <database_fa> <output_directory>"
	sys.exit(1)

# Dependencies
from variants import *
from alignment import *

# Determine output directory and output file names	
filename = os.path.basename(sys.argv[1])
index = (len(filename) - 1) - filename[::-1].index('.')
basename = filename[:index]
outputPath = sys.argv[6]
if not os.path.exists(outputPath):
	os.mkdir(outputPath)
	
# Determine file descriptor for input vcf file (compressed .gz or uncompressed)
inputPath = sys.argv[1]
if inputPath.endswith(".gz"):
	inputPipe = Popen(['zcat', inputPath], stdout=PIPE)
	inputFile = inputPipe.stdout
else:
	inputFile = open(inputPath)

# Read in annotation gtf and fasta files INTERVAL FILE NOT NEEDED?
print 'Loading reference files...'
annotationGtf = parse_annotation_gtf(open(sys.argv[2], 'r'))
annotationFasta = parse_annotation_fa(open(sys.argv[4], 'r'))

# Set up genericMapper pipeline to process input vcf
print 'Running VAT genericMapper on input vcf file...'
vatPipe = Popen(['/home2/jjl83/aloft/vat-bin/Linux_x86_64/genericMapper', sys.argv[3], 'unitary_pseudogene'], 
				stdin=inputFile, stdout=PIPE)
vatOutput = vatPipe.stdout

# Process genericMapper output file, also in vcf format
print 'Reading VAT genericMapper output vcf file...'
intersects = {}
intersectTranscripts = []
for var in vatOutput:
	if not var.startswith('#'):
		lineComponents = var.strip().split('\t')
		info = lineComponents[7].split(':')
		transcripts = [item for item in info if item.startswith('ENST')]
		for transcript in transcripts:
		#transcript = lineComponents[7].split(':')[-2]
			uniqTranscript = transcript + '_' + lineComponents[1]
			intersects[uniqTranscript] = [lineComponents, transcript]
			if transcript not in intersectTranscripts:
				intersectTranscripts.append(transcript)

print 'Calculating alternate sequences...'
newSeqFile = open(os.path.join(outputPath, basename + '.alternate.fa'), 'w')
for uniqTranscript, intersect in intersects.iteritems():
	var = intersect[0]
	#print var
	info = var[-1].split(':')
	transcript = intersect[1]
	transcriptIndex = info.index(transcript)
	transcriptLength = int(info[transcriptIndex + 1].split('_')[0])
	transcriptGenes = [x for x in info if x.startswith("ENSG")]
	transcriptStrand = info[info.index(transcriptGenes[0]) + 1]
	# Make sure to handle ###_###|# case (ENST00000379669.4)
	noStrandRelPos = int(info[transcriptIndex + 1].split('_')[-1].split('|')[0])
	if transcriptStrand == '+':
		relativePosition = noStrandRelPos
	elif transcriptStrand == '-':
		relativePosition = transcriptLength - noStrandRelPos + 1
	else:
		sys.stderr.write("BAD STRAND:" + transcriptStrand + '\t' + str(var) + '\n')
	#print '\t' + str(relativePosition)
	newSequence = sequence_substitute(transcript, annotationFasta[transcript][1], 
									var[1], var[3], var[4], relativePosition)
	#newtranscript = transcript + '_' + var[1]
	header = annotationFasta[transcript][0]
	dashes = list_duplicates_of(header, '|')
	newSeqFile.write(header[:int(dashes[0])] + '|' + uniqTranscript + header[int(dashes[1]):] + '\n')
	newSeqFile.write(newSequence + '\n')
newSeqFile.close()

print 'Running FASTX on original sequences...'
fasta2file(annotationFasta, intersectTranscripts, os.path.join(outputPath, basename + '.original.fa'))
run_fastx_alignment(os.path.join(outputPath, basename + '.original.fa'), sys.argv[5], 
	                os.path.join(outputPath, basename + '.original.fastx'))
	
print 'Running FASTX on alternate sequences...'
run_fastx_alignment(os.path.join(outputPath, basename + '.alternate.fa'), sys.argv[5], 
	                os.path.join(outputPath, basename + '.alternate.fastx'))

print 'Reading original FASTX alignments...'
fastxFile = open(os.path.join(outputPath, basename + '.original.fastx'), 'r')
alignOriginal = parse_alignment(fastxFile)
fastxFile.close()

print 'Reading alternate FASTX alignments...'
fastxFileAlt = open(os.path.join(outputPath, basename + '.alternate.fastx'), 'r')
alignAlt = parse_alignment(fastxFileAlt)
fastxFileAlt.close()

print 'Comparing alignments...'
outFile = open(os.path.join(outputPath, basename + '.unitary.vcf'), 'w')
outFile.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
for uniqTranscript in alignAlt.keys():
	transcript = uniqTranscript[:uniqTranscript.index('_')]
	commonOutput = (';UT=' + uniqTranscript + ':AQD=' + compute_query_density(alignAlt[uniqTranscript]) + 
		            ':ASD=' + compute_subject_density(alignAlt[uniqTranscript]))
	if not transcript in alignOriginal.keys():
		outFile.write('\t'.join(intersects[uniqTranscript][0]) + commonOutput + ':CHANGE=No_Original_Hit\n')
	if find_stop_fs(alignOriginal[transcript]) == True and find_stop_fs(alignAlt[uniqTranscript]) == False:
		outFile.write('\t'.join(intersects[uniqTranscript][0]) + commonOutput + ':CHANGE=Functional\n')
	elif find_stop_fs(alignOriginal[transcript]) == True and find_stop_fs(alignAlt[uniqTranscript]) == True:
		outFile.write('\t'.join(intersects[uniqTranscript][0]) + commonOutput + ':CHANGE=No_Effect\n')
	else:
		outFile.write('\t'.join(intersects[uniqTranscript][0]) + commonOutput + ':CHANGE=No_Original_StopFS\n')

outFile.close()
inputFile.close()
print 'Done.'
