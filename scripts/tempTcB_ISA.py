#!/usr/bin/env python
import csv
import random
import sys
import time
from operator import attrgetter

import Bio.Entrez
import Bio.Seq
import Bio.SeqIO
import Bio.SeqRecord
import colorama
import pysam
import runbin
import mappedreads
import twobitreader
from Bio.Alphabet import IUPAC

#import time
#import StringIO
#import ConfigParser
#import httplib

#Various functions:
mapreads=0
read_sam_file=0 #read input from samfile
read_csv=0 #read input from csv file
write_csv=0 # write csv (only relevant if reading from sam file)
getseqs=0#1#get sequences from either entrez or local TwoBit genome around locations indicated by genome_location_csv
sequencesource='LOCAL' #get sequences from 'LOCAL' TwoBit genome or 'REMOTE' entrez server
writeFASTA=0#1#write a fasta file of all sequences around IS
compressreads=0 #remove duplicate reads and reads shifted +/- 1 nt, number of reads are compressed in csv and fasta
writelogo=1#1#create a logo image of consensus sequence # requires getseqs and writeFASTA to be on
writevepfile=0#write a vep file to use with vep, either online or by setting getannotations to 1
getannotations=0#1 # uses VEP to get anotations for insertion sites
userandomIS=0#0 #useful for generating a random control data set chromosome matched to query set.

lwindow=70 #window on either side of indicated nt location to return
rwindow=30
#specified file names
rootname='merged'#'13-wind20-70'
samwindow=60
sam_file=f'samfiles/{rootname}.sam'#,'02-wind20-70.sam','03-wind20-70.sam','04-wind20-70.sam','05-wind20-70.sam','06-wind20-70.sam','07-wind20-70.sam','08-wind20-70.sam','09-wind20-70.sam','10-wind20-70.sam','11-wind20-70.sam','12-wind20-70.sam','13-wind20-70.sam',]
genome_location_csv=f'{rootname}2bit_IS_data.csv'#'IS_data.csv'
FASTAfile='merged_real.fasta'#f'{rootname}_retrieved_2bit.fasta'
outputVEPfile=f"{rootname}_VEPfile.csv"
logofile=f"{rootname}IS_logo.svg"
annotationsfile=f'{rootname}_IS_annotations.csv'

chromosome_ids = 'chromosomes.csv'
veplocation = '../ensembl-vep/vep'
bowtielocation = 'bowtie2'
weblogolocation = 'weblogo' #in PATH
twobitlocation = '../scripts/TwoBitToFa'
twobitgenomelocation='../genomes/GRCh38.2bit'
Bio.Entrez.email = "bryan@bmogen.com"     # Always tell NCBI who you are
Bio.Entrez.api_key = "679368bc3fef47bb440fb743c889befe4e09" #bryan's personal NCBI key (allows faster requests of their server)

#sanity checks:
if read_sam_file == read_csv == 1:
	sys.exit(colorama.Fore.RED + 'Specify only a sam or csv file not both')
if getseqs or writeFASTA or writevepfile:
	if read_sam_file == read_csv == 0:
		sys.exit(colorama.Fore.RED + 'Cannot get sequences, write a FASTA, or write a VEP file without being given a csv or sam file')
if writeFASTA == 1 and getseqs == 0:
	sys.exit(colorama.Fore.RED + "Can't write a fasta file without getting sequences.")
if sequencesource != 'LOCAL' and sequencesource != 'REMOTE' and getseqs == 1:
	sys.exit(colorama.Fore.RED + "Must specify 'LOCAL' or 'REMOTE' source to get sequences")
'''
if read_sam_file 
mapreads=0
read_sam_file=1 #read input from samfile
read_csv=0 #read input from csv file
write_csv=1 # write csv (only relevant if reading from sam file)
getseqs=0#1#download sequences from entrez around locations indicated by genome_location_csv
writeFASTA=0#1#write a fasta file of all sequences around IS
writelogo=0#1#create a logo image of consensus sequence # requires getseqs and writeFASTA to be on
writevepfile=0#write a vep file to use with vep, either online or by setting getannotations to 1
getannotations=0#1 # uses VEP to get anotations for insertion sites
userandomIS=0#0 #useful for generating a random control data set chromosome matched to query set.
'''

reads=[]
recordlist=[]        
chromIDS={}
chromNTS={}
programstart = time.time()
#make a dictionary of chromosome names and ID #s usich chromosomes.csv, which is the table copied from https://www.ncbi.nlm.nih.gov/genome?term=human&cmd=DetailsSearch
with open(chromosome_ids) as csv_file:
	csv_reader = csv.reader(csv_file, delimiter=",")
	line_count=0
	for row in csv_reader:
		name=row[1]
		ID=row[2]
		nts=row[4]
		print(f'chromosome ' + colorama.Fore.YELLOW + f'{name}'+colorama.Style.RESET_ALL + ' is called ' + colorama.Fore.YELLOW + f'{ID}'+colorama.Style.RESET_ALL + ' and has '+ colorama.Fore.YELLOW + f'{nts}'+colorama.Style.RESET_ALL + ' nucleotides')
		chromIDS[name]=ID
		chromNTS[name]=nts
if userandomIS:
	random.seed(a=None, version=2)
#parse a row from the read file and pull the sequence
class read_csv_line(object):
	def __init__(self, row):
		#tries = 1
		self.chrom=str(row[0])
		self.sense=row[1]
		if row[1]=="+":
			self.sensenum=1
		elif row[1]=="-":
			self.sensenum=2
		if userandomIS: #if random IS used, per "read" use a random number for insertion site
			self.loc=random.randrange(int(chromNTS[str(self.chrom)])) # random int returned in the range of the given chromosome range
		else: #if not random, use given values
			self.loc=int(row[2])
		self.gene=row[3]
		self.totalseqs=row[4]
		self.total=row[5]
		self.plasmid=row[6]
		self.sample1=row[7]
		self.sample2=row[8]
		self.sample3=row[9]
		self.sample4=row[10]
		self.record=None
		self.sequence='-'
class get_seq(object):
	def __init__(self, read_item):
		#tries = 1
		self.chrom=str(read_item.chrom)
		self.sense=read_item.sense
		self.sensenum=read_item.sensenum
		self.loc=read_item.loc
		self.gene=read_item.gene
		self.totalseqs=read_item.totalseqs
		self.total=read_item.total
		self.plasmid=read_item.plasmid
		self.sample1=read_item.sample1
		self.sample2=read_item.sample2
		self.sample3=read_item.sample3
		self.sample4=read_item.sample4
		self.record=read_item.record
		self.sequence=read_item.sequence
		#self=read_item
		if getseqs: #if get sequences, fetch sequence centered on given location, with a window on either side
			if sequencesource =='LOCAL':
				self.record=fetch_seq_from_twobit(self.chrom, self.sense, self.loc)
				#if self.record.seq[0:9]=="NNNNNNNNNN":
				#	print("WARNING, got invalid sequence, trying again")
				#	self.record=fetch_seq_from_twobit(self.chrom, self.sense, self.loc)
				#if self.record.seq[0:9]=="NNNNNNNNNN":
				#	print("WARNING, second attepmt also got invalid sequence, proceding with 'NNNN...' sequence.")
				self.sequence=self.record.seq
			elif sequencesource == 'REMOTE':
				handle = Bio.Entrez.efetch(db="nucleotide", #fetches the specified sequence as a fasta
					id=chromIDS[self.chrom], #id for indicated chromosome from refrence genome assembly
					rettype="fasta", #format
					strand=self.sensenum, #which strand
					#the behaviour of seq_start and seq_stop might vary depending on sense vs antisense
					seq_start=str(self.loc-lwindow), #range of sequence to return
					seq_stop=str(self.loc+rwindow-1))
				self.record = Bio.SeqIO.read(handle, "fasta") #reads the returned fasta into a sequence object
				handle.close()#close efetch 
				self.sequence=self.record.seq #adds the sequence property as a string 
#fetch indicated sequence from local 2bit genome. chrom given as number (or X/Y/NT), sense given as '+'/'-'.
def fetch_seq_from_twobit(chrom, sense, loc):
	if sense=='-':
		seq_start=int(loc+samwindow+lwindow) #range of sequence to return
		seq_stop=int(loc+samwindow-rwindow)
		sequence=genome[f'chr{chrom}'][seq_stop:seq_start]
		#print(f'reading from chromosome {chrom}, {sense} strand, from {seq_start} to {seq_stop}. \n sequence is {sequence}')
		record=Bio.SeqRecord.SeqRecord(Bio.Seq.reverse_complement(Bio.Seq.Seq(sequence, IUPAC.unambiguous_dna)), id=f'chr{chrom}{sense}:{seq_start}-{seq_stop}', description='')
		#return Bio.SeqRecord.SeqRecord(Bio.Seq.reverse_complement(Bio.Seq.Seq(genome[f'chr{sequencingread.chrom}'][seq_stop:seq_start], IUPAC.unambiguous_dna)), id=f'chr{sequencingread.chrom}{sequencingread.sense}:{seq_start}-{seq_stop}', description='')
		return record
	elif sense=='+':
		seq_start=int(loc-lwindow) #range of sequence to return
		seq_stop=int(loc+rwindow)
		sequence=genome[f'chr{chrom}'][seq_start:seq_stop]
		record=Bio.SeqRecord.SeqRecord(Bio.Seq.Seq(sequence, IUPAC.unambiguous_dna), id=f'chr{chrom}{sense}:{seq_start}-{seq_stop}', description='')
		#return Bio.SeqRecord.SeqRecord(Bio.Seq.Seq(genome[f'chr{sequencingread.chrom}'][seq_start:seq_stop], IUPAC.unambiguous_dna), id=f'chr{sequencingread.chrom}{sequencingread.sense}:{seq_start}-{seq_stop}', description='')
		return record
	
'''
combine python & bash script
if mapreads:
	bowtiecommand=f'{veplocation} -i {outputVEPfile} -o {annotationsfile} --buffer_size 100000 --merged --format "ensembl" --cache --nearest gene --distance 1 --offline --use_given_ref --fork 4 --pick' #runs vep using "--merged" ensembl and refseq genomes (what i currently have "--cached" on my machine). return the "--nearest gene". --distance from query to look for genes. maybe add --symbol --protein --ccds. Maybe add back: --numbers --domains --biotype --hgvs add --merged --gencode_basci
	print('running VEP analysis to annotate genomic positions using ensembl vep. Writting output to {annotationsfile} and {annotationsfile}_summary.html')
	print(bowtiecommand)
	bowtie=runbin.Command(vepcommand)
	vepout=vep.run(timeout=20000)
	Bash script:
	fa=".fasta"
	sm=".sam"
	lg=".log"
	for VARIABLE in  01-wind20-70 02-wind20-70 03-wind20-70 04-wind20-70 05-wind20-70 06-wind20-70   07-wind20-70   08-wind20-70   09-wind20-70   10-wind20-70   11-wind20-70   12-wind20-70   13-wind20-70
	do
		echo 	bowtie2 -f -x ./GRCh38.fna.bowtie_index/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index -U $VARIABLE$fa -S $VARIABLE$sm
		bowtie2 -f -x ./GRCh38.fna.bowtie_index/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index -U $VARIABLE$fa -S $VARIABLE$sm 2>&1 | tee $VARIABLE$lg
	done

'''
if sequencesource =='LOCAL':
#	GRCh38=Bio.SeqIO.index('../genomes/GRCh38.fasta',"fasta")
	genome=twobitreader.TwoBitFile(twobitgenomelocation) #reads genome into object
if read_sam_file:
	recordlist, unmapped = mappedreads.read_sam(sam_file)
	'''
	samfile = pysam.AlignmentFile(sam_file, "r")
	
	print(colorama.Fore.GREEN + f'Reading sam file: '+colorama.Style.RESET_ALL+f' {sam_file}')
	readslist=[]
	unmapped=[]
	for entry in samfile:
		if entry.reference_name and (entry.reference_name[3:] in chromIDS.keys()): #check if read has a refrence to a chromosome, otherwise ignore (i.e. unmapped reads)
			chrom=entry.reference_name[3:] #chromosome number of mapping
			address=entry.get_reference_positions()[0] #chromosome position of mapping
			if entry.is_reverse:
				sense="-"
			else:
				sense="+"
			readslist.append(read_csv_line([chrom, sense, address,'',1,'','','','','','']))
		else:
			unmapped.append(entry)
	readslist=sorted(readslist, key=attrgetter('chrom', 'loc')) #sort the list by chromosome# and location
	semicompressed_readslist=[]
	compressed_readslist=[]
	if compressreads:
		for i in range(len(readslist)):
			if i+1<len(readslist) and readslist[i].loc==readslist[i+1].loc and readslist[i].chrom == readslist[i+1].chrom: #if this read isn't the last one and if it's the same as the next one
				readslist[i+1].totalseqs=readslist[i+1].totalseqs+readslist[i].totalseqs #if they are the same chromosome and locus, add to the total
			else:
				semicompressed_readslist.append(readslist[i]) #if it's different than the next one, add it to the list (with the total count)
		#adjust loc +/- 1 to account for indels
		for i in range(len(semicompressed_readslist)):
			if i+1<len(semicompressed_readslist) and semicompressed_readslist[i].loc==(semicompressed_readslist[i+1].loc-1) and semicompressed_readslist[i].chrom == semicompressed_readslist[i+1].chrom: #if this read isn't the last one and if it's the same as the next one
				print(colorama.Style.RESET_ALL + f'merging adjacent read locations. read '+colorama.Fore.YELLOW + f'{i}'+colorama.Style.RESET_ALL + ', locations '+colorama.Fore.YELLOW + f'{semicompressed_readslist[i].loc} '+colorama.Style.RESET_ALL + ' and '+colorama.Fore.YELLOW + f'{semicompressed_readslist[i+1].loc}.')
				if semicompressed_readslist[i].totalseqs>semicompressed_readslist[i+1].totalseqs:
					semicompressed_readslist[i+1].loc=semicompressed_readslist[i].loc #if read i has more counts than i+1, then change i+1's location to i
				semicompressed_readslist[i+1].totalseqs=semicompressed_readslist[i+1].totalseqs+semicompressed_readslist[i].totalseqs #if they are the same chromosome and locus, add to the 
			else:
				print(colorama.Style.RESET_ALL + f'adding line '+colorama.Fore.YELLOW + f'{i}'+colorama.Style.RESET_ALL + ' to compressed readslist for addition to csv')
				compressed_readslist.append(semicompressed_readslist[i]) #if it's different than the next one, add it to the list (with the total count)
		readslist=compressed_readslist
	line_count=0
	for line in readslist:
		thisrow=get_seq(line)
		reads.append(thisrow)
		recordlist.append(thisrow.record)
		print(colorama.Style.RESET_ALL + f'Sequence number '+colorama.Fore.YELLOW + f'{line_count} '+colorama.Style.RESET_ALL + 'read on chromosome '+colorama.Fore.GREEN + f'{thisrow.chrom}'+colorama.Style.RESET_ALL + ' with '+colorama.Fore.YELLOW + f'{thisrow.sense}'+colorama.Style.RESET_ALL + ' sense at location '+colorama.Fore.YELLOW + f'{thisrow.loc} '+colorama.Style.RESET_ALL + 'with sequence '+colorama.Fore.BLUE + f'{thisrow.sequence}',end="\r" )
		line_count += 1 
	print(colorama.Style.RESET_ALL + f'\nProcessed '+colorama.Fore.GREEN + f'{line_count} lines.')
	'''

if write_csv:
	mappedreads.write_csv(recordlist,genome_location_csv)
	'''
	print(colorama.Style.RESET_ALL + f'Writing reads to csv file: '+colorama.Fore.YELLOW + f'{genome_location_csv}')
	with open(genome_location_csv, "w", newline="") as csv_file:
		csv_writer = csv.writer(csv_file, delimiter=',') #comma delimited csv file
		#write csv headder to match expected format
		csv_writer.writerow(['Chrom','Sense','Loc','Total # IS Sequences Found','Total # IS Found','Plasmid m995','Sample1','Sample2','Sample3','Sample4'])
		csv_writer.writerow([])
		csv_writer.writerow([])
		for i in readslist:
			csvline=[i.chrom,i.sense,i.loc,i.gene,i.totalseqs,i.total,i.plasmid,i.sample1,i.sample2,i.sample3,i.sample4]
			csv_writer.writerow(csvline)
	'''
#Chrom,Sense,Loc,Gene,Total # IS Sequences Found,Total # IS Found,Plasmid m995,Sample1,Sample2,Sample3,Sample4
if read_csv:
	readslist=mappedreads.read_csv(genome_location_csv)
	'''
	print(colorama.Style.RESET_ALL + f'reading data from '+colorama.Style.GREEN + f'{genome_location_csv}')
	with open(genome_location_csv) as csv_file: #open file
		csv_reader = csv.reader(csv_file, delimiter=',') #read lines
		line_count = 1
		for row in csv_reader:
			if line_count < 4:#first 4 lines are headders
				if line_count==1: print(colorama.Fore.YELLOW + f'Column names are {", ".join(row)}')
				line_count += 1
			else: #read data from each line
				thisrow=read_csv_line(row)
				thisrow=get_seq(thisrow)
				reads.append(thisrow)
				recordlist.append(thisrow.record)
				print(colorama.Style.RESET_ALL + f'Sequence number {line_count-3} '+colorama.Style.RESET_ALL + 'read on chromosome '+colorama.Style.GREEN + f'{thisrow.chrom}'+colorama.Style.RESET_ALL + ' with '+colorama.Style.ORANGE + f'{thisrow.sense}'+colorama.Style.RESET_ALL + ' sense at location '+colorama.Style.YELLOW + f'{thisrow.loc} '+colorama.Style.RESET_ALL + 'with sequence '+colorama.Style.BLUE + f'{thisrow.sequence}',end="\r" )
				line_count += 1 
		print(f'Processed {line_count} lines.')
	'''
if writeFASTA:#write FASTA formatted file with returned sequences
	print(colorama.Style.RESET_ALL + f'Writting fasta file:'+ colorama.Fore.GREEN + f' {FASTAfile}')
	Bio.SeqIO.write(recordlist, FASTAfile, "fasta")
if writelogo: #dependant on having sequences, optional to make logo plot
	weblogocommand=f'{weblogolocation} -f {FASTAfile} -D fasta -o {logofile} -F svg -A dna -F png --resolution 600 -s large -c classic -i {str(int(20-1*lwindow))} -l -10 -u 10 '#-l [lower bound] -u [upper bound]
	print(colorama.Style.RESET_ALL + f'Writting logo')
	print(colorama.Fore.YELLOW + f"'{weblogocommand}'")
	logocommand=runbin.Command(weblogocommand)
	logoout=logocommand.run(timeout=1800) 

if writevepfile:
	reads = sorted(reads, key = lambda x: (x.chrom, x.loc))
	print(colorama.Style.RESET_ALL + f'Writting file for VEP analysis: '+colorama.Fore.YELLOW + f'{outputVEPfile}')
	with open(outputVEPfile, 'w', newline='') as formatted_csv_file: #writting file formatted for VEP analysis
		csv_writer = csv.writer(formatted_csv_file, delimiter='\t') #tab delimited csv file
		for entry in reads:
			if entry.chrom == 23:
				entry.chrom ="X"
			elif entry.chrom == 24:
				entry.chrom ="Y"
			csv_writer.writerow([entry.chrom, str(entry.loc+1), entry.loc, "-/A", entry.sense])

if getannotations:
	vepcommand=f'{veplocation} -i {outputVEPfile} -o {annotationsfile} --force_overwrite --buffer_size 100000 --merged --format "ensembl" --cache --nearest gene --distance 1 --offline --use_given_ref --fork 4 --pick' #runs vep using "--merged" ensembl and refseq genomes (what i currently have "--cached" on my machine). return the "--nearest gene". --distance from query to look for genes. maybe add --symbol --protein --ccds. Maybe add back: --numbers --domains --biotype --hgvs add --merged --gencode_basci
	print(colorama.Style.RESET_ALL + f'Running VEP analysis to annotate genomic positions using ensembl vep.'+colorama.Fore.YELLOW + f'"{vepcommand}"')
	vep=runbin.Command(vepcommand)
	vepout=vep.run(timeout=20000)
	if len(vepout[1])>0:
		print(vepout[1]) #STDOUT
	if len(vepout[2])>0:
		print(vepout[2]) #ERRORS

programend = time.time()
message=(colorama.Fore.GREEN + f'Completed all tasks in {int(programend-programstart)} seconds. Exiting.'+colorama.Style.RESET_ALL)
print(message)
exit()
