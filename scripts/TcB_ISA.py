#!/btapps/miniconda3/bin/python3
import _mypath
import csv
import random
import sys
import time
import FASTQ
import Bio.SeqIO
import colorama
import os
#import pysam
import runbin
import mappedreads
import getseq
import twobitreader
import annotate

#inputs:
inputtype='fastq' #sam, csv, fasta, fastq
inputfile='05.fastq' #input file name
userandomIS=0#0 #useful for generating a random control data set chromosome matched to query set.
userandomSEQS=0 #generate random sequence data, usefull for calculating false positive rate when mapping reads w/ bowtie

#tasks:
getseqs=1#1#get sequences from either entrez or local TwoBit genome around locations indicated by genome_location_csv
compressreads=1 #remove duplicate reads and reads shifted +/- 1 nt, number of reads are compressed in csv and fasta
getannotations=0#1 # uses VEP to get anotations for insertion sites (needs writevepfile)

#parameters
barcode='ATCTGCGACG' #5' barcode sequence 05:ATCTGCGACG
lwindow=50 #window on either side of indicated nt location to return
rwindow=50
samwindow=0#60
sequencesource='LOCAL' #get sequences from 'LOCAL' TwoBit genome or 'REMOTE' entrez server
minimum_read_len=25 #minimum length of a read to try mapping. 25 will usually avoid any false positives in read sets of 200k reads
primer5='GGGTTCCGCCGGATGGC' #5' primer sequenc to remove from reads 
primer3='CCTAACTGCTGTGCCACT' #3' primer sequence to remove from reads
trim5=21 #additional (non-genomic) nts to trim off of 3' end of reads
trim3=16 #additional (non-genomic) nts to trim off of 5' end of reads (starts with TA)

#outputs:
write_csv=1 # write csv (only relevant if reading from sam file)
writeFASTA=1#1#write a fasta file of all sequences around IS
writelogo=1#1#create a logo image of consensus sequence # requires getseqs and writeFASTA to be on
writevepfile=0#write a vep file to use with vep, either online or by setting getannotations to 1

#specified file names
rootname=os.path.splitext(inputfile)[0]

#output file names
sam_file=f'{rootname}.sam'#,'02-wind20-70.sam','03-wind20-70.sam','04-wind20-70.sam','05-wind20-70.sam','06-wind20-70.sam','07-wind20-70.sam','08-wind20-70.sam','09-wind20-70.sam','10-wind20-70.sam','11-wind20-70.sam','12-wind20-70.sam','13-wind20-70.sam',]
genome_location_csv=f'{rootname}_IS_mappings.csv'#'IS_data.csv'
trimmedfastq=f'{rootname}_trimmed.fastq' # fastq file with adapters/primers/barcodes trimmed off with cutadapt
FASTAfile=f'{rootname}_retrieved_2bit.fasta'#output fasta file containing sequences surrounding mapped insertion site
outputVEPfile=f"{rootname}_VEPfile.csv" #VEP file containing mapped locations of insertions in VEP format so VEP can provide annotations
logofile=f"{rootname}IS_logo.svg" # logo file showing consensus integration site in logo format
annotationsfile=f'{rootname}_IS_annotations.csv'
logfile=f'{rootname}.log'
ISsamfilename=f'{rootname}IS.bam'#sam file name for single nt IS mappings

#refrence file locations
chromosome_ids = './refrence_datasets/chromosomes.csv'
veplocation = './ensembl-vep/vep'
bowtielocation = 'bowtie2'
bowtieindex = './refrence_datasets/genomes/GRCh38.fna.bowtie_index/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index'
weblogolocation = 'weblogo' #in PATH
twobitlocation = './scripts/TwoBitToFa'
twobitgenomelocation='./refrence_datasets/genomes/GRCh38.2bit'

#sanity checks:
if inputtype=="sam":
	read_sam_file=1
	mapreads=0
elif inputtype=="csv":
	read_csv=1
	mapreads=0
	if getseqs==0:
		sys.exit(colorama.Fore.RED + 'Given a csv file, but not asked to get seqs. Nothing to do.')
elif inputtype=="fasta":
	read_fasta=1
	mapreads=1
elif inputtype=="fastq":
	read_fastq=1
	mapreads=1
else:
	sys.exit(colorama.Fore.RED + 'Specify input file type as \'sam\', \'csv\', \'fasta\', or \'fastq\'')
#if getseqs or writeFASTA or writevepfile:
#	if read_sam_file == read_csv == 0:#
#		sys.exit(colorama.Fore.RED + 'Cannot get sequences, write a FASTA, or write a VEP file without being given a csv or sam file')
if (writeFASTA == 1 or writevepfile ==1 or writelogo ==1 or getannotations==1) and getseqs == 0:
	sys.exit(colorama.Fore.RED + "Can't write a fasta file, vep file, or logo file, or get annotations without getting sequences.")
if sequencesource != 'LOCAL' and sequencesource != 'REMOTE' and getseqs == 1:
	sys.exit(colorama.Fore.RED + "Must specify 'LOCAL' or 'REMOTE' source to get sequences")

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

if userandomIS or userandomSEQS: 
	random.seed(a=None, version=2) #initialized random number generator w/ no seed (use system time)

if sequencesource =='LOCAL':
#	GRCh38=Bio.SeqIO.index('../genomes/GRCh38.fasta',"fasta")
	genome=twobitreader.TwoBitFile(twobitgenomelocation) #reads genome into object
else:
	genome=None

#trim primers and adapters from raw fastqreads, return trimmed "genomic" sequences in fasta format.
if inputtype=="fastq":
	FASTQ.Trim(inputfile,trimmedfastq,barcode, primer5, primer3, trim3, trim5, minimum_read_len)
	#FASTQ.QtoA(trimmedfastq, f'{rootname}.fasta', Ltrim=0, trim=0)
	if userandomSEQS: #replace all real reads with random NT
		FASTQ.randomize(trimmedfastq,format='fastq')
	if mapreads: #run bowtie using trimmed fastq and quality scores (--phred33)
		bowtiecommand=f'{bowtielocation} --phred33 -p 4 -x {bowtieindex} -U {trimmedfastq} -S {rootname}.sam' #2>&1 | tee {rootname}_bowtie.log'
		print(f'mapping reads genome using bowtie2. Writing output to {rootname}.sam')
		print(bowtiecommand)
		bowtie=runbin.Command(bowtiecommand)
		bowtieout=bowtie.run(timeout=20000)
		bowtieout[1]
		print(bowtieout[1].decode())
		print(bowtieout[2].decode())	


#map reads from fasta file to genome, return sam file with genome locations
elif mapreads: #if mapreads, but not using fastq, run bowtie specifying fasta (-f)
	fastareads=f'{rootname}.fasta'
	if userandomSEQS:
		FASTQ.randomize(fastareads,format='fasta')
	#add bowtie mapping function to map plasmid sequences?
	bowtiecommand=f'{bowtielocation} -f -x {bowtieindex} -p 4 -U {fastareads} -S {rootname}.sam' #2>&1 | tee {rootname}_bowtie.log'
	print(f'mapping reads genome using bowtie2. Writing output to'+colorama.Fore.YELLOW+f' {rootname}.sam'+colorama.Style.RESET_ALL+f' and '+colorama.Fore.YELLOW+f'{rootname}_bowtie.log'+colorama.Style.RESET_ALL)
	print(colorama.Fore.CYAN+f'{bowtiecommand}'+colorama.Style.RESET_ALL)
	bowtie=runbin.Command(bowtiecommand)
	bowtieout=bowtie.run(timeout=20000)
	print(bowtieout[1].decode())
	print(colorama.Fore.RED+f'{bowtieout[2].decode()}'+colorama.Style.RESET_ALL)

if inputtype=="sam" or mapreads:
	readslist, unmapped = mappedreads.read_sam(sam_file,chromIDS,ISsamfilename, compressreads=compressreads,random=userandomIS)

if write_csv:
	mappedreads.write_csv(readslist,genome_location_csv)

#Chrom,Sense,Loc,Gene,Total # IS Sequences Found,Total # IS Found,Plasmid m995,Sample1,Sample2,Sample3,Sample4

#read_csv should always be followed by getseqs, which will genererate readslist
#if read_csv:
#	readslist=mappedreads.read_csv(genome_location_csv)#reads genome location csv file

if getseqs:
	recordlist, readslist=getseq.get_seqs(readslist, lwindow, rwindow, samwindow,genome, source=sequencesource)

if writeFASTA:#write FASTA formatted file with returned sequences
	print(colorama.Style.RESET_ALL + f'\rWriting fasta file:'+ colorama.Fore.GREEN + f' {FASTAfile}')
	Bio.SeqIO.write(recordlist, FASTAfile, "fasta")
if writelogo: #dependant on having sequences, optional to make logo plot
	weblogocommand=f'{weblogolocation} -f {FASTAfile} -D fasta -o {logofile} -F svg -A dna -F png --resolution 600 -s large -c classic -i {str(int(1-1*lwindow))} -l -10 -u 10 '#-l [lower bound] -u [upper bound]
	print(colorama.Style.RESET_ALL + f'Writing logo')
	print(colorama.Fore.YELLOW + f'{weblogocommand}')
	logocommand=runbin.Command(weblogocommand)
	logoout=logocommand.run(timeout=1800) 

if writevepfile:
	reads = sorted(readslist, key = lambda x: (x.chrom, x.loc))
	print(colorama.Style.RESET_ALL + f'Writing file for VEP analysis: '+colorama.Fore.YELLOW + f'{outputVEPfile}')
	with open(outputVEPfile, 'w', newline='') as formatted_csv_file: #writing file formatted for VEP analysis
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
