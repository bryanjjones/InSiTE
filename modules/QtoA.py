#!/home/notal/anaconda3/bin/python3
#!/usr/bin/python3
#import _mypath
import os
import time
import Bio.SeqIO
import Bio.Entrez
import Bio.SeqRecord
import Bio.Seq
import sys
import csv
#import runbin
import random
#import time
#import StringIO
#import ConfigParser
#import httplib

def QtoA(inputfile, outputfile, Ltrim=0, trim=10000):
	seqs=[]#Bio.SeqIO.read(inputfile, "fastq")
	print(f'reading {inputfile}')
	for record in Bio.SeqIO.parse("./RAWfastq/"+inputfile, "fastq"):
		if trim:
			trimmed=record.seq[Ltrim:trim]
			record.letter_annotations = {}
			record.seq=trimmed
			print(f'{record.seq[:20]}...{record.seq[-20:]}  sequnce number {len(seqs):,}', end="\r")
		seqs.append(record)
	print(f'writing {outputfile}')
	Bio.SeqIO.write(seqs, "./inprocess/"+outputfile, "fasta")
	return seqs

filelist=("../../../../mnt/c/Users/Bryan.Jones/Documents/05genomiconly.fastq",)
'''
"../../../../mnt/c/Users/Bryan.Jones/Documents/05.fastq"
filelist=("ILM170_04x07x2019_CL_BMO01GW_CL0174u01_cPCR_LMPCR_CovarisLM_Sample1_24x04x2019_TcellgDNA_ILM170_ILM170_CL0175u01_1000_MiS5TXIR4x159_AGTATCTCGT_GGGTTCCGCCGGATGGC_hg_0_CCTAACTGCTGTGCCACT_TCTACTGAGTGA_TAGTCACT.fastq",
"ILM170_04x07x2019_CL_BMO01GW_CL0174u02_cPCR_LMPCR_CovarisLM_Sample1_24x04x2019_TcellgDNA_ILM170_ILM170_CL0175u02_1000_MiS5TXIR4x159_AGTATCTCGT_GGGTTCCGCCGGATGGC_hg_0_CCTAACTGCTGTGCCACT_TGCATCTAGCTC_TAGGAGCT.fastq",
"ILM170_04x07x2019_CL_BMO01GW_CL0174u03_cPCR_LMPCR_CovarisLM_Sample1_24x04x2019_TcellgDNA_ILM170_ILM170_CL0175u03_1000_MiS5TXIR4x159_AGTATCTCGT_GGGTTCCGCCGGATGGC_hg_0_CCTAACTGCTGTGCCACT_TGTAGTCTACTG_TAGCAGTA.fastq",
"ILM170_04x07x2019_CL_BMO01GW_CL0174u04_cPCR_LMPCR_CovarisLM_Sample2_24x04x2019_TcellgDNA_ILM170_ILM170_CL0175u04_1000_MiS5TXIR4x160_ATCTGCGACG_GGGTTCCGCCGGATGGC_hg_0_CCTAACTGCTGTGCCACT_ATCATCACACTC_TAGGAGTG.fastq",
"ILM170_04x07x2019_CL_BMO01GW_CL0174u05_cPCR_LMPCR_CovarisLM_Sample2_24x04x2019_TcellgDNA_ILM170_ILM170_CL0175u05_1000_MiS5TXIR4x160_ATCTGCGACG_GGGTTCCGCCGGATGGC_hg_0_CCTAACTGCTGTGCCACT_CGTGCATAGCTG_TAGCAGCT.fastq",
"ILM170_04x07x2019_CL_BMO01GW_CL0174u06_cPCR_LMPCR_CovarisLM_Sample2_24x04x2019_TcellgDNA_ILM170_ILM170_CL0175u06_1000_MiS5TXIR4x160_ATCTGCGACG_GGGTTCCGCCGGATGGC_hg_0_CCTAACTGCTGTGCCACT_TGCTATACGTGC_TAGGCACG.fastq",
"ILM170_04x07x2019_CL_BMO01GW_CL0174u07_cPCR_LMPCR_CovarisLM_Sample3_24x04x2019_TcellgDNA_ILM170_ILM170_CL0175u07_1000_MiS5TXIR4x161_CTCTCTGATG_GGGTTCCGCCGGATGGC_hg_0_CCTAACTGCTGTGCCACT_CTGTCTGACGAG_TAGCTCGT.fastq",
"ILM170_04x07x2019_CL_BMO01GW_CL0174u08_cPCR_LMPCR_CovarisLM_Sample3_24x04x2019_TcellgDNA_ILM170_ILM170_CL0175u08_1000_MiS5TXIR4x161_CTCTCTGATG_GGGTTCCGCCGGATGGC_hg_0_CCTAACTGCTGTGCCACT_TATAGTAGTCTC_TAGGAGAC.fastq",
"ILM170_04x07x2019_CL_BMO01GW_CL0174u09_cPCR_LMPCR_CovarisLM_Sample3_24x04x2019_TcellgDNA_ILM170_ILM170_CL0175u09_1000_MiS5TXIR4x161_CTCTCTGATG_GGGTTCCGCCGGATGGC_hg_0_CCTAACTGCTGTGCCACT_CTGCGTGACAGC_TAGGCTGT.fastq",
"ILM170_04x07x2019_CL_BMO01GW_CL0174u10_cPCR_LMPCR_CovarisLM_Sample4_24x04x2019_TcellgDNA_ILM170_ILM170_CL0175u10_1000_MiS5TXIR4x162_TATATAGCAC_GGGTTCCGCCGGATGGC_hg_0_CCTAACTGCTGTGCCACT_CATGTAGAGACT_TAGAGTCT.fastq",
"ILM170_04x07x2019_CL_BMO01GW_CL0174u11_cPCR_LMPCR_CovarisLM_Sample4_24x04x2019_TcellgDNA_ILM170_ILM170_CL0175u11_1000_MiS5TXIR4x162_TATATAGCAC_GGGTTCCGCCGGATGGC_hg_0_CCTAACTGCTGTGCCACT_CTCAGAGATGAT_TAGATCAT.fastq",
"ILM170_04x07x2019_CL_BMO01GW_CL0174u12_cPCR_LMPCR_CovarisLM_Sample4_24x04x2019_TcellgDNA_ILM170_ILM170_CL0175u12_1000_MiS5TXIR4x162_TATATAGCAC_GGGTTCCGCCGGATGGC_hg_0_CCTAACTGCTGTGCCACT_AGTCTCAGCACA_TAGTGTGC.fastq",
"ILM170_04x07x2019_CL_BMO01GW_CL0174u13_cPCR_LMPCR_CovarisLM_Plasmidm995_0_Plasmid_ILM170_ILM170_CL0175u13_1_MiS5TXIR4x163_CTCACGATCT_GGGTTCCGCCGGATGGC_Ctrl_0_CCTAACTGCTGTGCCACT_TCAGCGATCGAT_TAGATCGATCGCTGAC.fastq",
"ILM170_04x07x2019_CL_BMO01GW_CL0174u14_cPCR_LMPCR_CovarisLM_hgDNA_0_PBMC_ILM170_ILM170_CL0175u14_1000_MiS5TXIR4x164_TACTCTGCGT_GGGTTCCGCCGGATGGC_hg_Ctrl_CCTAACTGCTGTGCCACT_TCTGCGAGAGAG_TAGCTCTCTCGCAGACCTAAC.fastq")
outputfiles=("01_full.fasta","02_full.fasta","03_full.fasta","04_full.fasta","05_full.fasta","06_full.fasta","07_full.fasta","08_full.fasta","09_full.fasta","10_full.fasta","11_full.fasta","12_full.fasta","13_full.fasta","14_full.fasta")
../../../../mnt/c/Users/Bryan.Jones/Documents/05.fasta",
'''
outputfiles=("../../../../mnt/c/Users/Bryan.Jones/Documents/05genomiconly.fasta",)
for i in range(len(filelist)):
	QtoA(filelist[i],outputfiles[i])

exit()
'''
#Various functions:
#getseqs=0#1#download sequences from entrez
#writeFASTA=0#1#write a fasta file of all sequences around IS
writelogo=0#1#create a logo image of consensus sequence # requires getseqs and writeFASTA to be on
writevepfile=0#write a vep file to use with vep, either online or by setting getannotations to 1
getannotations=1#1 # uses VEP to get anotations for insertion sites
userandomIS=0#0 #useful for generating a random control data set chromosome matched to query set.

window=50 #window on either side of indicated nt location to return
#specified file names
inputfile='IS_data.csv'#'IS_data.csv'
outputfile="InsertionSequences.fasta"
outputVEPfile="rand_VEPfile.csv"
logofile="IS_logo.svg"
annotationsfile='rand_IS_annotations.csv'

veplocation = './ensembl-vep/vep'
weblogolocation = 'weblogo' #in PATH
Bio.Entrez.email = "bryan@bmogen.com"     # Always tell NCBI who you are
Bio.Entrez.api_key = "679368bc3fef47bb440fb743c889befe4e09" #bryan's personal NCBI key (allows faster requests of their server)
#maxtries=5
reads=[]
recordlist=[]        
chromIDS={}
chromNTS={}

#make a dictionary of chromosome names and ID #s usich chromosomes.csv, which is the table copied from https://www.ncbi.nlm.nih.gov/genome?term=human&cmd=DetailsSearch
with open('chromosomes.csv') as csv_file:
	csv_reader = csv.reader(csv_file, delimiter=",")
	line_count=0
	for row in csv_reader:
		name=row[1]
		ID=row[2]
		nts=row[4]
		print(f'chromosome {name} is called {ID} and has {nts} nucleotides')
		chromIDS[name]=ID
		chromNTS[name]=nts
if userandomIS:
	random.seed(a=None, version=2)
#parse a row from the read file and pull the sequence
class read(object):
	def __init__(self, row):
		tries = 1
		self.chrom=int(row[0])
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
		if getseqs: #if get sequences, fetch sequence centered on given location, with a window on either side
			handle = Bio.Entrez.efetch(db="nucleotide", #fetches the specified sequence as a fasta
				id=chromIDS[self.chrom], #id for indicated chromosome from refrence genome assembly
				rettype="fasta", #format
				strand=self.sensenum, #which strand
				seq_start=str(self.loc-window), #range of sequence to return
				seq_stop=str(self.loc+window))
			self.record = Bio.SeqIO.read(handle, "fasta") #reads the returned fasta into a sequence object
			handle.close()#close efetch 
			self.sequence=self.record.seq #adds the sequence property as a string 
		else: #if not getting sequences, just put empt bits.
			self.record=None
			self.sequence="-"
if getseqs or writeFASTA or writevepfile 
	print(f'reading data from {inputfile}')
	with open(inputfile) as csv_file: #open file
		csv_reader = csv.reader(csv_file, delimiter=',') #read lines
		line_count = 0
		for row in csv_reader:
			if line_count < 4:#first 4 lines are headders
				print(f'Column names are {", ".join(row)}')
				line_count += 1
			else: #read data from each line
				thisrow=read(row)
				reads.append(thisrow)
				recordlist.append(thisrow.record)
				print(f'\t {line_count} read on chromosome {thisrow.chrom} with {thisrow.sense} sense at location {thisrow.loc}. and has sequence {thisrow.sequence}')
				line_count += 1 
		print(f'Processed {line_count} lines.')

if writeFASTA:#write FASTA formatted file with returned sequences
	print(f'writting fasta file: {outputfile}')
	Bio.SeqIO.write(recordlist, outputfile, "fasta")
	if writelogo: #dependant on having sequences, optional to make logo plot
		weblogocommand=f'{weblogolocation} -f {outputfile} -D fasta -o {logofile} -F svg -A dna -F png --resolution 600 -s large -c clasic -i '+str(int(-1*window)) #-l [lower bound] -u [upper bound]
		print('writting logo')
		print(weblogocommand)
		logocommand=runbin.Command(weblogocommand)
		logoout=logocommand.run(timeout=1800) 

if writevepfile:
	reads = sorted(reads, key = lambda x: (x.chrom, x.loc))
	print(f'writing file for VEP analysis: {outputVEPfile}')
	with open(outputVEPfile, 'w', newline='') as formatted_csv_file: #writting file formatted for VEP analysis
		csv_writer = csv.writer(formatted_csv_file, delimiter='\t') #tab delimited csv file
		for entry in reads:
			if entry.chrom == 23:
				entry.chrom ="X"
			elif entry.chrom == 24:
				entry.chrom ="Y"
			csv_writer.writerow([entry.chrom, str(entry.loc+1), entry.loc, "-/A", entry.sense])

if getannotations:
	vepcommand=f'{veplocation} -i {outputVEPfile} -o {annotationsfile} --buffer_size 100000 --merged --format "ensembl" --cache --nearest gene --distance 1 --offline --use_given_ref --fork 4 --pick' #runs vep using "--merged" ensembl and refseq genomes (what i currently have "--cached" on my machine). return the "--nearest gene". --distance from query to look for genes. maybe add --symbol --protein --ccds. Maybe add back: --numbers --domains --biotype --hgvs add --merged --gencode_basci
	print('running VEP analysis to annotate genomic positions using ensembl vep. Writting output to {annotationsfile} and {annotationsfile}_summary.html')
	print(vepcommand)
	vep=runbin.Command(vepcommand)
	vepout=vep.run(timeout=20000)
	print(vepout[1]) #STDOUT
	print(vepout[2]) #ERRORS

#Do I need this or will teh VEP report do that for me? read consequence, and classify simply (intron/exon/intragenic/UTR/regulatory)
print(f'Completed all tasks. Exiting.')
exit()
'''