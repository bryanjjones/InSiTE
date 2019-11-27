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
import argparse
import logging

#inputs:
inputtype='fastq' #sam, csv, fasta, fastq
inputfile='05.fastq' #input file name
userandomIS=0#0 #useful for generating a random control data set chromosome matched to query set.
userandomSEQS=0 #generate random sequence data, usefull for calculating false positive rate when mapping reads w/ bowtie

#tasks:
getseqs=1#1#get sequences from either entrez or local TwoBit genome around locations indicated by genome_location_csv
compressreads=1 #remove duplicate reads and reads shifted +/- 1 nt, number of reads are compressed in csv and fasta
getannotations=1#1 # use annotate.py module to map insertion sites to anotations 

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
featurenames=['intron', 'exon','codingexon','transcript','TSS'] #feature names found in feature files to map reads to
featuredist=[False,False,False,False,True]# weather to map distance of each read, or only whether reads overlap with feature
distance=1000 # distance in bp to be considered close to feature

#outputs:
write_csv=1 # write csv (only relevant if reading from sam file)
writeFASTA=1#1#write a fasta file of all sequences around IS
writelogo=1#1#create a logo image of consensus sequence # requires getseqs and writeFASTA to be on
#writevepfile=0#write a vep file to use with vep, either online or by setting getannotations to 1

#refrence file locations
chromosome_ids = './refrence_datasets/chromosomes.csv'
#veplocation = './ensembl-vep/vep'
bowtielocation = 'bowtie2'
bowtieindex = './refrence_datasets/genomes/GRCh38.fna.bowtie_index/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index'
weblogolocation = 'weblogo' #in PATH
twobitlocation = './scripts/TwoBitToFa'
twobitgenomelocation='./refrence_datasets/genomes/GRCh38.2bit'
annotations=['./refrence_datasets/annotations/gencode.v32.introns.bed','./refrence_datasets/annotations/gencode.v32.exons.bed','./refrence_datasets/annotations/gencode.v32.codingexons.bed','./refrence_datasets/annotations/gencode.v32.transcripts.gtf','./refrence_datasets/annotations/gencode.v32.transcripts.gtf'] #annotations file in gff, gtf, or bed format,
if __name__ == "__main__":

	ap = argparse.ArgumentParser(prog=os.path.basename(sys.argv[0]),
	                             usage=__doc__)
	ap.add_argument('-q', '--fastq', help='fastq input file')
	ap.add_argument('-r', '--rand_is', default=False, action='store_true', 
	                help='use random integration sites matched to given query set instead of actual query set')
	ap.add_argument('-n', '--rand_nt', default=False, action='store_true',
	                help='use random sequences for mapping')
	ap.add_argument('-c','--csv', help='csv input file')
	ap.add_argument('-a','--fasta', help='fasta input file')
	ap.add_argument('-s','--sam', help='sam/bam input file')
	ap.add_argument('--no_seqs', default=False, action='store_true',
	                help='do not get sequences from either entrez or local TwoBit genome around locations indicated by genome_location_csv')
	ap.add_argument('-z', '--compress_reads', action='store_true', help='compress duplicate reads and reads shifted +/- 1nt, number of reads are compressed in csv, fasta, and mapping outputs')
	ap.add_argument('--no_annotate', default=False, action='store_true', help='do not map insertion sites to genome annotations')
	ap.add_argument('--barcode', default='ATCTGCGACG', help=' barcode sequence to trim off of reads')
	ap.add_argument('--lwindow', default=50, help='numebr of nucleotides upstream of integration site to return', type=int) #window on either side of indicated nt location to return
	ap.add_argument('--rwindow', default=50, help='number of nucleotides downstream of integration site to return', type=int) #window on either side of indicated nt location to return
	ap.add_argument('--samwindow', default=0, type=int, help='depreciated')#60
	ap.add_argument('--remote', action='store_true', help="get sequences from entrez server instead of 'LOCAL' TwoBit genome")
	ap.add_argument('--min', default=25, type=int, help='minimum length of a read to try mapping. default (25) will usually avoid any false positives in read sets of 200k reads')
	ap.add_argument('--primer5', default='GGGTTCCGCCGGATGGC', help="5' primer sequenc to remove from reads") 
	ap.add_argument('--primer3', default='CCTAACTGCTGTGCCACT', help="3' primer sequence to remove from reads")
	ap.add_argument('--trim5', default=21, type=int, help="additional (non-genomic) nts to trim off of 3' end of reads")
	ap.add_argument('--trim3', default=16, type=int, help="additional (non-genomic) nts to trim off of 5' end of reads (starts with TA)")
	ap.add_argument('--feature', action='append', help='feature names found in feature files to map reads to')#['intron', 'exon','codingexon','transcript','TSS'] #feature names found in feature files to map reads to
	ap.add_argument('--dist', action='append', help='weather to map distance of each read, or only whether reads overlap with feature. same number of distance variables must be given as features.')
	ap.add_argument('--close', default=1000, type=int, help='distance in bp to be considered close to feature')
	ap.add_argument('--chromosome_ids', default='./refrence_datasets/chromosomes.csv')
	ap.add_argument('--bowtielocation', default = 'bowtie2')
	ap.add_argument('--bowtieindex', default = './refrence_datasets/genomes/GRCh38.fna.bowtie_index/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index')
	ap.add_argument('--weblogolocation', default = 'weblogo')#in PATH
	ap.add_argument('--twobitlocation', default = './scripts/TwoBitToFa')
	ap.add_argument('--twobitgenomelocation', default='./refrence_datasets/genomes/GRCh38.2bit')
	ap.add_argument('--annotations', action='append', help='location of annotation file(s) (bed/gff/gtf), must be same number of files as features specified')#['./refrence_datasets/annotations/gencode.v32.introns.bed','./refrence_datasets/annotations/gencode.v32.exons.bed','./refrence_datasets/annotations/gencode.v32.codingexons.bed','./refrence_datasets/annotations/gencode.v32.transcripts.gtf','./refrence_datasets/annotations/gencode.v32.transcripts.gtf']
	ap.add_argument('--supress_csv', default=False, action='store_true')
	ap.add_argument('--supress_fasta',default=False, action='store_true')
	ap.add_argument('--supress_logo',default=False, action='store_true')
	ap.add_argument('--append_summary',default=False, action='store_true')
	args = ap.parse_args()

	if args.fastq:
	 	inputfile = args.fastq
	 	inputtype ='fastq'
	if args.csv:
		inputfile = args.csv
		inputtype = 'csv'
	if args.fasta:
		inputfile = args.fasta
		inputtype = 'fasta'
	if args.sam:
		inputfile = args.sam
		inputtipe = 'sam'
	if not bool(args.fasta)+bool(args.fastq)+bool(args.csv)+bool(args.sam) == 1:
		loggging.critical('please only provide a single input file (fasta/fastq/sam/bam/csv)')
		sys.exit('please only provide a single input file (fasta/fastq/sam/bam/csv)')
	userandomIS = args.rand_is
	userandomSEQS = args.rand_nt
	getseqs = not args.no_seqs
	compressreads = args.compress_reads #remove duplicate reads and reads shifted +/- 1 nt, number of reads are compressed in csv and fasta
	getannotations = not args.no_annotate
	barcode = args.barcode
	lwindow = args.lwindow
	rwindow = args.rwindow
	samwindow = args.samwindow
	if args.remote:
		sequencesource = 'REMOTE'
	else:
		sequencesorce = 'LOCAL'
	minimum_read_len = args.min
	primer5 = args.primer5
	primer3 = args.primer3
	trim5 = args.trim5
	trim3 = args.trim3
	if args.feature:
		if not len(args.feature)==len(args.dist)==len(args.annotations):
			logging.critical('must provide same number of features as distance and annotation files as they correspond to each other')
			sys.exit('must provide same number of features as distance and annotation files as they correspond to each other')
		featurenames = args.feature
		featuredist = args.dist
		annotations = args.annotations
	distance = args.close
	chromosome_ids = args.chromosome_ids
	bowtielocation = args.bowtielocation
	bowtieindex = args.bowtieindex
	weblogolocation = args.weblogolocation
	twobitlocation = args.twobitlocation
	twobitgenomelocation = args.twobitgenomelocation
	write_csv = not args.supress_csv
	writeFASTA = not args.supress_fasta
	writelogo = not args.supress_logo

#specified file names
rootname=os.path.splitext(inputfile)[0]

#output file names
sam_file=f'{rootname}.sam'#,'02-wind20-70.sam','03-wind20-70.sam','04-wind20-70.sam','05-wind20-70.sam','06-wind20-70.sam','07-wind20-70.sam','08-wind20-70.sam','09-wind20-70.sam','10-wind20-70.sam','11-wind20-70.sam','12-wind20-70.sam','13-wind20-70.sam',]
genome_location_csv=f'{rootname}_IS_mappings.csv'#'IS_data.csv'
trimmedfastq=f'{rootname}_trimmed.fastq' # fastq file with adapters/primers/barcodes trimmed off with cutadapt
FASTAfile=f'{rootname}_retrieved_2bit.fasta'#output fasta file containing sequences surrounding mapped insertion site
#outputVEPfile=f"{rootname}_VEPfile.csv" #VEP file containing mapped locations of insertions in VEP format so VEP can provide annotations
logofile=f"{rootname}IS_logo.svg" # logo file showing consensus integration site in logo format
annotationsfile=f'{rootname}_IS_annotations.csv' # file contoining summary of mapping to annotations
distancesfile=f'{rootname}_distances.csv' # file containing list of distances of each read to nearest TSS
logfile=f'{rootname}.log'
ISbamfilename=f'{rootname}IS.bam'#sam file name for single nt IS mappings

logging.basicConfig(filename=logfile, filemode='w', level=logging.DEBUG)

#sanity checks:
if inputtype=="sam":
	read_sam_file=1
	mapreads=0
elif inputtype=="csv":
	read_csv=1
	mapreads=0
	if getseqs==0:
		logging.critical('Given a csv file, but not asked to get seqs. Nothing to do.')
		sys.exit(colorama.Fore.RED + 'Given a csv file, but not asked to get seqs. Nothing to do.')
elif inputtype=="fasta":
	read_fasta=1
	mapreads=1
elif inputtype=="fastq":
	read_fastq=1
	mapreads=1
else:
	logging.critical("Specify input file type as 'sam', 'csv', 'fasta', or 'fastq'")
	sys.exit(colorama.Fore.RED + 'Specify input file type as \'sam\', \'csv\', \'fasta\', or \'fastq\'')
#if getseqs or writeFASTA or writevepfile:
#	if read_sam_file == read_csv == 0:#
#		sys.exit(colorama.Fore.RED + 'Cannot get sequences, write a FASTA, or write a VEP file without being given a csv or sam file')
if (writeFASTA == 1 or writelogo ==1 or getannotations==1) and getseqs == 0: #or writevepfile ==1 
	logging.critical("Can't write a fasta file or logo file, or get annotations without getting sequences.")
	sys.exit(colorama.Fore.RED + "Can't write a fasta file or logo file, or get annotations without getting sequences.") #, vep file, 
if sequencesource != 'LOCAL' and sequencesource != 'REMOTE' and getseqs == 1:
	logging.critical("Must specify 'LOCAL' or 'REMOTE' source to get sequences")
	sys.exit(colorama.Fore.RED + "Must specify 'LOCAL' or 'REMOTE' source to get sequences")
if getannotations and not len(annotations)==len(featurenames)==len(featuredist):
	logging.critical("If getting annotations, length of annotation files, feature names, and feature distance must be the same")
	sys.exit(colorama.Fore.RED + "If getting annotations, length of annotation files, feature names, and feature distance must be the same")
summary=[rootname,barcode]
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
	trimlog=FASTQ.Trim(inputfile,trimmedfastq,barcode, primer5, primer3, trim3, trim5, minimum_read_len)
	logging.info(trimlog)
	summary.append('raw sequences')
	summary.append(int(len(open(inputfile).readlines())/4))#4 lines per entry in fastq
	summary.append('trimmed sequences')
	summary.append(int(len(open(trimmedfastq).readlines())/4))
	#FASTQ.QtoA(trimmedfastq, f'{rootname}.fasta', Ltrim=0, trim=0)
	if len(open(trimmedfastq).readlines())==0:
		logging.critical(f'No sequences left after trimming. Exiting')
		sys.exit(colorama.Fore.RED+f'No sequences left after trimming. Exiting'+colorama.Style.RESET_ALL)
	if userandomSEQS: #replace all real reads with random NT
		trimlog=FASTQ.randomize(trimmedfastq,format='fastq')
		logging.info(trimlog)
	#run bowtie using trimmed fastq and quality scores (--phred33)
	bowtiecommand=f'{bowtielocation} --phred33 -p 4 -x {bowtieindex} -U {trimmedfastq} -S {rootname}.sam' #2>&1 | tee {rootname}_bowtie.log'
	print(f'mapping reads genome using bowtie2. Writing output to {rootname}.sam')
	print(bowtiecommand)
	logging.info(f'mapping reads genome using bowtie2. Writing output to {rootname}.sam')
	logging.info(bowtiecommand)
	bowtie=runbin.Command(bowtiecommand)
	bowtieout=bowtie.run(timeout=20000)
	bowtieout[1]
	print(bowtieout[1].decode())
	print(bowtieout[2].decode())	
	logging.info(bowtieout[1].decode())
	logging.info(bowtieout[2].decode())	


#map reads from fasta file to genome, return sam file with genome locations
elif mapreads: #if mapreads, but not using fastq, run bowtie specifying fasta (-f)
	fastareads=f'{rootname}.fasta'
	if userandomSEQS:
		FASTQ.randomize(fastareads,format='fasta')
	#add bowtie mapping function to map plasmid sequences?
	bowtiecommand=f'{bowtielocation} -f -x {bowtieindex} -p 4 -U {fastareads} -S {sam_file}' #2>&1 | tee {rootname}_bowtie.log'
	print(f'mapping reads genome using bowtie2. Writing output to'+colorama.Fore.YELLOW+f' {sam_file}'+colorama.Style.RESET_ALL+f' and '+colorama.Fore.YELLOW+f'{rootname}_bowtie.log'+colorama.Style.RESET_ALL)
	print(colorama.Fore.CYAN+f'{bowtiecommand}'+colorama.Style.RESET_ALL)
	logging.info(f'mapping reads genome using bowtie2. Writing output to'+f' {sam_file}'+f' and '+f'{rootname}_bowtie.log')
	logging.info(f'{bowtiecommand}')
	bowtie=runbin.Command(bowtiecommand)
	bowtieout=bowtie.run(timeout=20000)
	logging.info(bowtieout[1].decode())
	print(bowtieout[1].decode())
	logging.error(bowtieout[2].decode())
	print(colorama.Fore.RED+f'{bowtieout[2].decode()}'+colorama.Style.RESET_ALL)

if inputtype=="sam" or mapreads:
	readslist, unmapped, message = mappedreads.read_sam(sam_file,chromIDS,ISbamfilename, compressreads=compressreads,chromNTS=chromNTS,randomize=userandomIS)
	logging.info(message)
if write_csv:
	message=mappedreads.write_csv(readslist,genome_location_csv)
	logging.info(message)

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
	print(colorama.Fore.YELLOW + f'{weblogocommand}'+colorama.Style.RESET_ALL)
	logging.info( f'Writing logo')
	logging.info( f'{weblogocommand}')
	logocommand=runbin.Command(weblogocommand)
	logoout=logocommand.run(timeout=1800) 
'''
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
'''
if getannotations and len:
	with open(annotationsfile, 'w', newline='') as output_file:
		totalprinted=False
		output_writer = csv.writer(output_file, delimiter=",")
		for i in range(len(featurenames)):
			if featuredist[i]:
				print(colorama.Style.RESET_ALL+f'mapping insertion site distances to '+colorama.Fore.YELLOW+f'{featurenames[i]}'+colorama.Style.RESET_ALL+f' in '+colorama.Fore.YELLOW+f'{annotations[i]}'+colorama.Style.RESET_ALL)
				logging.info(f'mapping insertion site distances to '+f'{featurenames[i]}'+f' in '+f'{annotations[i]}')
				distances, average, standarddev, close, b , message= annotate.closest (ISbamfilename,annotations[i],featurename=featurenames[i],limit=distance,position="start")
				logging.info(message)
				try:
					output_writer.writerow([f'{featurenames[i]}(average distance)',average])
					output_writer.writerow([f'{featurenames[i]}(standard deviation)',standarddev])
					output_writer.writerow([f'integration events within {distance} bp of {featurenames[i]}',close])
					summary.append(featurenames[i])
					summary.append('average distanance')
					summary.append(average)
					summary.append('st dev')
					summary.append(standarddev)
					summary.append(f'within {distance}bp')
					summary.append(close)
				except:
					logging.error("can't write stats. Maybe there were no sequences")
					print("can't write stats. Maybe there were no sequences")
				with open(distancesfile, 'w', newline='') as distance_file:
					dist_writer = csv.writer(distance_file, delimiter=",")
					dist_writer.writerow(distances)
			else:
				print(colorama.Style.RESET_ALL+f'mapping insertion sites to'+colorama.Fore.YELLOW+f' {featurenames[i]}'+colorama.Style.RESET_ALL+f' in '+colorama.Fore.YELLOW+f'{annotations[i]}'+colorama.Style.RESET_ALL)
				logging.info(f'mapping insertion sites to'+f' {featurenames[i]}'+f' in '+f'{annotations[i]}')
				results, message = annotate.featuremap (annotations[i], ISbamfilename, featurenames=featurenames[i], single=True ,procs=3)
				if not totalprinted: #featuremap returns total reads, add this line only once as it should be the same for each feature.
					output_writer.writerow([f'total reads', results[1]])
					summary.append('reads mapped to genome')
					summary.append(results[1])
					totalprinted=True
				output_writer.writerow([f'reads mapped to {featurenames[i]}', results[0]])
				summary.append(featurenames[i])
				summary.append(results[0])
				logging.info(message)
				
	'''
	vepcommand=f'{veplocation} -i {outputVEPfile} -o {annotationsfile} --force_overwrite --buffer_size 100000 --merged --format "ensembl" --cache --nearest gene --distance 1 --offline --use_given_ref --fork 4 --pick' #runs vep using "--merged" ensembl and refseq genomes (what i currently have "--cached" on my machine). return the "--nearest gene". --distance from query to look for genes. maybe add --symbol --protein --ccds. Maybe add back: --numbers --domains --biotype --hgvs add --merged --gencode_basci
	print(colorama.Style.RESET_ALL + f'Running VEP analysis to annotate genomic positions using ensembl vep.'+colorama.Fore.YELLOW + f'"{vepcommand}"')
	vep=runbin.Command(vepcommand)
	vepout=vep.run(timeout=20000)
	if len(vepout[1])>0:
		print(vepout[1]) #STDOUT
	if len(vepout[2])>0:
		print(vepout[2]) #ERRORS
	'''
programend = time.time()
message=(colorama.Fore.GREEN + f'Completed all tasks in {int(programend-programstart)} seconds. Exiting.'+colorama.Style.RESET_ALL)
logging.info(f'Completed all tasks in {int(programend-programstart)} seconds. Exiting.')
if args.append_summary:
	with open('./summary.csv','a') as summary_file:
		csv.writer(summary_file, delimiter=',').writerow(summary)

print(message)
exit()
