#!/btapps/miniconda3/bin/python3
import csv
import pysam
import colorama
import Bio
import Bio.Entrez
import Bio.SeqIO
import random
from operator import attrgetter

#read take a given string, formatted as row read from a csv file, and return an object with useful properties
class read_csv_line(object):
	def __init__(self, row):#, userandomIS=False, chromNTS={}):
		#tries = 1
		self.chrom=str(row[0])
		self.sense=row[1]
		if row[1]=="+":
			self.sensenum=1
		elif row[1]=="-":
			self.sensenum=2
		#if userandomIS: #if random IS used, per "read" use a random number for insertion site
		#	self.loc=random.randrange(int(chromNTS[str(self.chrom)])) # random int returned in the range of the given chromosome range
		#else: #if not random, use given values
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

#if compressedreads:
#compress a sorted readlist by merging identical locations and, if adjacent=True, merging adjacent locations too (merges to location with greatest count)
def compress(readslist,adjacent=True):
	semicompressed_readslist=[]
	compressed_readslist=[]
	for i in range(len(readslist)):
		if i+1<len(readslist) and readslist[i].loc==readslist[i+1].loc and readslist[i].chrom == readslist[i+1].chrom: #if this read isn't the last one and if it's the same as the next one
			readslist[i+1].totalseqs=readslist[i+1].totalseqs+readslist[i].totalseqs #if they are the same chromosome and locus, add to the total
		else:
			semicompressed_readslist.append(readslist[i]) #if it's different than the next one, add it to the list (with the total count)
	#adjust loc +/- 1 to account for indels
	if adjacent:
		for i in range(len(semicompressed_readslist)):
			if i+1<len(semicompressed_readslist) and semicompressed_readslist[i].loc==(semicompressed_readslist[i+1].loc-1) and semicompressed_readslist[i].chrom == semicompressed_readslist[i+1].chrom: #if this read isn't the last one and if it's the same as the next one
				print(colorama.Style.RESET_ALL + f'merging adjacent read locations. read '+colorama.Fore.YELLOW + f'{i}'+colorama.Style.RESET_ALL + ', locations '+colorama.Fore.YELLOW + f'{semicompressed_readslist[i].loc} '+colorama.Style.RESET_ALL + ' and '+colorama.Fore.YELLOW + f'{semicompressed_readslist[i+1].loc}.')
				if semicompressed_readslist[i].totalseqs>semicompressed_readslist[i+1].totalseqs:
					semicompressed_readslist[i+1].loc=semicompressed_readslist[i].loc #if read i has more counts than i+1, then change i+1's location to i
				semicompressed_readslist[i+1].totalseqs=semicompressed_readslist[i+1].totalseqs+semicompressed_readslist[i].totalseqs #if they are the same chromosome and locus, add to the 
			else:
				print(colorama.Style.RESET_ALL + f'adding line '+colorama.Fore.YELLOW + f'{i}'+colorama.Style.RESET_ALL + ' to compressed readslist for addition to csv',end="\r")
				compressed_readslist.append(semicompressed_readslist[i]) #if it's different than the next one, add it to the list (with the total count)
		return compressed_readslist
	else:
		return semicompressed_readslist

def read_sam(sam_file, chromIDS, ISbamfilename, compressreads=False,chromNTS={},randomize=False):
#if read_sam_file:
	message=[]
	samfile = pysam.AlignmentFile(sam_file, "r")
	ISbamfile = pysam.AlignmentFile(ISbamfilename,'wb',template=samfile)
	print(colorama.Fore.GREEN + f'Reading sam file: '+colorama.Style.RESET_ALL+f' {sam_file}')
	message.append(f'Reading sam file: '+f' {sam_file}')
	if randomize:
		print(colorama.Fore.ORANGE + f'Randomizing sequence locations.'+colorama.Style.RESET_ALL)
		message.append(f'Randomizing sequence locations.')
	print(colorama.Fore.GREEN + f'Writing bam file: '+colorama.Style.RESET_ALL+f' {ISbamfilename}')
	message.append(f'Writing bam file: '+f' {ISbamfilename}')
	readslist=[]
	recordlist=[]
	unmapped=[]
	entries=0
	for entry in samfile:
		entries+=1
		if entry.reference_name and (entry.reference_name[3:] in chromIDS.keys()): #check if read has a refrence to a chromosome, otherwise ignore (i.e. unmapped reads)
			chrom=entry.reference_name[3:] #chromosome number of mapping
			if entry.is_reverse:
				sense="-"
				address=entry.get_reference_positions()[-1] #chromosome position of mapping
				q=entry.query_qualities#copy quality scores
				entry.cigarstring='1M'
				entry.pos=address
				entry.query_sequence=entry.query_sequence[-1:]
				entry.query_qualities = q[-1:]
			else:
				sense="+"
				address=entry.get_reference_positions()[1] #chromosome position of mapping
				q=entry.query_qualities#copy quality scores
				entry.cigarstring='1M'
				entry.pos=address
				entry.query_sequence=entry.query_sequence[:1]
				entry.query_qualities = q[:1]
			if randomize:
				entry.pos=address=random.randrange(int(chromNTS[str(entry.reference_name[3:])]))
			ISbamfile.write(entry)#add entry to ISbamfile with 1 nt sequence
			readslist.append(read_csv_line([chrom, sense, address,'',1,'','','','','','']))#,chromNTS=chromNTS, userandomIS=random))
		else:
			unmapped.append(entry)
	print(colorama.Fore.YELLOW+f'{entries}'+colorama.Style.RESET_ALL+f' reads in sam file. '+colorama.Fore.YELLOW+f'{len(readslist)}'+colorama.Style.RESET_ALL+f' were mapped to a chromosome, '+colorama.Fore.YELLOW+f'{len(unmapped)} '+colorama.Style.RESET_ALL+f'did not map to a chromosome.')
	message.append(f'{entries} reads in sam file. {len(readslist)} were mapped to a chromosome, {len(unmapped)} did not map to a chromosome.')
	readslist=sorted(readslist, key=attrgetter('chrom', 'loc')) #sort the list by chromosome# and location
	if compressreads:
		readslist=compress(readslist)
	samfile.close()
	ISbamfile.close()
	return readslist, unmapped, "\n".join(message) #ISbamfile
#if write_csv:
def write_csv(readslist, mapped_csv_file):
	message=[]
	print(colorama.Style.RESET_ALL + f'Writing reads to csv file: '+colorama.Fore.YELLOW + f'{mapped_csv_file}')
	message.append(f'Writing reads to csv file: '+f'{mapped_csv_file}')
	with open(mapped_csv_file, "w", newline="") as csv_file:
		csv_writer = csv.writer(csv_file, delimiter=',') #comma delimited csv file
		#write csv headder to match expected format
		csv_writer.writerow(['Chrom','Sense','Loc','Total # IS Sequences Found','Total # IS Found','Plasmid m995','Sample1','Sample2','Sample3','Sample4'])
		csv_writer.writerow([])
		csv_writer.writerow([])
		for i in readslist:
			csvline=[i.chrom,i.sense,i.loc,i.gene,i.totalseqs,i.total,i.plasmid,i.sample1,i.sample2,i.sample3,i.sample4]
			csv_writer.writerow(csvline)
	return "\n".join(message)
#Chrom,Sense,Loc,Gene,Total # IS Sequences Found,Total # IS Found,Plasmid m995,Sample1,Sample2,Sample3,Sample4
#if read_csv:
def read_csv(mapped_csv_file,random=False):
	message=[]
	print(colorama.Style.RESET_ALL + f'reading data from '+colorama.Style.GREEN + f'{mapped_csv_file}',end="\r")
	message.append(f'reading data from '+f'{mapped_csv_file}')
	reads=[]
	with open(mapped_csv_file) as csv_file: #open file
		csv_reader = csv.reader(csv_file, delimiter=',') #read lines
		line_count = 1
		for row in csv_reader:
			if line_count < 4:#first 4 lines are headders
				if line_count==1: print(colorama.Fore.YELLOW + f'Column names are {", ".join(row)}')
				line_count += 1
			else: #read data from each line
				thisrow=read_csv_line(row,userandomIS=random)
				reads.append(thisrow)
				recordlist.append(thisrow.record)
				print(colorama.Style.RESET_ALL + f'Sequence {line_count-3} '+colorama.Style.RESET_ALL + 'found at chrom '+colorama.Style.GREEN + f'{thisrow.chrom}'+colorama.Style.RESET_ALL + ' with '+colorama.Style.ORANGE + f'{thisrow.sense}'+colorama.Style.RESET_ALL + ' sense '+colorama.Style.YELLOW + f'{thisrow.loc} '+colorama.Style.RESET_ALL + 'with sequence '+colorama.Style.BLUE + f'{thisrow.sequence[:5]}...{thisrow.sequence[-5:]}',end="\r" )
				line_count += 1 
		print(f'Processed {line_count} lines.')
		message.append(f'Processed {line_count} lines.')
	return reads, message