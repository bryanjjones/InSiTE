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
import pysam
import runbin
import mappedreads
import getseq
import annotate
import twobitreader
import argparse
import logging

#TODO check each logging event and move many to verbose only
#TODO check and fix or remove csv input type


# inputs:

userandomIS = 0  # 0 #useful for generating a random control data set chromosome matched to query set.
userandomSEQS = 0  # generate random sequence data, usefull for calculating false positive rate when mapping reads w/
# bowtie

# tasks:
getseqs = 1  # 1#get sequences from either entrez or local TwoBit genome around locations indicated by
# genome_location_csv
compressreads = 1  # remove duplicate reads and reads shifted +/- 1 nt, number of reads are compressed in csv and fasta
getannotations = 1  # 1 # use annotate.py module to map insertion sites to anotations

# parameters
barcode = ''  # 5' barcode sequence 05:ATCTGCGACG
lwindow = 50  # window on either side of indicated nt location to return
rwindow = 50
samwindow = 0  # 60
sequencesource = 'LOCAL'  # get sequences from 'LOCAL' TwoBit genome or 'REMOTE' entrez server
minimum_read_len = 25  # minimum length of a read to try mapping. 25 was found to usually avoid any false positives
# in read sets of 200k reads
primer5 = ''  # 5' primer sequenc to remove from reads
primer3 = ''  # 3' primer sequence to remove from reads
trim5 = 0  # additional (non-genomic) nts to trim off of 3' end of reads #21
trim3 = 0  # additional (non-genomic) nts to trim off of 5' end of reads (starts with TA) #16
featurenames = ['intron', 'exon', 'coding exon', 'transcript',
                'Transcriptional_Start_Site']  # feature names found in feature files to map reads to
featuredist = [False, False, False, False,
               True]  # weather to map distance of each read, or only whether reads overlap with feature
distance = [1000, 2000, 4000, 8000, 16000, 32000, 64000]  # distance in bp to be considered close to feature

# outputs:
write_csv = 1  # write csv (only relevant if reading from sam file)
writeFASTA = 1  # 1#write a fasta file of all sequences around IS
writelogo = 1  # 1#create a logo image of consensus sequence # requires getseqs and writeFASTA to be on


# reference file locations
chromosome_ids = './reference_datasets/chromosomes.csv'
bbmerge_location = './bins/bbmap/bbmerge-auto.sh'
bowtie_location = '/btapps/miniconda3/bin/bowtie2' #'bowtie2'
bowtie_index_file = './reference_datasets/genomes/GRCh38.fna.bowtie_index/GCA_000001405.15_GRCh38_no_alt_analysis_set.' \
              'fna.bowtie_index'
weblogo_location = 'weblogo'  # in PATH
twobit_location = './scripts/TwoBitToFa'
twobit_genome_location = './reference_datasets/genomes/GRCh38.2bit'
bbduk_location = './bins/bbmap/bbduk.sh'
annotations = ['./reference_datasets/annotations/refseq.introns.bed',
               './reference_datasets/annotations/refseq.exons.bed',
               './reference_datasets/annotations/refseq.codingexons.bed',
               './reference_datasets/annotations/refseq.transcripts.bed',
               './reference_datasets/annotations/refseq.transcripts.bed']  # ['./reerence_datasets/annotations/gencode
# .v32.introns.bed','./reference_datasets/annotations/gencode.v32.exons.bed',
# './reference_datasets/annotations/gencode.v32.codingexons.bed',
# './reference_datasets/annotations/gencode.v32.transcripts.gtf',
# './reference_datasets/annotations/gencode.v32.transcripts.gtf'] #[
# './reference_datasets/annotations/refseq.introns.bed','./reference_datasets/annotations/refseq.exons.bed',
# './reference_datasets/annotations/refseq.codingexons.bed','./reference_datasets/annotations/refseq.transcripts.bed',
# './reference_datasets/annotations/refseq.transcripts.bed'] annotations file in gff, gtf, or bed format,

verbose = False
pairedfile = None
bowtiecommand = None
readslist = None
append_summary = False
vector_tallies = []

#def process_fastq():

if __name__ == "__main__":

    ap = argparse.ArgumentParser(prog=os.path.basename(sys.argv[0]),
                                 usage=__doc__)
    ap.add_argument('-q', '--fastq', help='fastq input file') #TODO default=None
    ap.add_argument('-r', '--rand_is', default=False, action='store_true',
                    help='use random integration sites matched to given query set instead of actual query set')
    ap.add_argument('-n', '--rand_nt', default=False, action='store_true',
                    help='use random sequences for mapping')
    ap.add_argument('-c', '--csv', help='csv input file. WARNING: This feature may not work.')
    ap.add_argument('-a', '--fasta', help='fasta input file')
    ap.add_argument('-s', '--sam', help='sam/bam input file')
    ap.add_argument('--no_seqs', default=False, action='store_true',
                    help='do not get sequences from either entrez or local TwoBit genome around locations indicated '
                         'by genome_location_csv')
    ap.add_argument('-u', '--uncompress_reads', action='store_true', default=False,
                    help='compress duplicate reads and reads shifted +/- 1nt, number of reads are compressed in csv, '
                         'fasta, and mapping outputs')
    ap.add_argument('-p', '--pairs',
                    help="specify file with paired reads for paired end reads (used in conjunction with '-q' or '-a')") #TODO rever to default=None
    ap.add_argument('--no_annotate', default=False, action='store_true',
                    help='do not map insertion sites to genome annotations')
    ap.add_argument('--barcode', default='', metavar='NNNNN', help=' barcode sequence to trim off of reads')
    ap.add_argument('--vectors', help='FASTA file containing sequences te exclude from mapping, for '
                                                    'example plasmids used in the experiment') #TODO revert to None
    ap.add_argument('--lwindow', default=50, help='numebr of nucleotides upstream of integration site to return',
                    type=int)  # window on either side of indicated nt location to return
    ap.add_argument('--rwindow', default=50, help='number of nucleotides downstream of integration site to return',
                    type=int)  # window on either side of indicated nt location to return
    ap.add_argument('--samwindow', default=0, type=int, help='depreciated')  # 60
    ap.add_argument('--remote', action='store_true',
                    help="get sequences from entrez server instead of 'LOCAL' TwoBit genome")
    ap.add_argument('--min', default=25, type=int,
                    help='minimum length of a read to try mapping. default (25) will usually avoid any false positives '
                         'in read sets of 200k reads')
    #  ap.add_argument('--RC', default = None, metavar='file1 file2', help='XXX files to find & remove collisions, then reprocess')
    ap.add_argument('--primer5', default='', metavar='NNNNN', help="5' primer sequence to remove from reads")
    ap.add_argument('--primer3', default='', metavar='NNNNN', help="3' primer sequence to remove from reads")
    ap.add_argument('--trim5', default=0, type=int, help="additional (non-genomic) nts to trim off of 3' end of reads")
    ap.add_argument('--trim3', default=0, type=int,
                    help="additional (non-genomic) nts to trim off of 5' end of reads")  # (starts with TA)")
    ap.add_argument('--feature', action='append', metavar='intron/exon/transcript/TSS/etc',
                    help='feature names found in feature files to map reads to, e.g. "exon"')  # ['intron', 'exon',
    # 'codingexon','transcript','TSS'] #feature names found in feature files to map reads to
    ap.add_argument('--dist', action='append', metavar='True/False',
                    help='weather to map distance of each read, or only whether reads overlap with feature. '
                         'Same number of distance variables must be given as features.')
    ap.add_argument('--close', action='append', type=int, help='distance in bp to be considered close to feature')
    ap.add_argument('--chromosome_ids', metavar='/path/to/chromosomes.csv',
                    default='./reference_datasets/chromosomes.csv')
    ap.add_argument('--bbmergelocation', metavar='path/to/bbmerge', default='./bins/bbmap/bbmerge-auto.sh')
    ap.add_argument('--bbduklocation', metavar='path/to/bbduk', default='./bins/bbmap/bbduk.sh')
    ap.add_argument('--bowtielocation', metavar='/path/to/bowtie2', default='/btapps/miniconda3/bin/bowtie2')
    ap.add_argument('--bowtieindex', metavar='/path/to/bowtieindex',
                    default='./reference_datasets/genomes/GRCh38.fna.bowtie_index/'
                            'GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index')
    ap.add_argument('--bowtiethreads', default='4')
    ap.add_argument('--weblogolocation', metavar='/path/to/weblogo', default='weblogo')  # in PATH
    ap.add_argument('--twobitlocation', metavar='/path/to/TwoBitToFa', default='./scripts/TwoBitToFa')
    ap.add_argument('--twobitgenomelocation', metavar='/path/to/genome.2bit',
                    default='./reference_datasets/genomes/GRCh38.2bit')
    ap.add_argument('--annotations', action='append', metavar='/path/to/annotation_file.bed/gff/gtf',
                    help='location of annotation file(s) (bed/gff/gtf), must be same number of files as features '
                         'specified')  # ['./reference_datasets/annotations/gencode.v32.introns.bed',
    # './reference_datasets/annotations/gencode.v32.exons.bed',
    # './reference_datasets/annotations/gencode.v32.codingexons.bed',
    # './reference_datasets/annotations/gencode.v32.transcripts.gtf',
    # './reference_datasets/annotations/gencode.v32.transcripts.gtf']
    ap.add_argument('--supress_csv', default=False, action='store_true', help='do not output csv file')
    ap.add_argument('--supress_fasta', default=False, action='store_true', help='do not ouptput fasta file')
    ap.add_argument('--supress_logo', default=False, action='store_true', help='do not output logo')
    ap.add_argument('--append_summary', default=False, action='store_true',
                    help='add summary metrics to summary.csv file')
    ap.add_argument('-v', '--verbose', default=False, action='store_true', help='verbose output and logging')
    args = ap.parse_args()
    if args.fastq:
        inputfile = args.fastq
        inputtype = 'fastq'
    if args.csv:
        inputfile = args.csv
        inputtype = 'csv'
    if args.fasta:
        inputfile = args.fasta
        inputtype = 'fasta'
    if args.sam:
        inputfile = args.sam
        inputtype = 'sam'
    if not bool(args.fasta) + bool(args.fastq) + bool(args.csv) + bool(args.sam) == 1: #+ bool(args.RC)
        logging.critical('please provide a single input file (fasta/fastq/sam/bam/csv)')
        sys.exit('please provide a single input file (fasta/fastq/sam/bam/csv)')
    vector_file = args.vectors
    userandomIS = args.rand_is
    userandomSEQS = args.rand_nt
    getseqs = not args.no_seqs
    compressreads = not args.uncompress_reads  # remove duplicate reads and reads shifted +/- 1 nt, number of reads are
    # compressed in csv and fasta
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
        if not len(args.feature) == len(args.dist) == len(args.annotations):
            logging.critical(
                'must provide same number of features as distance and annotation files as they '
                'correspond to each other')
            sys.exit(
                'must provide same number of features as distance and annotation files as they '
                'correspond to each other')
        featurenames = args.feature
        featuredist = args.dist
        annotations = args.annotations
    if args.close:
        distance = args.close
    chromosome_ids = args.chromosome_ids
    bowtie_location = args.bowtielocation
    bowtie_index_file = args.bowtieindex
    bbmerge_location = args.bbmergelocation
    bbduk_location = args.bbduklocation
    weblogo_location = args.weblogolocation
    twobit_location = args.twobitlocation
    twobit_genome_location = args.twobitgenomelocation
    write_csv = not args.supress_csv
    writeFASTA = not args.supress_fasta
    writelogo = not args.supress_logo
    if args.pairs:
        pairedfile = args.pairs
    verbose = args.verbose
    append_summary = args.append_summary
    # if args.RC:
    #     remove_collisions_list = args.RC
    # specified file names
    if '.gz' in inputfile:
      rootname = inputfile.split('.fastq')[0]
    else:
      rootname = os.path.splitext(inputfile)[0]

    # output file names
    sam_file = f'{rootname}.sam'

    genome_location_csv = f'{rootname}_IS_mappings.csv'  # 'IS_data.csv'
    trimmedfastq = f'{rootname}_trimmed.fastq'  # fastq file with adapters/primers/barcodes trimmed off with cutadapt and
    # short sequences removed
    trimmedfasta = f'{rootname}_trimmed.fasta'  # fasta file with adapters/primers/barcodes trimmed off with cutadapt and
    vector_free_seqs = f'{rootname}_vector_free.fasta'  #fasta file with vector sequences removed
    # short sequences removed
    if '.gz' in inputfile:
      merged_reads = f'{rootname}_merged_reads.fastq'
    else:
      merged_reads = f'{rootname}_merged_reads{os.path.splitext(inputfile)[1]}'
    trimmedfastqpaired = f'{rootname}_trimmed_pairs.fastq'  # fastq file with adapters/primers/barcodes trimmed off with
    # cutadapt and short sequences removed
    trimmedfastapaired = f'{rootname}_trimmed_pairs.fasta'  # fastq file with adapters/primers/barcodes trimmed off with
    # cutadapt and short sequences removed
    FASTAfile = f'{rootname}_retrieved_2bit.fasta'  # output file containing sequences surrounding mapped insertion site
    logofile = f"{rootname}IS_logo.svg"  # logo file showing consensus integration site in logo format
    annotationsfile = f'{rootname}_IS_annotations.csv'  # file contoining summary of mapping to annotations
    distancesfile = f'{rootname}_distances.csv'  # file containing list of distances of each read to nearest TSS
    logfile = f'{rootname}.log'
    ISbamfilename = f'{rootname}IS.bam'  # sam file name for single nt IS mappings
    abundantfilename = f'{rootname}_abundantsort.sam'

    warnings = []
    # sanity checks:
    # TODO check if given files exist (sequence files and vector sequence file)
    if not os.path.isfile(inputfile):
        warnings.append(f"{inputfile} does not exist")
        sys.exit(
            colorama.Fore.RED + f"{inputfile} does not exist")
    if vector_file:
        if not os.path.isfile(vector_file):
            warnings.append(f"{vector_file} does not exist")
            sys.exit(colorama.Fore.RED + f"{vector_file} does not exist")
    if inputtype == "sam":
        read_sam_file = 1
        mapreads = 0
    elif inputtype == "csv":
        read_csv = 1
        mapreads = 0
        if getseqs == 0:
            warnings.append('Given a csv file, but not asked to get seqs. Nothing to do.')
            sys.exit(colorama.Fore.RED + 'Given a csv file, but not asked to get seqs. Nothing to do.')
    elif inputtype == "fasta":
        read_fasta = 1
        mapreads = 1
    elif inputtype == "fastq":
        read_fastq = 1
        mapreads = 1
        vector_free_seqs = f'{rootname}_vector_free.fastq'  # change vector_free_seqs file to fastq to maintain quality scores
    else:
        warnings.append("Specify input file type as 'sam', 'csv', 'fasta', or 'fastq'")
        sys.exit(colorama.Fore.RED + 'Specify input file type as \'sam\', \'csv\', \'fasta\', or \'fastq\'')
    if (writeFASTA == 1 or writelogo == 1) and getseqs == 0:
        warnings.append("Can't write a fasta file or logo file, or get annotations without getting sequences.")
        sys.exit(
            colorama.Fore.RED + "Can't write a fasta file or logo file, or get annotations without getting sequences.")
    if sequencesource != 'LOCAL' and sequencesource != 'REMOTE' and getseqs == 1:
        warnings.append("Must specify 'LOCAL' or 'REMOTE' source to get sequences")
        sys.exit(colorama.Fore.RED + "Must specify 'LOCAL' or 'REMOTE' source to get sequences")
    if getannotations and not len(annotations) == len(featurenames) == len(featuredist):
        warnings.append(
            "If getting annotations, length of annotation files, feature names, and feature distance must be the same")
        sys.exit(
            colorama.Fore.RED + "If getting annotations, length of annotation files, feature names, and feature distance "
                                "must be the same")
    if verbose:
        logging.basicConfig(filename=logfile, filemode='w', level=logging.DEBUG)
    else:
        logging.basicConfig(format='%(message)s', filename=logfile, filemode='w', level=logging.INFO)
    summary = {"file": rootname, "barcode sequence": barcode}
        #[rootname, barcode]
    reads = []
    recordlist = []
    chromIDS = {}
    chromNTS = {}
    programstart = time.time()
    if len(warnings) > 0:
        logging.critical("\n".join(warnings))
    # make a dictionary of chromosome names and ID #s using chromosomes.csv, which is the table copied from
    # https://www.ncbi.nlm.nih.gov/genome?term=human&cmd=DetailsSearch
    with open(chromosome_ids) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=",")
        line_count = 0
        for row in csv_reader:
            name = row[1]
            ID = row[2]
            nts = row[4]
            chromIDS[name] = ID
            chromNTS[name] = nts

    if userandomIS or userandomSEQS:
        random.seed(a=None, version=2)  # initialized random number generator w/ no seed (use system time)

    if sequencesource == 'LOCAL':
        genome = twobitreader.TwoBitFile(twobit_genome_location)  # reads genome into object
    else:
        genome = None

    #  if paired reads, merge what can be merged, and use the rest as individual reads. If insertion site is in a gap,
    #  it is not useful, thus gappd pairs are more useful if treated as unpaired reads.
    if pairedfile: # and not os.path.isfile(merged_reads):
        print(f'Merging paired reads files {inputfile} and {pairedfile} using bbmerge.')
        FASTQ.merge_pairs(inputfile, pairedfile, merged_reads, bbmerge_location=bbmerge_location)
        inputfile=merged_reads


    # trim primers and adapters and filter for length from raw fastqreads, return trimmed "genomic" sequences in fasta
    # format.
    if inputtype == "fastq":
        print(f'FASTQ format input file: ' + colorama.Fore.YELLOW + f'{inputfile}' + colorama.Style.RESET_ALL)
        trimlog = FASTQ.Trim(inputfile, trimmedfastq, barcode, primer5, primer3, trim3, trim5, minimum_read_len)
        logging.debug(trimlog)
        if '.gz' in inputfile:
          summary["raw reads"] = int(trimlog.split('\n')[3].split('\t')[1])
          summary["trimmed reads"] = int(trimlog.split('\n')[3].split('\t')[6])
        else:
          summary["raw reads"] = (int(sum(1 for line in open(inputfile)) / 4))  # 4 lines per entry in fastq
          summary["trimmed reads"] = (int(sum(1 for line in open(trimmedfastq)) / 4))

        if summary["trimmed reads"] == 0:
            logging.critical(f'No sequences left after trimming. Exiting')
            sys.exit(colorama.Fore.RED + f'No sequences left after trimming. Exiting' + colorama.Style.RESET_ALL)
        if userandomSEQS:  # replace all real reads with random NT
            trimlog = FASTQ.randomize(trimmedfastq, format='fastq')
            logging.debug(trimlog)
        # run bowtie using trimmed fastq and quality scores (--phred33)
        if vector_file:
            # if pairedfile:
            #     logging.critical(f'Can not weed-out vector mappings while using paired end reads. Exiting')
            #     sys.exit(colorama.Fore.RED + f'Can not weed-out vector mappings while using paired end reads. Exiting'
            #              + colorama.Style.RESET_ALL)
            #vector_tallies, non_vector_reads = mappedreads.filter_vector(vector_file, trimmedfastq, filetype="fastq",
            #                                                             trimmed_file=trimmedfastq)
            bbduk_command = f'{bbduk_location} -Xmx1g in={trimmedfastq} out={vector_free_seqs} ref={vector_file} ' \
                            f'k={min(minimum_read_len,31)} hdist=1 stats=stats.txt findbestmatch=t nzo=f'
            runbin.Command(bbduk_command)
            logging.debug(f'Weeding out vector sequences using bbduk. Writing output to' + f' {vector_free_seqs}')
            logging.debug(f'{bbduk_command}')
            bbduk = runbin.Command(bbduk_command)
            bbduk_out = bbduk.run(timeout=20000)
            print(bbduk_out[1].decode())
            print(f'{bbduk_out[2].decode()}')
            with open('./stats.txt', 'r', newline='') as bbduk_stats:
                vector_tallies = []
                bbduk_stats_csv = csv.reader(bbduk_stats, delimiter="\t")
                matched = total = 0
                for ln in bbduk_stats_csv:
                    if ln[0][0] != "#":
                        vector_tallies.append([ln[0], ln[1]])
                    elif ln[0] == "#Total":
                        total = ln[1]
                    elif ln[0] == "#Matched":
                        matched = ln[1]
                non_vector_reads = int(total) - int(matched)
            os.remove('./stats.txt')
            for i in vector_tallies:
                message = f'{str(i[1])} reads mapped to vector {str(i[0])}'
                logging.info(message)
                print(colorama.Fore.YELLOW + message + colorama.Style.RESET_ALL)
            message = f'{str(non_vector_reads)} reads that did not map to any vector'
            summary["reads after filtering vector sequence"] = non_vector_reads
            trimmedfastq = vector_free_seqs
        # if pairedfile:
        #     print(f'using paired reads from ' + colorama.Fore.YELLOW + f'{inputfile}'+colorama.Style.RESET_ALL + ' and ' +
        #           colorama.Fore.YELLOW + f'{pairedfile}' + colorama.Style.RESET_ALL)
        #     bowtiecommand = f'{bowtie_location} --phred33 -p 4 -x {bowtie_index_file} -1 {trimmedfastq} -2 ' \
        #                     f'{trimmedfastqpaired} --fr --no-unal -S {sam_file}'
        # else:
        bowtiecommand = f'{bowtie_location} --phred33 -p {args.bowtiethreads} -x {bowtie_index_file} -U {trimmedfastq} --no-unal ' \
                            f'-S {sam_file}'
    # trim priers and adapters and filter for length from fasta reads, return trimmed "genomic" sequences in fasta format.
    elif inputtype == "fasta":
        print(f'FASTA format input file: ' + colorama.Fore.YELLOW + f'{inputfile}' + colorama.Style.RESET_ALL)
        trimlog = FASTQ.Trim(inputfile, trimmedfasta, barcode, primer5, primer3, trim3, trim5, minimum_read_len,
                             filetype='fasta')  #, paired=pairedfile, pairedoutputfile=trimmedfastapaired)
        logging.debug(trimlog)
        with open(inputfile, "r") as fi:
            sequencecounter = 0
            for ln in fi:
                if ln.startswith(">"):
                    sequencecounter += 1
            summary["raw reads"] = sequencecounter
        with open(trimmedfasta, "r") as fi:
            sequencecounter = 0
            for ln in fi:
                if ln.startswith(">"):
                    sequencecounter += 1
            summary["trimmed reads"] = sequencecounter
        if sum(1 for line in open(trimmedfasta)) == 0:
            logging.critical(f'No sequences left after trimming.')
            sys.exit(colorama.Fore.RED + f'No sequences left after trimming. Exiting' + colorama.Style.RESET_ALL)
        if userandomSEQS:  # replace all real reads with random NT
            trimlog = FASTQ.randomize(trimmedfastq, format='fasta')
            logging.debug(trimlog)
        if vector_file:
            # if pairedfile:
            #     logging.critical(f'Can not weed-out vector mappings while using paired end reads. Exiting')
            #     sys.exit(colorama.Fore.RED + f'Can not weed-out vector mappings while using paired end reads. Exiting'
            #              + colorama.Style.RESET_ALL)
            #     #TODO add bbduk function with paired files
            bbduk_command = f'{bbduk_location} -Xmx1g in={trimmedfasta} out={vector_free_seqs} ref={vector_file} k={min(minimum_read_len,31)} hdist=1 stats=stats.txt findbestmatch=t nzo=f'
            runbin.Command(bbduk_command)
            logging.debug(f'Weeding out vector sequences using bbduk. Writing output to' + f' {vector_free_seqs}')
            logging.debug(f'{bowtiecommand}')
            bbduk = runbin.Command(bbduk_command)
            bbduk_out = bbduk.run(timeout=20000)
            print(bbduk_out[1].decode())
            print(f'{bbduk_out[2].decode()}')
            with open('./stats.txt', 'r', newline='') as bbduk_stats:
                vector_tallies = []
                bbduk_stats_csv = csv.reader(bbduk_stats, delimiter="\t")
                matched = total = 0
                for ln in bbduk_stats_csv:
                    if ln[0][0] != "#":
                        vector_tallies.append([ln[0], ln[1]])
                    elif ln[0] == "#Total":
                        total = ln[1]
                    elif ln[0] == "#Matched":
                        matched = ln[1]
                non_vector_reads = int(total) - int(matched)
            os.remove('./stats.txt')
            #vector_tallies, non_vector_reads = mappedreads.filter_vector(vector_file, trimmedfasta, filetype="fasta", trimmed_file=trimmedfasta)
            for i in vector_tallies:
                message = f'{str(i[1])} reads mapped to vector {str(i[0])}'
                logging.info(message)
                print(colorama.Fore.YELLOW + message + colorama.Style.RESET_ALL)
            trimmedfasta = vector_free_seqs
            message = f'{str(non_vector_reads)} reads that did not map to any vector'
            summary["read after filtering vector sequence"] = non_vector_reads
        # run bowtie using trimmed fastq and quality scores (--phred33)
        # if pairedfile:
        #     print(f'using paired reads from ' + colorama.Fore.YELLOW + f'{inputfile}' + colorama.Style.RESET_ALL +
        #           ' and ' + colorama.Fore.YELLOW + f'{pairedfile}' + colorama.Style.RESET_ALL)
        #     bowtiecommand = f'{bowtie_location} --phred33 -p 4 -x {bowtie_index_file} -1 {trimmedfasta} ' \
        #                     f'-2 {trimmedfastapaired} --no-unal -S {sam_file}'
        # else:
        bowtiecommand = f'{bowtie_location} --phred33 -p 4 -f -x {bowtie_index_file} -U {trimmedfasta} --no-unal ' \
                            f'-S {sam_file}'

    # map reads from fasta file to genome, return sam file with genome locations
    elif mapreads:  # if mapreads, but not using fastq, run bowtie specifying fasta (-f)
        fastareads = f'{rootname}.fasta'
        if userandomSEQS:
            FASTQ.randomize(fastareads, format='fasta')
        # add bowtie mapping function to map plasmid sequences?
        bowtiecommand = f'{bowtie_location} -f -x {bowtie_index_file} -p 4 -U {fastareads} --no-unal -S {sam_file}'
    if mapreads:
        print(f'mapping reads genome using bowtie2. Writing output to' + colorama.Fore.YELLOW + f' {sam_file}' +
              colorama.Style.RESET_ALL + f' and ' + colorama.Fore.YELLOW + f'{rootname}_bowtie.log' +
              colorama.Style.RESET_ALL)
        print(colorama.Fore.CYAN + f'{bowtiecommand}' + colorama.Style.RESET_ALL)
        logging.debug(f'mapping reads genome using bowtie2. Writing output to' + f' {sam_file}' + f' and ' +
                      f'{rootname}_bowtie.log')
        logging.debug(f'{bowtiecommand}')
        bowtie = runbin.Command(bowtiecommand)
        bowtieout = bowtie.run(timeout=20000)
        print(bowtieout[1].decode())
        print(f'{bowtieout[2].decode()}')
        logging.debug(bowtieout[1].decode())
        logging.debug(bowtieout[2].decode())

    if inputtype == "sam" or mapreads:
        readslist, unmapped, message, topclone, topten = mappedreads.read_sam(sam_file, chromIDS, ISbamfilename, compressreads=compressreads,
                                                            chromNTS=chromNTS, randomize=userandomIS,
                                                            abundant=abundantfilename)
        # count total mapped sequences
        sam = pysam.AlignmentFile(sam_file, 'r')
        alignedcount = 0
        for i in sam:
            alignedcount += 1
        summary['aligned reads'] = alignedcount
        logging.info(message)
        summary['most abundant clone reads'] = topclone
        summary['Top 10 most abundant clone reads'] = topten
    if write_csv:
        message = mappedreads.write_csv(readslist, genome_location_csv)
        logging.debug(message)

    if getseqs:
        recordlist, readslist = getseq.get_seqs(readslist, lwindow, rwindow, samwindow, genome, source=sequencesource)

    if writeFASTA:  # write FASTA formatted file with returned sequences
        print(colorama.Style.RESET_ALL + f'\nWriting fasta file:' + colorama.Fore.YELLOW + f' {FASTAfile}' +
              colorama.Style.RESET_ALL)
        Bio.SeqIO.write(recordlist, FASTAfile, "fasta")
    if writelogo:  # dependant on having sequences, optional to make logo plot
        weblogocommand = f'{weblogo_location} -f {FASTAfile} -D fasta -o {logofile} ' \
                         f'-F svg -A dna -F png --resolution 600 -s large -c classic -i {str(int(1 - 1 * lwindow))} ' \
                         f'-l -10 -u 10 '
        print(colorama.Style.RESET_ALL + f'Writing logo')
        print(colorama.Fore.CYAN + f'{weblogocommand}' + colorama.Style.RESET_ALL)
        logging.debug(f'Writing logo')
        logging.debug(f'{weblogocommand}')
        logocommand = runbin.Command(weblogocommand)
        logoout = logocommand.run(timeout=1800)

    if getannotations and len:
        with open(annotationsfile, 'w', newline='') as output_file:
            totalprinted = False
            output_writer = csv.writer(output_file, delimiter=",")
            for i in range(len(featurenames)):
                if featuredist[i]:
                    print(f'mapping insertion site distances to {featurenames[i]} in ' + colorama.Fore.YELLOW +
                          f'{annotations[i]}' + colorama.Style.RESET_ALL)
                    logging.debug(
                        f'mapping insertion site distances to ' + f'{featurenames[i]}' + f' in ' + f'{annotations[i]}')
                    distances, average, standarddev, median, distancebins, b, message = annotate.closest(ISbamfilename,
                                                                                                         annotations[i],
                                                                                                         featurename=
                                                                                                         featurenames[i],
                                                                                                         distances=distance,
                                                                                                         position="start")
                    logging.info(message)
                    try:
                        output_writer.writerow([f'{featurenames[i]}(average distance)', average])
                        output_writer.writerow([f'{featurenames[i]}(standard deviation)', standarddev])
                        output_writer.writerow([f'{featurenames[i]}(median)', median])
                        summary['average distance'] = average
                        summary['standard deviation distance'] = standarddev
                        summary['median distance'] = median
                        for j in range(len(distance)):
                            output_writer.writerow(
                                [f'integration events within {distance[j]} bp of {featurenames[i]}', distancebins[j]])
                            summary[f'{distance[j]} bp of {featurenames[i]}'] = distancebins[j]
                    except:
                        logging.error("can't write stats. Maybe there were no sequences")
                        print(colorama.Fore.RED + "can't write stats. Maybe there were no sequences" +
                              colorama.Style.RESET_ALL)
                    with open(distancesfile, 'a+', newline='') as distance_file:
                        dist_writer = csv.writer(distance_file, delimiter=",")
                        dist_writer.writerow(distances)
                else:
                    print(f'mapping insertion sites to {featurenames[i]} in ' + colorama.Fore.YELLOW + f'{annotations[i]}' +
                          colorama.Style.RESET_ALL)
                    logging.debug(f'mapping insertion sites to' + f' {featurenames[i]}' + f' in ' + f'{annotations[i]}')
                    results, message = annotate.featuremap(annotations[i], ISbamfilename, featurenames=featurenames[i],
                                                           single=True, procs=3)
                    if not totalprinted:  # featuremap returns total reads, add this line only once as it should be the
                        # same for each feature.
                        output_writer.writerow([f'total mapped reads', results[1]])
                        summary['total mapped sequences'] = results[1]
                        totalprinted = True
                    output_writer.writerow([f'reads mapped to {featurenames[i]}', results[0]])
                    summary[featurenames[i]] = results[0]
                    logging.info(message)


    if append_summary:  # append summary, using headers at the top of existing summary or if not all categories to write
        # are listed in the header, write header specific to this line.
        print('Writing Summary')
        new_summary_header = False
        if os.path.exists('./summary.csv'):
            print("appending to existing summary file")
            with open('./summary.csv', newline='') as summary_file:
                reader = csv.reader(summary_file)
                summary_header = next(reader)
        else:
            summary_header = []
        summary_line = [""] * len(summary_header)
        for i in range(len(summary_header)):
            if summary_header[i] in summary:
                summary_line[i] = summary[summary_header[i]]
        for element in summary:
            if element in summary_header:
                summary_line[summary_header.index(element)] = summary[element]
            else:
                new_summary_header = True
                summary_header.append(element)
                summary_line.append(summary[element])
        with open('./summary.csv', 'a') as summary_file:
            if new_summary_header:
                csv.writer(summary_file, delimiter=',').writerow(summary_header)
            csv.writer(summary_file, delimiter=',').writerow(summary_line)
    programend = time.time()
    message = (colorama.Fore.GREEN + f'Completed all tasks for {rootname} in {int(programend - programstart)} seconds. '
                                     f'Exiting.' + colorama.Style.RESET_ALL)
    logging.debug(f'Completed all tasks for {rootname} in {int(programend - programstart)} seconds. Exiting.')
    print(message)
    exit()
