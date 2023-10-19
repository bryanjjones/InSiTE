[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/bryanjjones/InSiTE/blob/master/Site_Integration_Analysis.ipynb)
# InSiTE
## Integration Site Analysis of Transposable Elements
  
- InSiTE is designed to analyze raw sequencing results of integration sites from viral or transposon systems.  
- InSiTE will provide useful analysis of integration trends and biases that may affect safety profiles, including count integration sites in annotated regions (e.g. intron/exon/transcript) and distances to specified features (e.g. transcriptional start site).  
  
## Explanation of usage  
### Example of typical usage:   
`python3 ./scripts/InSiTE.py -q ./raw_sequences.fastq -z --primer5 'ACTGACTG' --primer3 'GTCAGTCA'`  
  
### Input file types:   
- .fastq or .fasta (Also .csv or .sam/.bam, but these may still have bugs). Specified by `-q`/`-a`/`-c`/`-s`.  
- Sequence files can be provided with barcodes or primer sequences (`--barcode`, `--primer5`, `--primer3`). Paired reads can also be processed (`-p`).  

### Multiple mappings:  
Reads can be compressed with `-z`. This will merge reads mapping to the same position in the genome (+/- 1 nt) and count the number of reads. The counts and % of the most frequent and top 10 most frequent integration sites is reported in the log file and a .bam file is provided sorted by read counts. Read counts are stored in the .bam file under the 'IH:i:" flag. If reads are compressed, all stats regarding mapping to genome features weight all identical IS reads (+/-1nt) as a single integration. If reads are not compressed, all identical reads are counted when calculating integration into/near genome features--this may cause biased results if duplicate integration site reads are due to clonal out-growth or PCR bias and not due to independant integration events.  
  
### Random control data sets:  
As controls, randomized data can also be generated in two ways:   
- `-n` will generate random nucleotides to replace reads provided (number of nt read and quality scores will match input files). This is useful for calibrating non-specific mappings and adjusting --min to adjust number of reads and number of likely false alignments.  
- `-r` will generate random integration sites. The number of sites and sites per chromosome will be matched to the provided real dataset. The sites will be randomly distributed throughout the chromosome. This is useful for comparing actual data set and integration in/near genome features to a truely random distribution so that biases can be calculated.   
  
### Outputs:  
A number of output files are generated, as well as useful information printed in the terminal (& log file).  
- `[root_file]_trimmed.fast[a/q]` sequences after trimming off adapters/barcodes  
- `[root_file]_retrieved_2bit.fasta` fasta file of integration sites with `--rwindow` and `--lwindow` up and down-stream of integration sites, used for creating logo plot.  
- `[root_file]IS_logo.svg` image showing consensus logo of integration site  
- `[root_file].sam` sam file of mapped reads after trimming. mapped to genome location where possible.   
- `[rootfile]IS.bam`  bam file containing mapped single nt insertion sites. Sorted by position in genome (Chromosome & position).  
- `[root_file]_abundantsort.bam`  bam file containing single nt insertion sites, Sorted by frequency of reads to that site (count indicated by IH:i:## tag).  
- `[root_file]_IS_mappings.csv` file listing integration locations matching the format provided by GeneWerk  
- `[root_file]_distances.csv` file containing lists of distances for each read to features (e.g. TSS) for which distance mapping is indicated.  
- `[root_file]_IS_annotations.csv` summary file showing where reads mapped to for each indicated genome feature (exon, intron, TSS, etc)  
- `[root_file].log` log file of run  

## Requirements:  
### Directory structure:   
parent directory for contains '/scripts', '/refrence_datasets', and /modules folders. Main script is /scripts/InSiTE.py supporting modules are in /modules and the list of cromosomes names is in /refrence_datasets.  
  
### Manually downloaded local copies of genomes and genome annotation files:  
Annotation bed files for intron/exon/etc (these can be aquired from places like https://genome.ucsc.edu/cgi-bin/hgTables)  
.2bit genome file (available at https://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.2bit) (optional but increses performance, if not present specify `--remote`)  
bowtie2 genome index file (e.g. ftp://ftp.ncbi.nlm.nih.gov/genomes/archive/old_genbank/Eukaryotes/vertebrates_mammals/Homo_sapiens/GRCh38/seqs_for_alignment_pipelines/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index.tar.gz)  
  
### Dependencies:  
- bowtie2
- bbtools
- weblogo  
- TwoBitToFa  
- Python3 with packages:  
  - Biopython  
  - colorama  
  - pysam  
  - pybedtools  
  - twobitreader  
  - if python3 is at /btapps/miniconda3/bin/python3, script can be run directly, otherwise it needs to be run via python3.  
## Detailed Usage  
```
usage: InSiTE.py [-h] [-q FASTQ] [-r] [-n] [-c CSV] [-a FASTA] [-s SAM]
                 [--no_seqs] [-u] [-p PAIRS] [--no_annotate] [--barcode NNNNN]
                 [--vectors VECTORS] [--lwindow LWINDOW] [--rwindow RWINDOW]
                 [--samwindow SAMWINDOW] [--remote] [--min MIN]
                 [--primer5 NNNNN] [--primer3 NNNNN] [--trim5 TRIM5]
                 [--trim3 TRIM3] [--feature intron/exon/transcript/TSS/etc]
                 [--dist True/False] [--close CLOSE]
                 [--chromosome_ids /path/to/chromosomes.csv]
                 [--bowtielocation /path/to/bowtie2]
                 [--bowtieindex /path/to/bowtieindex]
                 [--weblogolocation /path/to/weblogo]
                 [--twobitlocation /path/to/TwoBitToFa]
                 [--twobitgenomelocation /path/to/genome.2bit]
                 [--annotations /path/to/annotation_file.bed/gff/gtf]
                 [--supress_csv] [--supress_fasta] [--supress_logo]
                 [--append_summary] [-v]

optional arguments:
  -h, --help            show this help message and exit
  -q FASTQ, --fastq FASTQ
                        fastq input file
  -r, --rand_is         use random integration sites matched to given query
                        set instead of actual query set
  -n, --rand_nt         use random sequences for mapping
  -c CSV, --csv CSV     csv input file. WARNING: This feature may not work.
  -a FASTA, --fasta FASTA
                        fasta input file
  -s SAM, --sam SAM     sam/bam input file
  --no_seqs             do not get sequences from either entrez or local
                        TwoBit genome around locations indicated by
                        genome_location_csv
  -u, --uncompress_reads
                        compress duplicate reads and reads shifted +/- 1nt,
                        number of reads are compressed in csv, fasta, and
                        mapping outputs
  -p PAIRS, --pairs PAIRS
                        specify file with paired reads for paired end reads
                        (used in conjunction with '-q' or '-a')
  --no_annotate         do not map insertion sites to genome annotations
  --barcode NNNNN       barcode sequence to trim off of reads
  --vectors VECTORS     FASTA file containing sequences te exclude from
                        mapping, for example plasmids used in the experiment
  --lwindow LWINDOW     numebr of nucleotides upstream of integration site to
                        return
  --rwindow RWINDOW     number of nucleotides downstream of integration site
                        to return
  --samwindow SAMWINDOW
                        depreciated
  --remote              get sequences from entrez server instead of 'LOCAL'
                        TwoBit genome
  --min MIN             minimum length of a read to try mapping. default (25)
                        will usually avoid any false positives in read sets of
                        200k reads
  --primer5 NNNNN       5' primer sequence to remove from reads
  --primer3 NNNNN       3' primer sequence to remove from reads
  --trim5 TRIM5         additional (non-genomic) nts to trim off of 3' end of
                        reads
  --trim3 TRIM3         additional (non-genomic) nts to trim off of 5' end of
                        reads
  --feature intron/exon/transcript/TSS/etc
                        feature names found in feature files to map reads to,
                        e.g. "exon"
  --dist True/False     weather to map distance of each read, or only whether
                        reads overlap with feature. Same number of distance
                        variables must be given as features.
  --close CLOSE         distance in bp to be considered close to feature
  --chromosome_ids /path/to/chromosomes.csv
  --bowtielocation /path/to/bowtie2
  --bowtieindex /path/to/bowtieindex
  --weblogolocation /path/to/weblogo
  --twobitlocation /path/to/TwoBitToFa
  --twobitgenomelocation /path/to/genome.2bit
  --annotations /path/to/annotation_file.bed/gff/gtf
                        location of annotation file(s) (bed/gff/gtf), must be
                        same number of files as features specified
  --supress_csv         do not output csv file
  --supress_fasta       do not ouptput fasta file
  --supress_logo        do not output logo
  --append_summary      add summary metrics to summary.csv file
  -v, --verbose         verbose output and logging
  --supress_fasta       do not ouptput fasta file
  --supress_logo        do not output logo
  --append_summary      add summary metrics to summary.csv file
  ```
