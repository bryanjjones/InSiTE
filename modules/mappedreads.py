#!/btapps/miniconda3/bin/python3
import csv
import pysam
import colorama
import Bio
import Bio.Entrez
import Bio.SeqIO
import random
from operator import attrgetter


# read take a given string, formatted as row read from a csv file, and return an object with useful properties
class read_csv_line(object):
    def __init__(self, row):  # , userandomIS=False, chromNTS={}):
        # tries = 1
        self.chrom = str(row[0])
        self.sense = row[1]
        if row[1] == "+":
            self.sensenum = 1
        elif row[1] == "-":
            self.sensenum = 2
        self.loc = int(row[2])
        self.gene = row[3]
        self.totalseqs = row[4]
        self.total = row[5]
        self.plasmid = row[6]
        self.sample1 = row[7]
        self.sample2 = row[8]
        self.sample3 = row[9]
        self.sample4 = row[10]
        self.record = None
        self.sequence = '-'


# if compressedreads:
# compress a sorted readlist by merging identical locations and, if adjacent=True, merging adjacent locations too
# (merges to location with greatest count)
def compress(bamfile, compressedbam=None, adjacent=True):  # compressedbam=f'compressed_{bamfile}',
    if not compressedbam:  # if not given a new file name, overwrite existing bam file
        compressedbam = bamfile
    print(f'compressing duplicate reads in {bamfile}. Writing compressed list with counts to {compressedbam}.')
    semicompressed_readslist = []
    compressed_readslist = []
    expanded = pysam.AlignmentFile(bamfile, 'rb')

    previous = None
    for entry in expanded:  # merge
        entry.set_tag('IH', 1)
        if previous:
            if previous.pos == entry.pos:  # if this entry matches the previous one, add one to previous, and
                # ignore this one.
                count = previous.get_tag('IH') + entry.get_tag('IH')
                previous.set_tag('IH', count)
            else:  # if previous exists and is different than this one, write prevous to compressed and set this
                # entry to previous.
                semicompressed_readslist.append(previous)
                # compressed.write(previous)
                previous = entry
        else:
            previous = entry
    semicompressed_readslist.append(previous)  # write the last read to compressed
    compressed = pysam.AlignmentFile(compressedbam, 'wb', template=expanded)
    expanded.close()
    if adjacent:
        compressedcount = 0
        for i in range(len(semicompressed_readslist)):
            if i + 1 < len(semicompressed_readslist) and semicompressed_readslist[i].pos == (
                    semicompressed_readslist[i + 1].pos - 1) and semicompressed_readslist[i].reference_name == \
                    semicompressed_readslist[
                        i + 1].reference_name:  # if this read isn't the last one and if it's the same as the next one
                print(
                    colorama.Style.RESET_ALL + f'merging adjacent read locations. read ' + colorama.Fore.YELLOW +
                    f'{i}' + colorama.Style.RESET_ALL + ', locations ' + colorama.Fore.YELLOW +
                    f'{semicompressed_readslist[i].pos} ' + colorama.Style.RESET_ALL + ' and ' +
                    colorama.Fore.YELLOW + f'{semicompressed_readslist[i + 1].pos}.',
                    end="\r")
                if semicompressed_readslist[i].get_tag('IH') > semicompressed_readslist[i + 1].get_tag('IH'):
                    semicompressed_readslist[i + 1].pos = semicompressed_readslist[
                        i].pos  # if read i has more counts than i+1, then change i+1's location to i
                summedreads = semicompressed_readslist[i + 1].get_tag('IH') + semicompressed_readslist[i].get_tag(
                    'IH')  # if they are the same chromosome and locus, add to the
                semicompressed_readslist[i + 1].set_tag('IH', summedreads)
            else:
                print(
                    colorama.Style.RESET_ALL + f'adding line ' + colorama.Fore.YELLOW + f'{i}' +
                    colorama.Style.RESET_ALL + f' to compressed {compressedbam}',
                    end="\r")
                compressedcount += 1
                compressed.write(semicompressed_readslist[i])
            # compressed_readslist.append(semicompressed_readslist[i]) #if it's different than the next one,
        # add it to the list (with the total count)
        print(f'{compressedcount} sequences written to {compressedbam}')
        compressed.close()
    # return compressed_readslist
    else:
        for i in semicompressed_readslist:
            compressed.write(i)
        compressed.close()
    # return semicompressed_readslist


def read_sam(sam_file, chromIDS, ISbamfilename, compressreads=True, chromNTS={}, randomize=False, expandbam=False,
             abundant=None):  # sorted_file=f'{sam_file}_sorted.bam',
    # if read_sam_file:
    message = []
    samfile = pysam.AlignmentFile(sam_file, 'r')
    ISbamfile = pysam.AlignmentFile(ISbamfilename, 'wb', template=samfile)
    print(colorama.Fore.GREEN + f'Reading sam file: ' + colorama.Style.RESET_ALL + f' {sam_file}')
    if randomize:
        print(colorama.Fore.MAGENTA + f'Randomizing sequence locations.' + colorama.Style.RESET_ALL)
        message.append(f'Randomizing sequence locations.')
    print(colorama.Fore.GREEN + f'Writing bam file: ' + colorama.Style.RESET_ALL + f' {ISbamfilename}')
    readslist = []
    recordlist = []
    unmapped = []
    entries = 0
    # adjust entries in samfile to be to singl mapped nt and only mapped reads

    for entry in samfile:
        entries += 1
        if entry.reference_name and (entry.reference_name[
                                     3:] in chromIDS.keys()):  # check if read has a refrence to a chromosome,
            # otherwise ignore (i.e. unmapped reads)
            if entry.is_reverse:
                sense = "-"
                address = entry.get_reference_positions()[-1]  # chromosome position of mapping
                q = entry.query_qualities  # copy quality scores
                entry.cigarstring = '1M'
                entry.pos = address
                entry.query_sequence = entry.query_sequence[-1:]
                entry.query_qualities = q[-1:]

            else:
                sense = "+"
                address = entry.get_reference_positions()[1]  # chromosome position of mapping
                q = entry.query_qualities  # copy quality scores
                entry.cigarstring = '1M'
                entry.pos = address
                entry.query_sequence = entry.query_sequence[:1]
                entry.query_qualities = q[:1]
            if randomize:
                entry.pos = random.randrange(int(chromNTS[str(entry.reference_name[3:])]))
            ISbamfile.write(entry)  # add entry to ISbamfile with 1 nt sequence
        else:
            unmapped.append(entry)
    samfile.close()
    ISbamfile.close()
    # sort those bamfile entries by entry.refrence_name and entry.pos
    pysam.sort("-o", ISbamfilename, ISbamfilename)

    # compress bamfile
    if compressreads:
        compress(ISbamfilename)
    # convert bamfile entries to csv
    ISbamfile = pysam.AlignmentFile(ISbamfilename, 'rb')
    abundantlist = []
    for entry in ISbamfile:
        if entry.is_reverse:
            sense = "-"
        else:
            sense = "+"
        chrom = entry.reference_name[3:]  # chromosome number of mapping
        address = entry.pos
        if compressreads:
            count = entry.get_tag('IH')
            if abundant:
                abundantlist.append((count, entry))
        else:
            count = 1
        readslist.append(read_csv_line(
            [chrom, sense, address, '', count, '', '', '', '', '', '']))
    print(colorama.Fore.YELLOW + f'{entries}' + colorama.Style.RESET_ALL + f' reads in sam file. '
          + colorama.Fore.YELLOW + f'{len(readslist)}' + colorama.Style.RESET_ALL + f' were mapped to a chromosome, '
          + colorama.Fore.YELLOW + f'{len(unmapped)} ' + colorama.Style.RESET_ALL + f'did not map to a chromosome.')
    message.append(f'{entries} reads in sam file. {entries - len(unmapped)} were mapped to a chromosome, '
                   f'{len(unmapped)} did not map to a chromosome.')
    if abundant and compressreads:  # if abundant filename given, sort reads by most abundant,
        # and write to that file as bam file
        abundantbamfile = pysam.AlignmentFile(abundant, 'wb', template=ISbamfile)
        abundantlist.sort(key=lambda tup: tup[0], reverse=True)
        mostover = round((100 * abundantlist[0][0] / entries), 3)
        print(
            colorama.Style.RESET_ALL + f'Most abundant clone found ' + colorama.Fore.YELLOW + f'{abundantlist[0][0]}' +
            colorama.Style.RESET_ALL + f' times, ' + colorama.Fore.YELLOW + f'{mostover}% ' + colorama.Style.RESET_ALL +
            f'of total reads.')
        message.append(f'Most abundant clone found {abundantlist[0][0]} times, {mostover}% of total reads.')
        topten = 0
        for i in abundantlist[0:10]:
            topten += i[0]
        toptenpct = round(100 * (topten / entries), 3)
        print(colorama.Style.RESET_ALL + f'Top ten most abundant clones found a total of ' + colorama.Fore.YELLOW +
              f'{topten}' + colorama.Style.RESET_ALL + f' times, ' + colorama.Fore.YELLOW + f'{toptenpct}% ' +
              colorama.Style.RESET_ALL + f'of total reads.')
        print(
            colorama.Fore.GREEN + f'Writing bam file sorted by read count: ' + colorama.Style.RESET_ALL + f'{abundant}')
        message.append(f'Top ten most abundant clones found a total of {topten} times, {toptenpct}% of total reads.')
        #message.append(f'Writing bam file sorted by read count: {abundant}')
        for i in abundantlist:
            abundantbamfile.write(i[1])
        abundantbamfile.close()
    elif abundant:
        print(colorama.Fore.RED + f'Need to count reads to calculate most abundant sequences. If you want abundant '
                                  f'sequences, please specify compressreads too.' + colorama.Style.RESET_ALL)

    ISbamfile.close()
    return readslist, unmapped, "\n".join(message)  # ISbamfile


# if write_csv:
def write_csv(readslist, mapped_csv_file):
    message = []
    print(colorama.Style.RESET_ALL + f'Writing reads to csv file: ' + colorama.Fore.YELLOW + f'{mapped_csv_file}')
    message.append(f'Writing reads to csv file: ' + f'{mapped_csv_file}')
    with open(mapped_csv_file, "w", newline="") as csv_file:
        csv_writer = csv.writer(csv_file, delimiter=',')  # comma delimited csv file
        # write csv headder to match expected format
        csv_writer.writerow(
            ['Chrom', 'Sense', 'Loc', 'Total # IS Sequences Found', 'Total # IS Found', 'Plasmid m995', 'Sample1',
             'Sample2', 'Sample3', 'Sample4'])
        csv_writer.writerow([])
        csv_writer.writerow([])
        for i in readslist:
            csvline = [i.chrom, i.sense, i.loc, i.gene, i.totalseqs, i.total, i.plasmid, i.sample1, i.sample2,
                       i.sample3, i.sample4]
            csv_writer.writerow(csvline)
    return "\n".join(message)


# Chrom,Sense,Loc,Gene,Total # IS Sequences Found,Total # IS Found,Plasmid m995,Sample1,Sample2,Sample3,Sample4
# if read_csv:
def read_csv(mapped_csv_file, random=False):
    message = []
    print(colorama.Style.RESET_ALL + f'reading data from ' + colorama.Style.GREEN + f'{mapped_csv_file}', end="\r")
    message.append(f'reading data from ' + f'{mapped_csv_file}')
    reads = []
    recordlist = []
    with open(mapped_csv_file) as csv_file:  # open file
        csv_reader = csv.reader(csv_file, delimiter=',')  # read lines
        line_count = 1
        for row in csv_reader:
            if line_count < 4:  # first 4 lines are headders
                if line_count == 1: print(colorama.Fore.YELLOW + f'Column names are {", ".join(row)}')
                line_count += 1
            else:  # read data from each line
                thisrow = read_csv_line(row, userandomIS=random)
                reads.append(thisrow)
                recordlist.append(thisrow.record)
                print(
                    colorama.Style.RESET_ALL + f'Sequence {line_count - 3} ' + colorama.Style.RESET_ALL +
                    'found at chrom ' + colorama.Style.GREEN + f'{thisrow.chrom}' + colorama.Style.RESET_ALL +
                    ' with ' + colorama.Style.ORANGE + f'{thisrow.sense}' + colorama.Style.RESET_ALL + ' sense ' +
                    colorama.Style.YELLOW + f'{thisrow.loc} ' + colorama.Style.RESET_ALL + 'with sequence ' +
                    colorama.Style.BLUE + f'{thisrow.sequence[:5]}...{thisrow.sequence[-5:]}',
                    end="\r")
                line_count += 1
        print(f'Processed {line_count} lines.')
        message.append(f'Processed {line_count} lines.')
    return reads, message
