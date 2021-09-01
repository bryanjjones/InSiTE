#!/btapps/miniconda3/bin/python3
import _mypath
import csv
import sys
import colorama
import os
import pysam
import mappedreads
import argparse
import time

def import_chromosomes():
    chromosome_ids = './reference_datasets/chromosomes.csv'
    with open(chromosome_ids) as csv_file:
        chrom_ids = {}
        chrom_nts = {}
        csv_reader = csv.reader(csv_file, delimiter=",")
        line_count = 0
        for row in csv_reader:
            name = row[1]
            id_name = row[2]
            nts = row[4]
            chrom_ids[name] = id_name
            chrom_nts[name] = nts
    return chrom_ids, chrom_nts


def collision_remover(file_list, chromIDS, collisions_file_name):

    # mapped_locations_sets = []
    bam_file_list = []
    file_binary_dict = {}
    print("make temporary bam files from sam files")
    #TODO sum IH counts for collisions, add HI field to give bitwise occourance in each file
    # no easy way to add metadata
    # if previous:
    #     if previous.pos == entry.pos:  # if this entry matches the previous one, add one to previous, and
    #         # ignore this one.
    #         count = previous.get_tag('IH') + entry.get_tag('IH')
    #         previous.set_tag('IH', count)
    #     else:  # if previous exists and is different than this one, write prevous to compressed and set this
    #         # entry to previous.
    #         semicompressed_readslist.append(previous)
    #         # compressed.write(previous)
    #         previous = entry
    binary_counter = 1
    for i in range(len(file_list)):  #samfiles
        temp_bam_file = os.path.splitext(file_list[i])[0] + "_TEMP.bam"
        bam_file_list.append(temp_bam_file)
        file_binary_dict[file_list[i]] = binary_counter
        binary_counter = binary_counter * 2
        sam_file = pysam.AlignmentFile(file_list[i], 'r')
        int_site_bam_file = pysam.AlignmentFile(temp_bam_file, 'wb', template=sam_file)
        print(f'Reading sam file: ' + colorama.Fore.YELLOW + f' {file_list[i]}' + colorama.Fore.GREEN + f' . . . done.'
              + colorama.Style.RESET_ALL)
        print(f'Writing bam file: ' + colorama.Fore.YELLOW + f' {temp_bam_file}' + colorama.Style.RESET_ALL, end=' ')
        reads_list = []
        unmapped = []
        # entries = 0
        # adjust entries in sam_file to be to single mapped nt and only mapped reads
        for entry in sam_file:
            # entries += 1
            if entry.reference_name and (entry.reference_name[
                                         3:] in chromIDS.keys()):  # check if read has a refrence to a chromosome,
                # otherwise ignore (i.e. unmapped reads)
                if entry.is_reverse:
                    address = entry.get_reference_positions()[-1]  # chromosome position of mapping
                    q = entry.query_qualities  # copy quality scores
                    entry.cigarstring = '1M'
                    entry.pos = address
                    entry.query_sequence = entry.query_sequence[-1:]
                    entry.query_qualities = q[-1:]

                else:
                    address = entry.get_reference_positions()[1]  # chromosome position of mapping
                    q = entry.query_qualities  # copy quality scores
                    entry.cigarstring = '1M'
                    entry.pos = address
                    entry.query_sequence = entry.query_sequence[:1]
                    entry.query_qualities = q[:1]
                int_site_bam_file.write(entry)  # add entry to int_site_bam_file with 1 nt sequence
            else:
                unmapped.append(entry)
        sam_file.close()
        int_site_bam_file.close()
        print(colorama.Fore.GREEN + f' . . .  done.' + colorama.Fore.RESET)

    # sort those temporary bamfile entries by entry.refrence_name and entry.pos, then find all locations from all bam
    # files
    temp_collisions_file_name = os.path.splitext(collisions_file_name)[0] + "_TEMP.bam"
    print(f'Writing temporary collisions file: {temp_collisions_file_name}')
    collisions_file = pysam.AlignmentFile(temp_collisions_file_name, "wb",
                                          template=pysam.AlignmentFile(file_list[0], 'r'))
    all_mapped_locations = set()  #
    #  collisions = set()
    collisions = dict()
    for temp_bam_file in bam_file_list:
        mappedreads.compress(temp_bam_file, adjacent=False)
        #print(f'compressed duplicate reads in {temp_bam_file}.')
        int_site_bam_file = pysam.AlignmentFile(temp_bam_file, 'rb')
        print(f'Adding collisions from {temp_bam_file}', end=" ")
        for entry in int_site_bam_file:
            if entry.is_reverse:
                orientation = "-"
            else:
                orientation = "+"
            # int_site_sam_file.write(entry)
            name = entry.reference_name+orientation+str(entry.pos)
            if name in all_mapped_locations:
                if name in collisions:
                   pass
                else:
                    collisions[name] = [0, 0]
                    entry.set_tag('IH', 0)
                    entry.set_tag('HI', 0) # bitwise strig of h
                    collisions_file.write(entry)
            else:
                all_mapped_locations.add(name)
                #     mapped_locations_sets[i].append(entry)
        # int_site_bam_file.close()
        # int_site_sam_file.close()
        print(colorama.Fore.GREEN + f' . . .  done.' + colorama.Fore.RESET)
        print(f'Removing {temp_bam_file}.', end=' ')
        os.remove(temp_bam_file)
        print(colorama.Fore.GREEN + f' . . .  done.' + colorama.Fore.RESET)
    #  write collisions file
    collisions_file.close()
    print('all collisions added.')
    # remove collisions
    for i in range(len(file_list)):  # samfiles

        sam_file = pysam.AlignmentFile(file_list[i], 'r')
        no_collision_sam_file = pysam.AlignmentFile(os.path.splitext(file_list[i])[0] + "_no_collisions.sam", "w",
                                                    template=sam_file)
        print(f'Finding collisions in {file_list[i]}.', end=' ')
        # int_site_bam_file = pysam.AlignmentFile(temp_bam_file, 'wb', template=sam_file)
        # print(f'Reading sam file: ' + colorama.Fore.YELLOW + f' {file_list[i]}' + colorama.Style.RESET_ALL)
        # print(f'Writing bam file: ' + colorama.Fore.YELLOW + f' {temp_bam_file}' + colorama.Style.RESET_ALL)
        # no_collision_reads_list = []
        # unmapped = []
        # entries = 0
        # adjust entries in sam_file to be to single mapped nt and only mapped reads
        for entry in sam_file:
            if entry.reference_name and (entry.reference_name[3:] in chromIDS.keys()):
                if entry.is_reverse:
                    address = entry.get_reference_positions()[-1]  # chromosome position of mapping
                    orientation = "-"
                else:
                    address = entry.get_reference_positions()[1]  # chromosome position of mapping
                    orientation = "+"
                # int_site_sam_file.write(entry)
                name = entry.reference_name + orientation + str(address)
                if name in collisions:
                    collisions[name][0] += 1
                    if file_binary_dict[file_list[i]] > collisions[name][1]:
                        collisions[name][1] += file_binary_dict[file_list[i]]
                elif name in all_mapped_locations:
                    no_collision_sam_file.write(entry)
                else:
                    print("we have a problem not in all_mapped_locations")
                    print(f'Name: {entry.query_name}\nChromosome: {entry.reference_name}\nIs Reverse: {entry.is_reverse}\n'
                          f'Position: {str(entry.pos)}')
                    no_collision_sam_file.close()
                    all_seqs_dump = pysam.AlignmentFile("./all_seq_dump.sam", "w", template=sam_file)
                    sam_file.close()
                    for j in all_mapped_locations:
                        all_seqs_dump.write(j)
                    all_seqs_dump.close()
                    exit()
        print(colorama.Fore.GREEN + f' . . .  done.' + colorama.Fore.RESET)
        sam_file.close()
        no_collision_sam_file.close()
    print('Summing up all occurrences of each collision', end=' ')
    temp_collisions = pysam.AlignmentFile(temp_collisions_file_name, 'rb')
    collisions_file = pysam.AlignmentFile(collisions_file_name, 'w', template=temp_collisions)
    for collision in temp_collisions:
        if collision.is_reverse:
            orientation = "-"
        else:
            orientation = "+"
        name = collision.reference_name+orientation+str(collision.pos)
        assert name in collisions
        collision.set_tag('IH', collisions[name][0], "i")
        collision.set_tag('HI', str(collisions[name][1]), "Z")
        collisions_file.write(collision)
    temp_collisions.close()
    os.remove(temp_collisions_file_name)
    print('... Done.')
    print('Binary Index:')
    for file_name in file_list:
        print(file_name, end='\t')
        print(file_binary_dict[file_name])


if __name__ == "__main__":
    ap = argparse.ArgumentParser(prog=os.path.basename(sys.argv[0]),
                                 usage=__doc__)
    ap.add_argument('-f', '--files', metavar='filename.xxx', type=str, nargs="+", help='input list of files to check for collisions')
    # parser.add_argument('integers', metavar='N', type=int, nargs='+',
    #                     help='an integer for the accumulator')
    # ap.add_argument('-r', '--rand_is', default=False, action='store_true',
    #                 help='use random integration sites matched to given query set instead of actual query set')
    args = ap.parse_args()
    files = args.files
    print(files)
    programstart = time.time()
    chromIDS, a = import_chromosomes()
    collision_remover(files, chromIDS, 'collisions.sam')
    programstop = time.time()
    print(f'Completed all tasks in {str(int(programstop-programstart))} seconds.')
