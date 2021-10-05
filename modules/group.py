#!/btapps/miniconda3/bin/python3
from __future__ import print_function
import csv
import argparse
import os
import sys
import multiprocessing
import statistics
import colorama
import Bio.SeqIO
import Bio.Entrez
import Bio.SeqRecord
import Bio.Seq
import runbin
import random
import annotate

aligner = Bio.Align.PairwiseAligner()
aligner.mode = 'global'
aligner.open_gap_score = -2
aligner.extend_gap_score = -2
aligner.mismatch_score = -1
score_threshold = 15
chromosomes = {'1': 1, '2': 2, '3': 3, '4': 4, '5': 5, '6': 6, '7': 7, '8': 8, '9': 9, '10': 10, '11': 11, '12': 12,
               '13': 13, '14': 14, '15': 15, '16': 16, '17': 17, '18': 18, '19': 19, '20': 20, '21': 21, '22': 22,
               '23': 23, '24': 24, "X": 25, "Y": 26, "MT": 27}

class Locus(object):
    def __init__(self, row, locus_names):
        self.chrom = str(row[0])
        self.sense = row[1]
        if row[1] == "+":
            self.sensenum = 1
            self.loc = int(int(row[2]) + 3)  #TODO fix 3bp offset
        elif row[1] == "-":
            self.sensenum = 2
            self.loc = int(int(row[2]) - 3)  #TODO fix 3bp offset

        [self.gene, self.ingene, self.dist_to_gene] = locus_names[f"chr{str(self.chrom)}:{str(self.loc)}"]
        self.totalreads = row[4]
        self.totalreadsboth = row[4]
        self.similar_loci = []
        self.complement_loci = None
        self.similar_complement_loci = []
        self.sequence = '-'
        self.name = "chr" + str(self.chrom) + self.sense + str(self.loc)


class LocusCluster(object):
    def __init__(self, primary_locus, loci_list):
        self.primary = primary_locus
        self.totalreads = 0
        self.complement_loci_list = []
        for locus in loci_list:
            self.totalreads += int(locus.totalreads)
            if locus.complement_loci:
                self.complement_loci_list.append(locus.complement_loci)
        self.loci_list = loci_list
        self.complement_primary = primary_locus.complement_loci


def group(fastafile, csv_file, loci_names, outfile=None, percent=0, filteredfile=None):
    if not outfile:
        outfile = f"{os.path.splitext(os.path.realpath(csv_file))[0]}_grouped.csv"

    loci = []
    seqs = []
    for record in Bio.SeqIO.parse(fastafile, "fasta"):
        seqs.append(record)
        if seqs[-1].name.split(":")[0][-1] == "+":
            seqs[-1].id = seqs[-1].name.split(":")[0] + str(int(min(seqs[-1].name.split(":")[1].split(
                "-"))) + 53)  # [seqs[-1].name.split(":")[0][3:-1], seqs[-1].name.split(":")[0][-1], int(min(seqs[-1].name.split(":")[1].split("-")))+53]  #TODO fix 3bp offset
        elif seqs[-1].name.split(":")[0][-1] == "-":
            seqs[-1].id = seqs[-1].name.split(":")[0] + str(int(min(seqs[-1].name.split(":")[1].split(
                "-"))) + 47) #TODO fix 3bp offset
    with open(csv_file) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=",")
        line_count = 0
        for row in csv_reader:
            if len(row) > 0:
                if row[0] == "Chrom":  # headder row
                    pass
                else:
                    loci.append(Locus(row, loci_names))
                    for seq in seqs:
                        if seq.id == loci[-1].name:
                            loci[-1].sequence = seq.seq
                            break
            else:
                pass
    # connect complement loci
    loci_with_complement = []
    for locus in loci:
        loci_with_complement.append(locus)
        if loci_with_complement[-1].sense == "+":
            complement_sense = "-"
        elif loci_with_complement[-1].sense == "-":
            complement_sense = "+"
        complement_name = "chr{0}{1}{2}".format(str(loci_with_complement[-1].chrom), complement_sense,
                                                str(loci_with_complement[-1].loc))
        for other_locus in loci:
            if other_locus.name == complement_name:
                loci_with_complement[-1].complement_loci = other_locus
                loci_with_complement[-1].totalreadsboth = loci_with_complement[-1].totalreads + other_locus.totalreads
    loci = loci_with_complement
    # group similar loci
    groupped_loci = []
    for locus in loci:
        matched = False
        for i in range(len(groupped_loci)):
            if aligner.align(locus.sequence[53:], groupped_loci[i][0].sequence[53:]).score >= score_threshold:
                matched = True
                if int(locus.chrom) == 17:
                    print(f'complement? {locus.complement_loci}. total reads:{locus.totalreadsboth}. Current champ: {groupped_loci[i][0].complement_loci}. total reads:{groupped_loci[i][0].totalreadsboth}.')
                if locus.complement_loci and (not groupped_loci[i][0].complement_loci or groupped_loci[i][0].totalreadsboth < locus.totalreadsboth):  # add locus to the front of the loci group if it has the most reads
                            groupped_loci[i].insert(0, locus)
                elif not groupped_loci[i][0].complement_loci and groupped_loci[i][0].totalreadsboth < locus.totalreadsboth:  # add locus to the front of the loci group if it has the most reads
                            groupped_loci[i].insert(0, locus)
                else:
                    groupped_loci[i].append(locus)
                break
        if not matched:
            groupped_loci.append([locus])
    clustered_loci = []
    for group in groupped_loci:
        clustered_loci.append(LocusCluster(group[0], group))
    # for i in range
    #     #if locus.
    #     #[53:] 0.75 target
    #     for seq in seqs:
    #         if seq.id == loci[-1].name:
    #             loci[-1].sequence = seq.seq
    # Bio.SeqIO.write(seqs, outputfile, "fasta")
    # groupped_loci.sort(key=lambda x: x.primary.chrom, reverse=True)
    # sorted(groupped_loci, key=trial_dict.get)
    clustered_loci = sorted(clustered_loci, key=lambda x: (chromosomes[x.primary.chrom], x.primary.loc))
    for cluster in clustered_loci:
        pass
    total_mapped=0
    for group in clustered_loci:
        total_mapped += group.totalreads
    with open(outfile, "w", newline="") as csv_file:
        csv_writer = csv.writer(csv_file, delimiter=',')  # comma delimited csv file
        # write csv headder to match expected format
        csv_writer.writerow(
            ['Chrom', 'Sense', 'Loc', 'Total # IS Found', "fraction of total", "complement read", "Alternate loci", "alternate complement loci", 'Nearest Gene', "In gene", "Distance to Gene (bp)"])
        if percent > 0 and filteredfile:
            filtered = open(filteredfile, "w", newline="")
            filtered_writer = csv.writer(filtered, delimiter=',')  # comma delimited csv file
            # write csv headder to match expected format
            filtered_writer.writerow(
                ['Chrom', 'Sense', 'Loc', 'Total # IS Found', "fraction of total", "complement read", "Alternate loci",
                 "alternate complement loci", 'Nearest Gene', "In gene", "Distance to Gene (bp)"])
        for cluster in clustered_loci:
            loci_list = []
            for alt in cluster.loci_list[1:]:
                loci_list.append(f"chr{alt.chrom}{alt.sense}:{alt.loc}")
            comp_list = []
            for alt in cluster.complement_loci_list[1:]:
                comp_list.append(f"chr{alt.chrom}{alt.sense}:{alt.loc}")
            if cluster.complement_primary:
                csvline = [cluster.primary.chrom, cluster.primary.sense, cluster.primary.loc, cluster.totalreads, int(cluster.totalreads)/int(total_mapped),
                       f"chr{cluster.complement_primary.chrom}{cluster.complement_primary.sense}"
                       f":{cluster.complement_primary.loc}", loci_list, comp_list, cluster.primary.gene, cluster.primary.ingene,
                       cluster.primary.dist_to_gene]
            else:
                csvline = [cluster.primary.chrom, cluster.primary.sense, cluster.primary.loc, cluster.totalreads, int(cluster.totalreads)/int(total_mapped),
                           f"None", loci_list, comp_list, cluster.primary.gene,
                           cluster.primary.ingene,
                           cluster.primary.dist_to_gene]
            if int(cluster.totalreads)/int(total_mapped) > percent:
                csv_writer.writerow(csvline)
            elif filteredfile:
                filtered_writer.writerow(csvline)


if __name__ == "__main__":
    in_fasta = "../examples/PGK-High-C1_R1_001_retrieved_2bit.fasta"
    in_csv = "../examples/PGK-High-C1_R1_001_IS_mappings.csv"
    in_bam = "../examples/PGK-High-C1_R1_001IS.bam"
    transcripts = "../reference_datasets/annotations/refseq.transcripts.bed"
    root_name = os.path.splitext(os.path.realpath(in_csv))[0]
    out_csv = f"{root_name}_grouped.csv"
    filtered_csv = f"{root_name}_filtered.csv"
    loci_names = annotate.map_locus(transcripts, in_bam)
    group(in_fasta, in_csv, loci_names, out_csv, 0.005, filtered_csv)