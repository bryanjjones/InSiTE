#!/btapps/miniconda3/bin/python3
from __future__ import print_function
import csv
import pybedtools
import os


file_list = ["../30-572263308/00_fastq/PGK-High-A2_R1_001_IS_mappings_grouped.csv",
             "../30-572263308/00_fastq/PGK-Low-A3_R1_001_IS_mappings_grouped.csv",
             "../30-572263308/00_fastq/UBC-Low-B3_R1_001_IS_mappings_grouped.csv",
             "../30-572263308/00_fastq/PGK-High-C1_R1_001_IS_mappings_grouped.csv",
             "../30-572263308/00_fastq/PGK-Low-B1_R1_001_IS_mappings_grouped.csv",
             "../30-572263308/00_fastq/UBC-Low-C2_R1_001_IS_mappings_grouped.csv",
             "../30-572263308/00_fastq/PGK-High-D1_R1_001_IS_mappings_grouped.csv",
             "../30-572263308/00_fastq/PGK-Low-B2_R1_001_IS_mappings_grouped.csv",
             "../30-572263308/00_fastq/UBC-Low-D2_R1_001_IS_mappings_grouped.csv",
             "../30-572263308/00_fastq/PGK-High-D2_R1_001_IS_mappings_grouped.csv",
             "../30-572263308/00_fastq/PGK-Low-C1_R1_001_IS_mappings_grouped.csv",
             "../30-572263308/00_fastq/UBC-Med-A3_R1_001_IS_mappings_grouped.csv",
             "../30-572263308/00_fastq/PGK-High-D3_R1_001_IS_mappings_grouped.csv",
             "../30-572263308/00_fastq/PGK-Low-C2_R1_001_IS_mappings_grouped.csv",
             "../30-572263308/00_fastq/UBC-Med-B2_R1_001_IS_mappings_grouped.csv",
             "../30-572263308/00_fastq/PGK-Low-A1_R1_001_IS_mappings_grouped.csv",
             "../30-572263308/00_fastq/UBC-High-B4_R1_001_IS_mappings_grouped.csv",
             "../30-572263308/00_fastq/PGK-Low-A2_R1_001_IS_mappings_grouped.csv",
             "../30-572263308/00_fastq/UBC-Low-A2_R1_001_IS_mappings_grouped.csv"]


def startsite(feature):
    if feature.strand == "+":
        feature.stop = feature.start
    elif feature.strand == "-":
        feature.start = feature.stop
    return feature


def closest_single(files, ref='../reference_datasets/annotations/refseq.transcripts.bed', write=True):
    for file in files:
        output = None
        if write:
            out_css = f"{os.path.splitext(os.path.realpath(file))[0]}_wTSSdist.csv"
            outfile = open(out_css, "w", newline="")
            output = csv.writer(outfile, delimiter=',')
        reference = pybedtools.BedTool(ref).sort()
        reference = reference.each(startsite).sort().saveas()
        with open(file) as csv_file:
            csv_reader = csv.reader(csv_file, delimiter=",")
            line_count = 0
            for row in csv_reader:
                line_count += 1
                if len(row) > 0:
                    if row[0] == "Chrom":  # headder row
                        row.append("distance to nearest TSS")
                        if output:
                            output.writerow(row)
                    else:
                        # locus = f'chr{row[0]} {row[2]} {row[2]}'
                        locus = pybedtools.BedTool(f'chr{row[0]} {row[2]} {row[2]}', from_string=True)
                        print(locus)
                        b = locus.closest(reference, d=True, t="first")
                        for i in b:
                            print(i.count)
                            row.append(i.count)
                        if output:
                            output.writerow(row)

if __name__ == "__main__":
    closest_single(file_list, ref='../reference_datasets/annotations/refseq.transcripts.bed', write=True)

#
#         loci.append(Locus(row, loci_names))
#         if complex_loci_name:
#             loci[
#                 -1].gene = f"{loci[-1].gene.split('_')[0]}_{loci[-1].gene.split('_')[1]}"  # bed files from table viewere sometimes have complex names, this will extract the Gene ID from the name.
#         loci[-1].gene_definition = annotate.retrieve_gene_definition(loci[-1].gene)
#         for seq in seqs:
#             if seq.id == loci[-1].name:
#                 loci[-1].sequence = seq.seq
#                 break
# else:
#     pass
#


# reference = pybedtools.BedTool('../../reference_datasets/annotations/refseq.transcripts.bed').sort()
# reference = reference.each(startsite).sort().saveas()
# csv_file = open("./PGK-High-A2_R1_001_IS_mappings_grouped.csv")
# csv_reader = csv.reader(csv_file, delimiter=",")
# loci = ""
# for row in csv_reader:
#     if row[0] == "Chrom":
#         pass
#     else:
#         # loci = loci + f'\nchr{row[0]} {row[2]} {row[2]}'
#         locus = pybedtools.BedTool(f'chr{row[0]} {row[2]} {row[2]}', from_string=True)
#         # print(locus)
#         b = locus.closest(reference, d=True, t="first")
#         for i in b:
#             row.append(i.count)
#             # print(i.count)
#         print(row)
# loci_bed = pybedtools.BedTool(loci, from_string=True)
# loci_bed = loci_bed.sort().saveas()
# print(loci_bed)
# b = loci_bed.closest(reference, d=True, t="first")
# for i in b:
#     print(i.count)
#
