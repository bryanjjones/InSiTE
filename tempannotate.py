#!/btapps/miniconda3/bin/python3
from __future__ import print_function
import pybedtools
import argparse
import os
import sys
import multiprocessing


gff = './GENCODE32.bed'
bam = './05IS.bam'
stranded = False

def featuretype_filter(feature, featuretype):
    """
    Only passes features with the specified *featuretype*
    """
    if feature[2] == featuretype:
        return True
    return False


def subset_featuretypes(featuretype):
    """
    Returns the filename containing only `featuretype` features.
    """
    return g.filter(featuretype_filter, featuretype).saveas().fn


def count_reads_in_features(features):
    """
    Callback function to count reads in features
    """
    return (
        pybedtools.BedTool(bam)
        .intersect(
            features,
            s=stranded,
            bed=True,
            stream=True,
        )
    ).count()

# Some GFF files have invalid entries -- like chromosomes with negative coords
# or features of length = 0.  This line removes them (the `remove_invalid`
# method) and saves the result in a tempfile
g = pybedtools.BedTool(gff).remove_invalid().saveas()
#g = g.sort

# Set up pool of workers
pool = multiprocessing.Pool(processes=3)

# Get separate files for introns and exons in parallel
featuretypes = ['intron', 'exon']
introns, exons = pool.map(subset_featuretypes, featuretypes)

# Since `subset_featuretypes` returns filenames, we convert to BedTool objects
# to do intersections below.
introns = pybedtools.BedTool(introns)
exons = pybedtools.BedTool(exons)

# Identify unique and shared regions using bedtools commands subtract, merge,
# and intersect.
introns =introns.sort()
exons = exons.sort()
exon_only = exons.subtract(introns).merge()
intron_only = introns.subtract(exons).merge()
intron_and_exon = (
    exons
    .intersect(introns)
    .merge()
)

    # Do intersections with BAM file in parallel. Note that we're passing filenames
    # to multiprocessing.Pool rather than BedTool objects.
features = (exon_only.fn, intron_only.fn, intron_and_exon.fn)

# Run count_reads_in_features in parallel over features
results = pool.map(count_reads_in_features, features)

labels = ('exon_only',
          'intron_only',
          'intron_and_exon')

for label, reads in zip(labels, results):
    print('{0}\t{1}'.format(label, reads))