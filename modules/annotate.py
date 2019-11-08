#!/btapps/miniconda3/bin/python3
from __future__ import print_function
import pybedtools
import argparse
import os
import sys
import multiprocessing




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
def features (featuretypes, procs, gff, bam)
    # Some GFF files have invalid entries -- like chromosomes with negative coords
    # or features of length = 0.  This line removes them (the `remove_invalid`
    # method) and saves the result in a tempfile
    g = pybedtools.BedTool(gff).remove_invalid().saveas()
    #g = g.sort

    # Set up pool of workers
    pool = multiprocessing.Pool(processes=procs)

    # Get separate files for introns and exons in parallel

    a, b = pool.map(subset_featuretypes, featuretypes)

    # Since `subset_featuretypes` returns filenames, we convert to BedTool objects
    # to do intersections below.
    a = pybedtools.BedTool(a)
    b = pybedtools.BedTool(b)

    # Identify unique and shared regions using bedtools commands subtract, merge,
    # and intersect.
    a =a.sort()
    b = b.sort()
    a_only = a.subtract(b).merge()
    b_only = b.subtract(a).merge()
    a_and_b = (
        a
        .intersect(b)
        .merge()
    )
    print(f'{len(pybedtools.BedTool(bam))} reads')
    print(f'{len(g)} annotations')
    print(f'{len(a)} {featuretypes[0]}')
    print(f'{len(b)} {featuretypes[1]}')
    print(f'{len(a_only)} {featuretypes[0]} onlys')
    print(f'{len(b_only)} {featuretypes[1]} onlys')
    print(f'{len(a_and_b)} {featuretypes[0]} and {featuretypes[1]}')

    # Do intersections with BAM file in parallel. Note that we're passing filenames
    # to multiprocessing.Pool rather than BedTool objects.
    features = (a_only.fn, b_only.fn, a_and_b.fn)

    # Run count_reads_in_features in parallel over features
    results = pool.map(count_reads_in_features, features)

    
    return results, labels


if __name__ == "__main__":

    ap = argparse.ArgumentParser(prog=os.path.basename(sys.argv[0]),
                                 usage=__doc__)
    ap.add_argument('--gff', required=True,
                    help='GFF or GTF file containing annotations')
    ap.add_argument('--bam', required=True,
                    help='BAM file containing reads to be counted')
    ap.add_argument('--stranded', action='store_true',
                    help='Use strand-specific merging and overlap. '
                         'Default is to ignore strand')
    ap.add_argument('--processes', default=1, type=int,
                    help='Number of processes to use in parallel.')
    ap.add_argument('-v', '--verbose', action='store_true',
                    help='Verbose (goes to stderr)')
    args = ap.parse_args()

    gff = args.gff
    bam = args.bam
    stranded = args.stranded
    procs = args.processes
    featuretypes = ['intron', 'exon']


    if args.processes > 3:
        print(
            "Only need 3 processes (one each for exon, intron, both), so "
            "resetting processes from {0} to 3".format(args.processes)
        )
        args.processes = 3
    features (featuretypes, procs, gff, bam)

    labels = (f'{featuretypes[0]}_only',
              f'{featuretypes[1]}_only',
              f'{featuretypes[0]}_and_{featuretypes[1]}')

    for label, reads in zip(labels, results):
    print('{0}\t{1}'.format(label, reads))