#!/btapps/miniconda3/bin/python3
from __future__ import print_function
import pybedtools
import argparse
import os
import sys
import multiprocessing
import statistics
import colorama


def featuretype_filter(feature, featuretype):
    """
    Only passes features with the specified *featuretype*
    """
    if feature[2] == featuretype:
        return True
    return False


def subset_featuretypes(featuretype, g):
    """
    Returns the filename containing only `featuretype` features.
    """
    return g.filter(featuretype_filter, featuretype).saveas().fn


def count_reads_in_features(featurefile, bam, features=None):
    """
    Callback function to count reads in features
    """

    if features is not None:
        if featurefile:
            print(colorama.Fore.RED + f'Given both file name and feature set: Using feature set, ignoring file: ' +
                  colorama.Fore.YELLOW + f'{featurefile}' + colorama.Style.RESET_ALL)
    else:
        features = pybedtools.BedTool(
            featurefile).sort().merge()  # sort and merge all features in feature file
        # to categorize whole genome as feature or non-feature
    bambed = pybedtools.BedTool(bam)
    mapped = (bambed.intersect(features, s=False, bed=True, stream=True, )).count()
    # total = pybedtools.BedTool(bam).count()
    return mapped  # , total


def map_locus(featurefile, bam):
    """
    function to map individual reads to features
    """
    bambed = pybedtools.BedTool.bam_to_bed(pybedtools.BedTool(bam)).sort()
    features = pybedtools.BedTool(featurefile).sort()
    locus_names = {}
    for nearest in pybedtools.BedTool(bambed).closest(features):

        ID = nearest.fields[9]
        # correction factor of +/- 3 should be removed, figure out where the error is coming from.
        if nearest.fields[5] == "+":
            loc = int(int(nearest.fields[1]) + 3)
        elif nearest.fields[5] == "-":
            loc = int(int(nearest.fields[1]) - 3)
        else:
            print(nearest.fields[5])
            print('is not expected + or - for orientation')
            sys.exit()
        # loc = int(nearest.fields[1])
        feature_end = int(nearest.fields[8])
        feature_start = int(nearest.fields[7])
        chrom = (nearest.fields[0])

        if feature_start <= loc <= feature_end:
            infeature = True
            distance = 0
        else:
            infeature = False
            distance = min(abs(loc - feature_end), abs(loc - feature_start))
        locus_names[f"{str(chrom)}:{str(loc)}"] = [ID, infeature, distance]
    return locus_names


    # if features is not None:
    #     if featurefile:
    #         print(colorama.Fore.RED + f'Given both file name and feature set: Using feature set, ignoring file: ' +
    #               colorama.Fore.YELLOW + f'{featurefile}' + colorama.Style.RESET_ALL)
    # else:
    #     features = pybedtools.BedTool(
    #         featurefile).sort().merge()  # sort and merge all features in feature file
    #     # to categorize whole genome as feature or non-feature
    # bambed = pybedtools.BedTool(bam)
    # mapped = (bambed.intersect(features, s=False, bed=True, stream=True, )).count()
    # # total = pybedtools.BedTool(bam).count()
    # return mapped  # , total


def intersection(a, b):
    # takes two BedTool objects and returns three BedTool objects of overlap, and non-overlapping regions of each.
    a_only = a.subtract(b).merge()
    b_only = b.subtract(a).merge()
    a_and_b = (
        a
            .intersect(b)
            .merge()
    )
    return a_only, b_only, a_and_b


def featuresets(gff, featuretypes=None, save=False):
    # Takes a gff file and list of feature types. Returns a list of BedTool objects for each specified featuretype
    # Some GFF files have invalid entries -- like chromosomes with negative coords
    # or features of length = 0.  This line removes them (the `remove_invalid`
    # method) and saves the result in a tempfile
    if featuretypes is None:
        featuretypes = ['intron', 'exon']
    g = pybedtools.BedTool(gff).remove_invalid().saveas()
    g = g.sort()
    sets = []
    for featuretype in featuretypes:
        set = pybedtools.BedTool(subset_featuretypes(featuretype, g))
        set = set.sort().merge()  # process sets by sorting and merging overlapping features
        if save:
            set.saveas(f'{featuretype}.gtf')
        sets.append(set.sort())
    return sets


# returns BedTool containing only the starting nt position when given a BedTool including strand orientation
def startsite(feature):
    if feature.strand == "+":
        feature.stop = feature.start
    elif feature.strand == "-":
        feature.start = feature.stop
    return feature


# returns list of distances and BedTool object containing list of distances between each feature in a given queryfile
# (bam) and the closest feature in the given refereence file. given featurename used for labeling output.
# If position=start, distance to startting nt of reference features will be used to calculate distance.
def closest(queryfilename, reference_filename, featurename='TSS', position="start", distances=[1000], quiet=False):
    message = []
    distancelist = []
    closedistances = []
    distancebins = []
    for i in range(len(distances)):
        closedistances.append([])
    query = pybedtools.BedTool(queryfilename).bam_to_bed().sort()
    reference = pybedtools.BedTool(reference_filename).sort()
    # if position is start, make a duplicate reference site with only startsite nt
    if position == 'start':
        reference = reference.each(startsite).sort().saveas()
    b = query.closest(reference, d=True, t="first")
    for i in b:  # make list of distances
        distancelist.append(i.count)
        for j in range(len(distances)):
            if i.count <= distances[j]:
                closedistances[j].append(i.count)
    average = statistics.mean(distancelist)
    median = statistics.median(distancelist)
    if len(distancelist) > 1:
        standarddev = statistics.stdev(distancelist)
    else:
        standarddev = 0
    if not quiet:
        print(
            f'average distance to ' + colorama.Fore.GREEN + f'{featurename}\t\t\t{round(average):,}\t' +
            colorama.Style.RESET_ALL + f'bp')
        print(
            f'standard deviation distance to ' + colorama.Fore.GREEN + f'{featurename}\t{round(standarddev):,}\t' +
            colorama.Style.RESET_ALL + f'bp')
        print(
            f'median distance to ' + colorama.Fore.GREEN + f'{featurename}\t{round(median):,}\t' +
            colorama.Style.RESET_ALL + f'bp')
        message.append(f'average distance to {featurename}\t\t\t{round(average):,}\tbp')
        message.append(f'standard deviation distance to {featurename}\t{round(standarddev):,}\tbp')
        message.append(f'median distance to {featurename}\t\t\t{round(average):,}\tbp')
    for i in range(len(closedistances)):
        distancebins.append(len(closedistances[i]))
        print(
            f'sites within ' + colorama.Fore.GREEN + f'{distances[i]}' + colorama.Style.RESET_ALL + f' bp of ' +
            colorama.Fore.GREEN + f'{featurename}\t\t{len(closedistances[i]):,}' + colorama.Style.RESET_ALL)
        message.append(
            f'sites within {distances[i]} bp of {featurename}\t\t{len(closedistances[i]):,}')
    return distancelist, average, standarddev, median, distancebins, b, "\n".join(message)


# main function to map features in IS bam file to feature file (gff/gtf/bed). if single=Falso, map to given
# featurenames list in feature file. if single=True, all features in feature file are treated as the target feature.
# If save=True, each subset of features is saved as new file.
def featuremap(gff, bam, featurenames=['intron', 'exon'], single=False, save=False, procs=3):
    pool = multiprocessing.Pool(processes=procs)
    results = []
    message = []
    if single:
        # single feature type in reference gff/gtf/bed file
        results.append(count_reads_in_features(gff, bam, features=None))
        featurenames = [featurenames]
    else:
        # multipse feature types in reference gff/gtf
        sets = featuresets(gff, featurenames, save=save)
        for feature in sets:
            results.append(count_reads_in_features(None, bam, features=feature))
    for label, reads in zip(featurenames, results):
        print(f'sites in ' + colorama.Fore.GREEN + f'{label}\t{reads:,}' + colorama.Style.RESET_ALL)
        message.append(f'sites mapped in {label}\t{reads:,}')
    # Run count_reads_in_features in parallel over features
    total = pybedtools.BedTool(bam).count()
    featurenames.append('total')
    results.append(total)
    return results, "\n".join(message)


# map bam reads to featuresets
if __name__ == "__main__":

    ap = argparse.ArgumentParser(prog=os.path.basename(sys.argv[0]),
                                 usage=__doc__)
    ap.add_argument('-a', required=True,
                    help='GFF or GTF file containing annotations')
    ap.add_argument('-m', '--gff', default=False, action='store_true',
                    help='specify if annotation file is gtf/gff file with multiple annotations instead of gtf/gff/bam '
                         'for single annotation type')
    ap.add_argument('-b', '--bam', required=True,
                    help='BAM file containing reads to be counted')
    ap.add_argument('-p', '--processes', default=3, type=int,
                    help='Number of processes to use in parallel.')
    ap.add_argument('-v', '--verbose', action='store_true',
                    help='Verbose (goes to stderr)')
    ap.add_argument('-s', '--save', action='store_true', help='save gtf of subset of feature from given gtf file')
    ap.add_argument('-f', '--feature', default=None, help='feature to map from gff/gtf file, default intron & exon')
    args = ap.parse_args()

    procs = args.processes
    gff = args.gff
    annotations = args.a
    bam = args.bam
    save = args.save
    if args.feature is None:
        featurenames = ['intron', 'exon']
    else:
        featurenames = [f'{args.feature}']

    if gff:
        featuremap(annotations, bam, featurenames=featurenames, save=save, procs=3)
    else:
        featuremap(annotations, bam, single=True, save=save, procs=3)
