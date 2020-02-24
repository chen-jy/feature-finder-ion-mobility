"""A feature comparison tool used to benchmark the LC-IMS-MS/MS feature finder.
"""

import argparse
import csv
from functools import singledispatch
from math import floor
from operator import itemgetter
from typing import Any, List

import pyopenms as ms

import common_utils_im as util


output_group = ''
num_common = 0
times_matched = [0, 0, 0]  # Zero matches, one match, multiple matches
im_threshold = 0.0


def reset_stats() -> None:
    """Resets the global variables."""
    global num_common, times_matched
    num_common = 0
    times_matched = [0, 0, 0]


def csv_to_list(input_filename: str) -> List[List[float]]:
    """Reads a csv file and extracts its feature data.

    The csv file must be formatted like this: Retention time, mass to charge, ion mobility index

    Keyword arguments:
    input_filename: the csv file to read from

    Returns: a list of lists, where each interior list represents a feature, holding its RT, m/z,
    ion mobility (1/K_0), and the value False, in that order.
    """
    csv_list, points = [], []
    with open(input_filename, 'r') as f:
        reader = csv.reader(f)
        csv_list = list(reader)
    for i in range(1, len(csv_list)):  # Skip the header
        points.append([float(x) for x in csv_list[i]])
        points[i - 1].append(False)  # If this feature is common
    return points


def reset_csv_list(csv_list: List[List[float]]) -> None:
    """Resets the common status of every feature in a list."""
    for i in range(len(csv_list)):
        csv_list[i][3] = False


def print_summary() -> None:
    """Prints and logs a summary of the most recently run comparison."""
    global output_group, num_common, times_matched
    with open(output_group + '-summary.txt', 'w') as f:
        print('Common features:', num_common)
        f.write('Common features: %d\n' % num_common)
        print('No matches:', times_matched[0])
        f.write('No matches: %d\n' % times_matched[0])
        print('One match:', times_matched[1])
        f.write('One match: %d\n' % times_matched[1])
        print('Multiple matches:', times_matched[2])
        f.write('Multiple matches: %d\n' % times_matched[2])

    # for ($i=0.005; $i -le 0.3; $i+=0.005) { python37 ...; } <- A Windows for loop
    with open(output_group + '-im_graph.csv', 'a') as f:
        f.write(f'{im_threshold},{num_common},{times_matched[0]},{times_matched[1]},{times_matched[2]}\n')


def convert_to_bidx(im: float) -> int:
    """DEPRECATED! Converts an ion mobility value (1/K_0) to its bin index.
    
    Keyword arguments:
    im: the ion mobility value to convert

    Returns: the bin index in which the IM value would exist.
    """
    idx = int((im - im_start) / ((im_stop - im_start) / num_bins))
    return idx if idx < num_bins else num_bins - 1


# A stupidly contrived way of achieving function overloading
@singledispatch
def cmp1(features2: ms.FeatureMap, features1: ms.FeatureMap) -> None:
    raise NotImplementedError


@cmp1.register
def _(features2: list, features1: ms.FeatureMap) -> None:
    global output_group, num_common, times_matched
    reset_stats()
    reset_csv_list(features2)

    common_features, missing_features = ms.FeatureMap(), ms.FeatureMap()

    for j in range(len(features2)):
        similar = []

        for feature in features1:
            if util.similar_features(feature, features2[j]):
                similar.append(feature)

        if len(similar) == 0: times_matched[0] += 1
        elif len(similar) == 1: times_matched[1] += 1
        else: times_matched[2] += 1

        if len(similar) > 0:
            max_area = util.polygon_area(similar[0].getConvexHull().getHullPoints())
            max_feature = similar[0]

            for f in similar:
                area = util.polygon_area(f.getConvexHull().getHullPoints())
                if area > max_area:
                    max_area = area
                    max_feature = f

            common_features.push_back(max_feature)
            num_common += 1
        else:
            f = ms.Feature()
            f.setRT(features2[j][0])
            f.setMZ(features2[j][1])
            f.setIntensity(features2[j][2])
            missing_features.push_back(f)

    common_features.setUniqueIds()
    ms.FeatureXMLFile().store(output_group + '-common.featureXML', common_features)
    missing_features.setUniqueIds()
    ms.FeatureXMLFile().store(output_group + '-missing.featureXML', missing_features)

    print_summary()


@singledispatch
def cmp2(features2: ms.FeatureMap, features1: list) -> None:
    raise NotImplementedError


@cmp2.register
def _(features2: list, features1: list) -> None:
    global output_group, num_common, times_matched
    reset_stats()
    reset_csv_list(features2)

    for j in range(len(features2)):
        similar = [features2[j]]

        for i in range(len(features1)):
            if util.similar_features_im(features1[i], features2[j], 5.0, 0.01, im_threshold):
                if not features1[i][3]:
                    features1[i][3] = True
                    num_common += 1
                similar.append(features1[i])

        if len(similar) == 1:
            times_matched[0] += 1
        elif len(similar) == 2:
            times_matched[1] += 1
        else:
            times_matched[2] += 1

    with open(output_group + '.csv', 'w') as f:
        f.write('RT,m/z,Intensity\n')
        for j in range(len(features2)):
            f.write(str.format('{0},{1},{2},{3}\n', features2[j][0], features2[j][1], features2[j][2],
                                'FOUND' if features2[j][3] else ''))

    print_summary()


def compare(features1: list, features2: list) -> None:
    global output_group, num_common, times_matched, im_threshold
    reset_stats()
    reset_csv_list(features2)

    for j in range(len(features2)):
        similar = []

        for i in range(len(features1)):
            if util.similar_features_im(features1[i], features2[j], 5.0, 0.01, im_threshold):
                similar.append(features1[i])

        if len(similar) == 0:
            times_matched[0] += 1
        elif len(similar) == 1:
            num_common += 1
            times_matched[1] += 1
        else:
            num_common += 1
            times_matched[2] += 1

    print_summary()


@singledispatch
def compare_features(features1: ms.FeatureMap, features2: Any) -> None:
    cmp1(features2, features1)


@compare_features.register
def _(features1: list, features2: Any) -> None:
    cmp2(features2, features1)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Feature comparison tool.')  # Currently in IM REQUIRED mode
    parser.add_argument('-i', '--in', action='store', required=True, type=str, dest='in_',
                        help='the input features (e.g. found by feature_finder_im)')
    parser.add_argument('-r', '--ref', action='store', required=True, type=str,
                        help='the reference features (e.g. found by MaxQuant)')
    parser.add_argument('-o', '--out', action='store', required=True, type=str,
                        help='the output group name (not a single filename)')
    parser.add_argument('--threshold', action='store', required=False, type=float, default=0.1)

    args = parser.parse_args()
    output_group = args.out
    im_threshold = args.threshold

    print('WARNING: this tool currently only supports im-im comparisons for list-list\n')

    input_mask, ref_mask = ms.FeatureMap(), ms.FeatureMap()
    input_is_csv = True if args.in_.endswith('.csv') else False
    ref_is_csv = True if args.ref.endswith('.csv') else False

    if input_is_csv: input_mask = csv_to_list(args.in_)
    else: ms.FeatureXMLFile().load(args.in_, input_mask)

    if ref_is_csv: ref_mask = csv_to_list(args.ref)
    else: ms.FeatureXMLFile().load(args.ref, ref_mask)

    compare(input_mask, ref_mask)
    #compare_features(ref_mask, input_mask)
