"""A feature comparison tool used to benchmark the LC-IMS-MS/MS feature finder.
"""

import argparse
import csv
from math import floor
from operator import itemgetter
from typing import Any, List

import pyopenms as ms

import common_utils_im as util


output_group = ''
multiple = None
num_common = 0
times_matched = [0, 0, 0, 0]  # Zero matches, one match, multiple matches
im_threshold = 0.0


def reset_stats() -> None:
    """Resets the global variables."""
    global num_common, times_matched
    num_common = 0
    times_matched = [0, 0, 0, 0]


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
        print('Missed:', times_matched[3])
        f.write('Missed: %d\n' % times_matched[3])

    # for ($i=0.005; $i -le 0.3; $i+=0.005) { python37 ...; } <- A Windows for loop
    #with open(output_group + '-im_graph.csv', 'a') as f:
        #f.write(f'{im_threshold},{num_common},{times_matched[0]},{times_matched[1]},{times_matched[2]}\n')


def convert_to_bidx(im: float) -> int:
    """DEPRECATED! Converts an ion mobility value (1/K_0) to its bin index.
    
    Keyword arguments:
    im: the ion mobility value to convert

    Returns: the bin index in which the IM value would exist.
    """
    idx = int((im - im_start) / ((im_stop - im_start) / num_bins))
    return idx if idx < num_bins else num_bins - 1


def compare(features1: list, features2: list) -> None:
    global output_group, num_common, times_matched, im_threshold
    reset_stats()
    reset_csv_list(features2)
    features2 = sorted(features2, key=itemgetter(2))

    for j in range(len(features2)):
        similar = []

        attempted = False
        for i in range(len(features1)):
            if util.similar_features_im(features1[i], features2[j], 5.0, 0.01, im_threshold):
                similar.append(features1[i])
            elif util.similar_features(features1[i], features2[j]):
                attempted = True

        if len(similar) == 0:
            times_matched[0] += 1
            if attempted:
                times_matched[3] += 1
        elif len(similar) == 1:
            num_common += 1
            times_matched[1] += 1
        else:
            num_common += 1
            times_matched[2] += 1

            multiple.write('{:.6f}, {:.6f}, {:.6f}\n'.format(*features2[j]))
            for feature in similar:
                multiple.write('{:.6f}, {:.6f}, {:.6f}\n'.format(*feature))
            multiple.write('==================================================\n')

    print_summary()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Feature comparison tool.')
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

    input_mask, ref_mask = ms.FeatureMap(), ms.FeatureMap()
    input_is_csv = True if args.in_.endswith('.csv') else False
    ref_is_csv = True if args.ref.endswith('.csv') else False

    if not input_is_csv or not ref_is_csv:
        print('Error! compare_features is currently in csv-csv only mode')
        quit()

    if input_is_csv: input_mask = csv_to_list(args.in_)
    else: ms.FeatureXMLFile().load(args.in_, input_mask)

    if ref_is_csv: ref_mask = csv_to_list(args.ref)
    else: ms.FeatureXMLFile().load(args.ref, ref_mask)

    multiple = open(output_group + '-multiple.txt', 'w')

    compare(input_mask, ref_mask)
    #compare(ref_mask, input_mask)

    multiple.close()
