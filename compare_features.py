"""A feature comparison tool used to benchmark the LC-IMS-MS/MS feature finder.
"""

import argparse
import csv
from operator import itemgetter
from typing import Any, List

import pyopenms as ms

import common_utils_im as util


# For file writing
output_file = ''
# Global statistics
times_matched = [0, 0, 0]  # Zero matches, one match, multiple matches
# Program parameters
thresholds = []


def reset_stats() -> None:
    """Resets the global variables."""
    global times_matched
    times_matched = [0, 0, 0, 0]


def csv_to_list(input_filename: str) -> List[List[float]]:
    """Reads a csv file and extracts its feature data.

    The csv file must be formatted like this: Retention time, mass to charge, ion mobility

    Keyword arguments:
    input_filename: the csv file to read from

    Returns: a list of lists, where each interior list represents a feature, holding its RT, m/z,
    ion mobility, and the value False, in that order.
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
    global output_file, times_matched
    with open(output_file, 'w') as f:
        num_common = times_matched[1] + times_matched[2]
        print('Common features:', num_common)
        f.write('Common features: %d\n' % num_common)
        print('No matches:', times_matched[0])
        f.write('No matches: %d\n' % times_matched[0])
        print('One match:', times_matched[1])
        f.write('One match: %d\n' % times_matched[1])
        print('Multiple matches:', times_matched[2])
        f.write('Multiple matches: %d\n' % times_matched[2])


def compare_features(features1: Any, features2: list, im_mode: bool = True, quiet: bool = True) -> None:
    """Compares reference features against a list of found features, checking how many times each
    reference feature maps to found features.

    The resulting statistics are written to the global statistics values for printing.

    Keyword arguments:
    features1: the found features (e.g. by feature_finder_im)
    features2: the list of reference features (e.g. converted from MaxQuant output)
    im_mode: if true, compare by IM data in addition to RT and m/z data
    quiet: if true, suppress output for current progress
    """
    global times_matched, thresholds
    reset_stats()
    reset_csv_list(features2)
    if type(features1) == list:
        features1 = sorted(features1, key=itemgetter(0))
    else:
        features1.sortByRT()

    for j in range(len(features2)):
        if not quiet and j % 100 == 0:
            print('Comparing feature', j + 1, 'of', len(features2))

        num_common = 0
        first_idx = util.binary_search_left_rt(features1, features2[j][0] - thresholds[0])

        for i in range(first_idx, len(features1) if type(features1) == list else features1.size()):
            f = features1[i]
            if (f[0] if type(features1) == list else f.getRT()) > features2[j][0] + thresholds[0]:
                break

            if (im_mode and util.similar_features_im(f, features2[j], *thresholds)) or (not im_mode and
                    util.similar_features(f, features2[j], thresholds[0], thresholds[1])):
                num_common += 1
                if num_common > 1:
                    break

        if num_common == 0:
            times_matched[0] += 1
        elif num_common == 1:
            times_matched[1] += 1
        else:
            times_matched[2] += 1

    print_summary()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Feature comparison tool.')
    parser.add_argument('-i', '--in', action='store', required=True, type=str, dest='in_',
                        help='the input features (e.g. found by feature_finder_im)')
    parser.add_argument('-r', '--ref', action='store', required=True, type=str,
                        help='the reference features (e.g. found by MaxQuant)')
    parser.add_argument('-o', '--out', action='store', required=True, type=str,
                        help='the output group name (not a single filename)')

    parser.add_argument('-t', '--rt', action='store', required=False, type=float, default=5.0,
                        help='the RT threshold to use')
    parser.add_argument('-m', '--mz', action='store', required=False, type=float, default=0.01,
                        help='the m/z threshold to use')
    parser.add_argument('-z', '--im', action='store', required=False, type=float, default=0.031,
                        help='the IM threshold to use')

    parser.add_argument('--no-im', action='store_true', required=False, default=False,
                        help='do not use IM comparisons (only use RT and m/z)')
    parser.add_argument('-q', '--quiet', action='store_true', required=False, default=False,
                        help='suppress progress output when comparing features')

    args = parser.parse_args()
    output_file = args.out
    thresholds = [args.rt, args.mz, args.im]

    input_mask, ref_mask = ms.FeatureMap(), ms.FeatureMap()
    input_is_csv = True if args.in_.endswith('.csv') else False
    ref_is_csv = True if args.ref.endswith('.csv') else False

    if not input_is_csv and not args.in_.endswith('featureXML'):
        print('Error: input features must be in csv or featureXML formats')
        exit(1)
    if not ref_is_csv:  # TODO: Add more comparison formats
        print('Error: comparison currently requires reference features in csv format')
        exit(1)

    if input_is_csv: input_mask = csv_to_list(args.in_)
    else: ms.FeatureXMLFile().load(args.in_, input_mask)

    if ref_is_csv: ref_mask = csv_to_list(args.ref)
    else: ms.FeatureXMLFile().load(args.ref, ref_mask)

    compare_features(input_mask, ref_mask, not args.no_im, args.quiet)
