import argparse
import csv

import pyopenms as ms

from operator import itemgetter
from math import floor

# Snippet from im_binning
def similar_features(feature1, feature2, rt_threshold=5, mz_threshold=0.01):
    return abs(feature1[0] - feature2[0]) < rt_threshold and \
        abs(feature1[1] - feature2[1]) < mz_threshold

# Snippet from cmp_features
def binary_search_leftmost(arr, idx, target):
    l, r = 0, len(arr)
    while l < r:
        m = floor((l + r) / 2)
        if arr[m][idx] < target:
            l = m + 1
        else:
            r = m
    return l if l < len(arr) else l - 1

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='csv feature comparison by intensity.')
    parser.add_argument('--found', action='store', required=True, type=str)
    parser.add_argument('--tsv', action='store', required=True, type=str)
    parser.add_argument('--output', action='store', required=True, type=str)

    args = parser.parse_args()

    # Read input
    csv_list1, points1 = [], []
    with open(args.found + '.csv', 'r') as f:
        reader = csv.reader(f)
        csv_list1 = list(reader)

    for i in range(1, len(csv_list1)):
        points1.append([float(x) for x in csv_list1[i]])
        points1[i - 1].append(False)

    csv_list2, points2 = [], []
    with open(args.tsv + '.csv', 'r') as f:
        reader = csv.reader(f)
        csv_list2 = list(reader)

    for i in range(1, len(csv_list2)):
        points2.append([float(x) for x in csv_list2[i]])

    # Sort by m/z instead of intensity (so we can binary search)
    points2 = sorted(points2, key=itemgetter(1))

    # Begin comparison
    print('Beginning comparisons')
    for i in range(len(points1)):
        if i % 50 == 0:
            print('Processing feature', i + 1, 'of', len(points1))
        idx = binary_search_leftmost(points2, 1, points1[i][1] - 0.01)

        for j in range(idx, len(points2)):
            if similar_features(points1[i], points2[j]):
                points1[i][3] = True
                break
            if points2[j][1] + 0.01 > points1[i][1]:
                break

    with open(args.output + '.csv', 'w') as f:
        f.write('RT,m/z,Intensity\n')
        for i in range(len(points1)):
            f.write(str.format('{0},{1},{2}    {3}\n', points1[i][0], points1[i][1],
                               points1[i][2],
                               '[FOUND]' if points1[i][3] else ''))

    print('Done.')
