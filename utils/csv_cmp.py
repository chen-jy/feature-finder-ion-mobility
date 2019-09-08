import argparse
import csv

import pyopenms as ms

from operator import itemgetter
from math import floor

# Modified snippet from im_binning
def similar_features(feature1, feature2, rt_threshold=5, mz_threshold=0.01):
    if isinstance(feature1, list) and isinstance(feature2, list):
        return (abs(feature1[0] - feature2[0]) < rt_threshold and
                abs(feature1[1] - feature2[1]) < mz_threshold)
    if isinstance(feature1, list) and isinstance(feature2, ms.Feature):
        return (abs(feature1[0] - feature2.getRT()) < rt_threshold and
                abs(feature1[1] - feature2.getMZ()) < mz_threshold)
    if isinstance(feature1, ms.Feature) and isinstance(feature2, list):
        return similar_features(feature2, feature1, rt_threshold, mz_threshold)

    return False

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

    if args.found[-3:] == 'csv':
        # Read input
        csv_list1, points1 = [], []
        with open(args.found, 'r') as f:
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
        num_common = 0

        for i in range(len(points1)):
            if i % 50 == 0:
                print('Processing feature', i + 1, 'of', len(points1))
            idx = binary_search_leftmost(points2, 1, points1[i][1] - 0.01)

            for j in range(idx, len(points2)):
                if similar_features(points1[i], points2[j]):
                    points1[i][3] = True
                    num_common += 1
                    break
                if points2[j][1] > points1[i][1] + 0.01:
                    break

        with open(args.output + '.csv', 'w') as f:
            f.write('RT,m/z,Intensity\n')
            for i in range(len(points1)):
                f.write(str.format('{0},{1},{2}    {3}\n', points1[i][0], points1[i][1],
                                   points1[i][2],
                                   '[FOUND]' if points1[i][3] else ''))

        print(num_common, 'features common.')

    else:
        found_features = ms.FeatureMap()
        ms.FeatureXMLFile().load(args.found, found_features)

        dict = []
        for i in range(found_features.size()):
            dict.append([found_features[i].getIntensity(), i])

        dict = sorted(dict, key=itemgetter(0), reverse=True)

        csv_list, points = [], []
        with open(args.tsv + '.csv', 'r') as f:
            reader = csv.reader(f)
            csv_list = list(reader)

        for i in range(1, len(csv_list)):
            points.append([float(x) for x in csv_list[i]])

        points = sorted(points, key=itemgetter(1))

        common_features = ms.FeatureMap()

        print('Beginning comparisons')
        num_common = 0

        for i in range(len(dict)):
            if i % 50 == 0:
                print('Processing feature', i + 1, 'of', len(dict))

            k = dict[i][1]
            idx = binary_search_leftmost(points, 1, found_features[k].getMZ() - 0.01)

            for j in range(idx, len(points)):
                if similar_features(points[j], found_features[k]):
                    common_features.push_back(found_features[k])
                    num_common += 1
                    break
                if points[j][1] > found_features[k].getMZ() + 0.01:
                    break

        common_features.setUniqueIds()
        ms.FeatureXMLFile().store(args.output + '.featureXML', common_features)
        print(num_common, 'features common.')
