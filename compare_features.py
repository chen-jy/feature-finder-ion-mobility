import argparse
import csv

import pyopenms as ms

from operator import itemgetter
from math import floor

def similar_features(feature1, feature2, rt_threshold=5, mz_threshold=0.01):
    if isinstance(feature1, list) and isinstance(feature2, list):
        return (abs(feature1[0] - feature2[0]) < rt_threshold and
                abs(feature1[1] - feature2[1]) < mz_threshold)
    if isinstance(feature1, ms.Feature) and isinstance(feature2, ms.Feature):
        return (abs(feature1.getRT() - feature2.getRT()) < rt_threshold and
                abs(feature1.getMZ() - feature2.getMZ()) < mz_threshold)
    if isinstance(feature1, list) and isinstance(feature2, ms.Feature):
        return (abs(feature1[0] - feature2.getRT()) < rt_threshold and
                abs(feature1[1] - feature2.getMZ()) < mz_threshold)
    if isinstance(feature1, ms.Feature) and isinstance(feature2, list):
        return similar_features(feature2, feature1, rt_threshold, mz_threshold)

    return False

def hull_area(hull):
    area = 0.0
    for i in range(len(hull)):
        area += hull[i][0] * hull[(i + 1) % len(hull)][1]
        area -= hull[i][1] * hull[(i + 1) % len(hull)][0]
    return abs(area) / 2.0

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Feature comparison tool.')
    parser.add_argument('--input', action='store', required=True, type=str)
    parser.add_argument('--output', action='store', required=True, type=str)
    parser.add_argument('--ref', action='store', required=True, type=str)
    args = parser.parse_args()

    num_common = 0
    matched = [0, 0, 0] # Zero times, once, multiple times

    #TODO: add the option for ref to be featureXML, not just csv.

    if args.input[-3:] == 'csv':
        # Read found features
        csv_list1, points1 = [], []
        with open(args.input, 'r') as f:
            reader = csv.reader(f)
            csv_list1 = list(reader)

        for i in range(1, len(csv_list1)):
            points1.append([float(x) for x in csv_list1[i]])

        # Read reference features
        csv_list2, points2 = [], []
        with open(args.ref, 'r') as f:
            reader = csv.reader(f)
            csv_list2 = list(reader)

        for i in range(1, len(csv_list2)):
            points2.append([float(x) for x in csv_list2[i]])
            # To indicate whether this feature was found or not
            points2[i - 1].append(False)

        print('Beginning comparisons')
        for i in range(len(points2)):
            if i % 50 == 0:
                print('Processing feature', i + 1, 'of', len(points2))

            num_similar = 0
            for j in range(len(points1)):
                if similar_features(points2[i], points1[j]):
                    if points2[i][3] == False:
                        points2[i][3] = True
                        num_common += 1
                    num_similar += 1

            if num_similar == 0:
                matched[0] += 1
            elif num_similar == 1:
                matched[1] += 1
            else:
                matched[2] += 1

        points2 = sorted(points2, key=itemgetter(2), reverse=True)

        with open(args.output, 'w') as f:
            f.write('RT,m/z,Intensity\n')
            for i in range(len(points2)):
                f.write(str.format('{0},{1},{2},{3}\n', points2[i][0], points2[i][1],
                                   points2[i][2],
                                   'COMMON' if points2[i][3] else 'MISSING'))

        print(num_common, ' features common. (', sum(matched), ')', sep='')

    elif args.input[-10:] == 'featureXML':
        found_features = ms.FeatureMap()
        ms.FeatureXMLFile().load(args.input, found_features)

        # Read reference features
        csv_list, points = [], []
        with open(args.ref, 'r') as f:
            reader = csv.reader(f)
            csv_list = list(reader)

        for i in range(1, len(csv_list)):
            points.append([float(x) for x in csv_list[i]])

        common_features, missing_features, unique_features, all_features = \
            ms.FeatureMap(), ms.FeatureMap(), ms.FeatureMap(), ms.FeatureMap()

        print('Beginning comparisons')
        for i in range(len(points)):
            if i % 50 == 0:
                print('Processing feature', i + 1, 'of', len(points))

            similar = []
            for j in range(found_features.size()):
                if similar_features(points[i], found_features[j]):
                    similar.append(found_features[j])

            if len(similar) == 0:
                matched[0] += 1
            elif len(similar) == 1:
                matched[1] += 1
            else:
                matched[2] += 1

            if len(similar) > 0:
                max_feature = similar[0]
                max_area = hull_area(max_feature.getConvexHull().getHullPoints())

                for f in similar:
                    ha = hull_area(f.getConvexHull().getHullPoints())
                    if ha > max_area:
                        max_area = ha
                        max_feature = f

                common_features.push_back(max_feature)
                all_features.push_back(max_feature)
                num_common += 1
            else:
                f = ms.Feature()
                f.setRT(points[i][0])
                f.setMZ(points[i][1])
                f.setIntensity(points[i][2])

                missing_features.push_back(f)
                all_features.push_back(f)

        for f in found_features:
            if f not in common_features:
                unique_features.push_back(f)
                
        all_features.setUniqueIds()
        ms.FeatureXMLFile().store(args.output + '-all.featureXML', all_features)
        missing_features.setUniqueIds()
        ms.FeatureXMLFile().store(args.output + '-missing.featureXML', missing_features)
        unique_features.setUniqueIds()
        ms.FeatureXMLFile().store(args.output + '-unique.featureXML', unique_features)

        common_features.setUniqueIds()
        ms.FeatureXMLFile().store(args.output + '-common.featureXML', common_features)
        print(num_common, ' features common. (', sum(matched), ')', sep='')

    else:
        print("Error: input file format must be either csv or featureXML")
        exit(1)

    print("No matches:", matched[0])
    print("One match:", matched[1])
    print("Multiple matches:", matched[2])
