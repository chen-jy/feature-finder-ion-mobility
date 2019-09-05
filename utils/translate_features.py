import argparse
import csv

import pyopenms as ms
from operator import itemgetter

if __name__ == '__main__':
    print('Starting feature translation', flush=True)
    parser = argparse.ArgumentParser(description='Feature translator.')
    parser.add_argument('--input', action='store', required=True, type=str)
    parser.add_argument('--output', action='store', required=True, type=str)
    parser.add_argument('--mode', action='store', required=True, type=int)

    args = parser.parse_args()
    
    if args.mode == 0:
        csv_list = []
        with open(args.input + '.csv', 'r') as f:
            reader = csv.reader(f)
            csv_list = list(reader)
        
        features = ms.FeatureMap()
        print('Processing feature 1 of', len(csv_list) - 1)
        for i in range(1, len(csv_list)):
            if i % 50 == 0:
                print('Processing feature', i, 'of', len(csv_list) - 1)

            f = ms.Feature()
            f.setRT(float(csv_list[i][0]))
            f.setMZ(float(csv_list[i][1]))
            f.setIntensity(float(csv_list[i][2]))
            features.push_back(f)

        print('Done.')
        features.setUniqueIds()
        ms.FeatureXMLFile().store(args.output + '.featureXML', features)

    else:
        features = ms.FeatureMap()
        ms.FeatureXMLFile().load(args.input + '.featureXML', features)

        data_points = []
        for feature in features:
            data_points.append([feature.getRT(), feature.getMZ(), feature.getIntensity()])

        if args.mode == 2:
            data_points = sorted(data_points, key=itemgetter(2), reverse=True)

        with open(args.output + '.csv', 'w') as f:
            f.write('RT,m/z,Intensity\n')
            print('Processing feature 1 of', features.size())

            for i in range(features.size()):
                if i % 50 == 0:
                    print('Processing feature', i, 'of', features.size())

                f.write(str.format('{0},{1},{2}\n', data_points[i][0], data_points[i][1], data_points[i][2]))

        print('Done.')
