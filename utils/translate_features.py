import argparse
import csv

import pyopenms as ms
from operator import itemgetter

def checkFloat(val):
    try:
        return float(val)
    except ValueError:
        return False

if __name__ == '__main__':
    print('Starting feature translation', flush=True)
    parser = argparse.ArgumentParser(description='Feature translator.')
    parser.add_argument('--input', action='store', required=True, type=str)
    parser.add_argument('--output', action='store', required=True, type=str)

    args = parser.parse_args()
    
    if args.input[-3:] == 'csv':
        csv_list = []
        with open(args.input, 'r') as f:
            reader = csv.reader(f)
            csv_list = list(reader)
        
        features = ms.FeatureMap()
        print('Processing feature 1 of', len(csv_list) - 1)
        for i in range(1, len(csv_list)):
            if i % 50 == 0:
                print('Processing feature', i, 'of', len(csv_list) - 1)
                
            rt = checkFloat(csv_list[i][0])
            mz = checkFloat(csv_list[i][1])
            intensity = checkFloat(csv_list[i][2])
            if rt is False or mz is False or intensity is False:
                continue
            
            f = ms.Feature()
            f.setRT(float(csv_list[i][0]))
            f.setMZ(float(csv_list[i][1]))
            f.setIntensity(float(csv_list[i][2]))
            features.push_back(f)

        features.setUniqueIds()
        ms.FeatureXMLFile().store(args.output, features)
        print('Done.')

    elif args.input[-10:] == 'featureXML':
        features = ms.FeatureMap()
        ms.FeatureXMLFile().load(args.input, features)

        data_points = []
        for feature in features:
            data_points.append([feature.getRT(), feature.getMZ(), feature.getIntensity()])

        with open(args.output, 'w') as f:
            f.write('RT,m/z,Intensity\n')
            print('Processing feature 1 of', features.size())

            for i in range(features.size()):
                if i % 50 == 0:
                    print('Processing feature', i + 1, 'of', features.size())

                f.write(str.format('{0},{1},{2}\n', data_points[i][0], data_points[i][1],
                                   data_points[i][2]))

        print('Done.')

    else:
        print("Error: input file format must be csv or featureXML")
