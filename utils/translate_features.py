import argparse
import csv

import pyopenms as ms

if __name__ == '__main__':
    print('Starting feature translation', flush=True)
    parser = argparse.ArgumentParser(description='Feature translator (csv to featureXML).')
    parser.add_argument('--input', action='store', required=True, type=str)
    parser.add_argument('--output', action='store', required=True, type=str)

    args = parser.parse_args()
    
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

    features.setUniqueIds()
    ms.FeatureXMLFile().store(args.output + '.featureXML', features)
