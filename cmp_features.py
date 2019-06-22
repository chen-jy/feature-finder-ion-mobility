import argparse

import pyopenms as ms
import numpy as np

import time
from operator import itemgetter

from im_binning import run_ff

def get_f_points(features):
    """Extracts the retention time, mass to charge, and intensity values of all data
    points in an experiment.

    Args:
        features (MSExperiment): An OpenMS MSExperiment object.

    Returns:
        list<list<float>, list<float>, list<float>>: A 3xN array of RT, MZ, and intensity
           data points, in that order.
    """
    arr_rt, arr_mz, arr_intensity = [], [], []
    for feature in features:
        arr_rt.append(feature.getRT())
        arr_mz.append(feature.getMZ())
        arr_intensity.append(feature.getIntensity())

    return [arr_rt, arr_mz, arr_intensity]

def compare_features(found_features_map, openms_features_map):
    """Compares two sets of features (namely, those found by the new feature finder and
    those found by OpenMS), and finds the number of intersecting features.

    Args:
        found_features_map (FeatureMap): An OpenMS FeatureMap object containing feaatures
            found by the new feature finder.
        openms_features_map (FeatureMap): An OpenMS FeatureMap object containing
            features found by OpenMS.
    """
    f1 = get_f_points(found_features_map)
    f2 = get_f_points(openms_features_map)

    # Sorting for easier comparisons
    found_features = sorted(f1, key=itemgetter(0, 1, 2))
    openms_features = sorted(f2, key=itemgetter(0, 1, 2))

    # Need to implement:
    # 1. Set a threshold for feature similarity
    # 2. Use binary search instead of O(N^2) passes => O(N lgN)

    num_intersecting = 0

    for ffeature in found_features:
        for ofeature in openms_features:
            if ffeature == ofeature:
                num_intersecting += 1

    print("found features:", len(found_features))
    print("openms features:", len(openms_features))
    print("intersecting features:", num_intersecting)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Feature comparison.')
    parser.add_argument('--found', action='store', required=True, type=str)
    parser.add_argument('--openms', action='store', required=False, type=str)
    parser.add_argument('--source', action='store', required=False, type=str)
    parser.add_argument('--out', action='store', required=True, type=str)

    args = parser.parse_args()
    found_features, openms_features = ms.FeatureMap(), ms.FeatureMap()

    # No OpenMS features are provided; run the feature finder on the raw data first
    if args.openms is None:
        if args.source is None:
            print("Error: a source mzML file is required for this argument set.")
            quit()

        exp = ms.MSExperiment()
        print('Loading mzML input file....................', end='', flush=True)
        ms.MzMLFile().load(args.source + '.mzML', exp)
        print('Done')

        openms_features = run_ff(exp, 'centroided')
    else:
        ms.FeatureXMLFile.load(args.openms + '.mzML', openms_features)

    ms.FeatureXMLFile.load(args.found + '.mzML', found_features)
    
    compare_features(found_features, openms_features)
