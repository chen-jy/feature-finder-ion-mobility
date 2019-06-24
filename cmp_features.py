import argparse

import pyopenms as ms
import numpy as np

import time
from operator import itemgetter

from im_binning import run_ff

# Globals
clean_exp = ms.MSExperiment()

def get_f_points(features):
    """Extracts the retention time, mass to charge, and intensity values of all data
    points in an experiment.

    Args:
        features (MSExperiment): An OpenMS MSExperiment object.

    Returns:
        list<list<float, float, float>>: A list of lists, where each interior list
            contains the RT, MZ, and intensity of a single data point. The exterior list
            is unsorted.
    """
    data_points = []
    for feature in features:
        data_points.append([feature.getRT(), feature.getMZ(), feature.getIntensity()])

    return data_points

def compare_features(found_features_map, openms_features_map):
    """Compares two sets of features (namely, those found by the new feature finder and
    those found by OpenMS), and finds the number of intersecting features.

    Args:
        found_features_map (FeatureMap): An OpenMS FeatureMap object containing feaatures
            found by the new feature finder.
        openms_features_map (FeatureMap): An OpenMS FeatureMap object containing
            features found by OpenMS.

    Returns:
        list<Feature>: A list of features found in both feature maps (within a certain
            threshold).
    """
    f1 = get_f_points(found_features_map)
    f2 = get_f_points(openms_features_map)

    # Sorting for easier comparisons
    found_features = sorted(f1, key=itemgetter(0, 1, 2))
    openms_features = sorted(f2, key=itemgetter(0, 1, 2))

    rt_threshold = 0.01
    mz_threshold = 0.01
    intensity_threshold = 0.01

    num_intersecting = 0
    common_features = []

    # Should implement a binary search in O(N lg N) instead of O(N^2)
    for ffeature in found_features:
        for ofeature in openms_features:
            if ffeature == ofeature:
                num_intersecting += 1
                common_features.append(ffeature)
            elif abs(ffeature[0] - ofeature[0]) < rt_threshold and \
                 abs(ffeature[1] - ofeature[1]) < mz_threshold and \
                 abs(ffeature[2] - ofeature[2]) < intensity_threshold:
                num_intersecting += 1
                common_features.append(ffeature)

    print("found features:", len(found_features))
    print("openms features:", len(openms_features))
    print("intersecting features:", num_intersecting)

    return common_features

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

        # Remove MS2 scans
        spectra = exp.getSpectra()
        for i in range(len(spectra)):
            spec = spectra[i]
            if spec.getMSLevel() == 1:
                clean_exp.addSpectrum(spec)

        openms_features = run_ff(clean_exp, 'centroided')
    else:
        ms.FeatureXMLFile().load(args.openms + '.featureXML', openms_features)

    ms.FeatureXMLFile().load(args.found + '.featureXML', found_features)
    
    common_features = compare_features(found_features, openms_features)
    ms.FeatureXMLFile().store(args.out + 'common_features.featureXML', common_features)
