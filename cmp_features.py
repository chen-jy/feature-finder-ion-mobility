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

def compare_features(found_features_map, openms_features_map, baseline_features_map):
    """Compares two sets of features (namely, those found by the new feature finder and
    those found by OpenMS), and finds the number of intersecting features.

    Args:
        found_features_map (FeatureMap): An OpenMS FeatureMap object containing feaatures
            found by the new feature finder.
        openms_features_map (FeatureMap): An OpenMS FeatureMap object containing
            features found by OpenMS.
        baseline_features_map (FeatureMap): An OpenMS FeatureMap object containing
            features found by the baseline (baseline.py).

    Returns:
        FeatureMap: A list of features found in all feature maps (within a certain
            threshold).
    """
    f1 = get_f_points(found_features_map)
    f2 = get_f_points(openms_features_map)
    f3 = get_f_points(baseline_features_map)

    # Sorting for easier comparisons
    found_features = sorted(f1, key=itemgetter(0, 1, 2))
    openms_features = sorted(f2, key=itemgetter(0, 1, 2))
    baseline_features = sorted(f3, key=itemgetter(0, 1, 2))

    rt_threshold = 5
    mz_threshold = 0.01

    num_12_intersecting = 0
    common_12_features = ms.FeatureMap()

    # Should implement a binary search in O(N lg N) instead of O(N^2)
    for ffeature in found_features:
        for ofeature in openms_features:
            f = ms.Feature()
            f.setRT(ffeature[0])
            f.setMZ(ffeature[1])
            f.setIntensity(ffeature[2])

            if ffeature == ofeature or (abs(ffeature[0] - ofeature[0]) < rt_threshold and
                                        abs(ffeature[1] - ofeature[1]) < mz_threshold):
                num_12_intersecting += 1
                common_12_features.push_back(f)

    num_intersecting = 0
    common_features = ms.FeatureMap()

    for cfeature in common_12_features:
        for bfeature in baseline_features:
            f = ms.Feature()
            f.setRT(cfeature.getRT())
            f.setMZ(cfeature.getMZ())
            f.setIntensity(cfeature.getIntensity())

            if cfeature == bfeature or (abs(cfeature.getRT() - bfeature[0]) < rt_threshold and
                                        abs(cfeature.getMZ() - bfeature[1]) < mz_threshold):
                num_intersecting += 1
                common_features.push_back(f)

    print("found features:", len(found_features))
    print("openms features:", len(openms_features))
    print("baseline features:", len(baseline_features))
    print("intersecting features:", num_intersecting)

    return common_features

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Feature comparison.')
    parser.add_argument('--found', action='store', required=True, type=str)
    parser.add_argument('--openms', action='store', required=False, type=str)
    parser.add_argument('--source', action='store', required=False, type=str)
    parser.add_argument('--baseline', action='store', required=True, type=str)
    parser.add_argument('--out', action='store', required=True, type=str)

    args = parser.parse_args()
    found_features, openms_features, baseline_features = ms.FeatureMap(), \
        ms.FeatureMap(), ms.FeatureMap()

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
        ms.FeatureXMLFile().store('openms_features.featureXML', openms_features)
    else:
        ms.FeatureXMLFile().load(args.openms + '.featureXML', openms_features)

    ms.FeatureXMLFile().load(args.found + '.featureXML', found_features)
    ms.FeatureXMLFile().load(args.baseline + '.featureXML', baseline_features)
    
    common_features = compare_features(found_features, openms_features, baseline_features)
    ms.FeatureXMLFile().store(args.out + '-common_features.featureXML', common_features)
