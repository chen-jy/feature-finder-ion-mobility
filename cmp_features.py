import argparse

import pyopenms as ms
import numpy as np

from operator import itemgetter
from im_binning import run_ff, similar_features

# Globals
clean_exp = ms.MSExperiment()

def get_f_points(features):
    """Extracts the retention time, mass to charge, and intensity values of all data
    points in an experiment.

    Args:
        features (MSExperiment): An OpenMS experiment.

    Returns:
        list<list<float, float, float>>: A list of lists, where each interior list
            contains the RT, MZ, and intensity of a single data point. The exterior list
            is unsorted.
    """
    data_points = []
    for feature in features:
        data_points.append([feature.getRT(), feature.getMZ(), feature.getIntensity()])

    return data_points

# Need to rewrite this to be far more efficient
@DeprecationWarning
def compare_features_old(found_fm, openms_fm, baseline_fm, truth_fm):
    """Compares two sets of features (namely, those found by the new feature finder and
    those found by OpenMS), and finds the number of intersecting features.
    Args:
        found_fm (FeatureMap): An OpenMS FeatureMap object containing feaatures found by
            the new feature finder.
        openms_fm (FeatureMap): An OpenMS FeatureMap object containing features found by
            OpenMS.
        baseline_fm (FeatureMap): An OpenMS FeatureMap object containing features found
            by the baseline (baseline.py).
        truth_fm (FeatureMap): An OpenMS FeatureMap object containing features extracted
            from the data's corresponding tsv (csv) file.
    Returns:
        FeatureMap: A list of features found in all feature maps (within a certain
            threshold).
    """
    # Let N be the greatest number of features in the input FeatureMap

    # O(4N) => O(N) (need to check if I even need to extract points!)
    f1 = get_f_points(found_fm)
    f2 = get_f_points(openms_fm)
    f3 = get_f_points(baseline_fm)
    f4 = get_f_points(truth_fm)

    # Sorting for easier comparisons
    # O(4(N lg N)) => O(N lg N)
    found_features = sorted(f1, key=itemgetter(0, 1, 2))
    openms_features = sorted(f2, key=itemgetter(0, 1, 2))
    baseline_features = sorted(f3, key=itemgetter(0, 1, 2))
    truth_features = sorted(f4, key=itemgetter(0, 1, 2))

    rt_threshold = 5
    mz_threshold = 0.01

    # Get the features similar between found and OpenMS
    num_intersecting = 0
    common_12_features = ms.FeatureMap()

    # O(N^2) : Use binary search to improve to O(N lg N)
    for ffeature in found_features:
        for ofeature in openms_features:
            f = ms.Feature()
            f.setRT(ffeature[0])
            f.setMZ(ffeature[1])
            f.setIntensity(ffeature[2])

            if ffeature == ofeature or (abs(ffeature[0] - ofeature[0]) < rt_threshold and
                                        abs(ffeature[1] - ofeature[1]) < mz_threshold):
                num_intersecting += 1
                common_12_features.push_back(f)

    num_intersecting = 0
    common_123_features = ms.FeatureMap()

    # O(N^2)
    for cfeature in common_12_features:
        for bfeature in baseline_features:
            f = ms.Feature()
            f.setRT(cfeature.getRT())
            f.setMZ(cfeature.getMZ())
            f.setIntensity(cfeature.getIntensity())

            if cfeature == bfeature or (abs(cfeature.getRT() - bfeature[0]) < rt_threshold and
                                        abs(cfeature.getMZ() - bfeature[1]) < mz_threshold):
                num_intersecting += 1
                common_123_features.push_back(f)

    num_intersecting = 0
    common_features = ms.FeatureMap()

    # O(N^2)
    for cfeature in common_123_features:
        for tfeature in truth_features:
            f = ms.Feature()
            f.setRT(cfeature.getRT())
            f.setMZ(cfeature.getMZ())
            f.setIntensity(cfeature.getIntensity())

            if cfeature == tfeature or (abs(cfeature.getRT() - tfeature[0]) < rt_threshold and
                                        abs(cfeature.getMZ() - tfeature[1]) < mz_threshold):
                num_intersecting += 1
                common_features.push_back(f)

    # Create the full venn diagram?
    print("found features:         ", len(found_features))
    print("openms features:        ", len(openms_features))
    print("baseline features:      ", len(baseline_features))
    print("truth features:         ", len(truth_features))

    print("# of fo intersecting:   ", len(common_12_features))
    print("# of fob intersecting:  ", len(common_123_features))
    print("total intersecting:     ", num_intersecting)

    # Total WC runtime: O(N + N lg N + N^2) => O(N(lg N + N))
    return common_features

def compare_features_slow(found_fm, openms_fm, baseline_fm, truth_fm):
    # Need to play around with these
    rt_threshold = 5
    mz_threshold = 0.01

    # Compare new feature finder to existing
    common_new = 0
    for ffeature in found_fm:
        for ofeature in openms_fm:
            if similar_features(ffeature, ofeature):
                common_new += 1

    # Compare Leon's feature finder to existing
    common_base = 0
    for bfeature in baseline_fm:
        for ofeature in openms_fm:
            if similar_features(bfeature, ofeature):
                common_base += 1

    # Compare new feature finder to real features
    common_newtruth = 0
    for ffeature in found_fm:
        for tfeature in truth_fm:
            if similar_features(ffeature, tfeature):
                common_newtruth += 1

    # Compare Leon's feature finder to real features
    common_basetruth = 0
    for bfeature in baseline_fm:
        for tfeature in truth_fm:
            if similar_features(bfeature, tfeature):
                common_basetruth += 1

    # Control: compare OpenMS feature finder to real features
    common_control = 0
    for ofeature in openms_fm:
        for tfeature in truth_fm:
            if similar_features(ofeature, tfeature):
                common_control += 1

    # Check if any features are found by all feature finders
    common_features = 0
    common = ms.FeatureMap()

    for ffeature in found_fm:
        for ofeature in openms_fm:
            if not similar_features(ffeature, ofeature):
                continue

            for bfeature in baseline_fm:
                if not similar_features(ffeature, bfeature):
                    continue

                for tfeature in truth_fm:
                    if not similar_features(ffeature, tfeature):
                        continue

                common_features += 1
                # Maybe choose the greater convex hull as in im_binning?
                common.push_back(ffeature)

    print("Found    - OpenMS:  ", common_new.size())
    print("Baseline - OpenMS:  ", common_base.size())
    print("Found    - tsv:     ", common_newtruth.size())
    print("Baseline - tsv:     ", common_basetruth.size())
    print("OpenMS   - tsv:     ", common_control.size())
    print("All common:         ", common_features)

    return common

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Feature comparison.')
    parser.add_argument('--found', action='store', required=True, type=str)
    parser.add_argument('--openms', action='store', required=False, type=str)
    parser.add_argument('--source', action='store', required=False, type=str)
    parser.add_argument('--baseline', action='store', required=True, type=str)
    parser.add_argument('--truth', action='store', required=True, type=str)
    parser.add_argument('--out', action='store', required=True, type=str)

    args = parser.parse_args()
    found_features, openms_features, baseline_features, truth_features = \
        ms.FeatureMap(), ms.FeatureMap(), ms.FeatureMap(), ms.FeatureMap()

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
    ms.FeatureXMLFile().load(args.truth + '.featureXML', truth_features)
    
    common_features = compare_features_slow(found_features, openms_features,
                                       baseline_features, truth_features)
    ms.FeatureXMLFile().store(args.out + '-common_features.featureXML', common_features)
