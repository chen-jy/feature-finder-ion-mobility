import argparse

import pyopenms as ms
import numpy as np

from operator import itemgetter
from math import floor
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

def get_fs_points(features, key):
    """Does the same thing as get_f_points(), but also sorts the data according to the
    key.
    """
    data_points = get_f_points(features)
    return sorted(data_points, key=key)

def get_fsi_points(features, key):
    """Does the same thing as get_f_points(), but also keeps a pointer (as an integer
    index) to the original OpenMS Feature object, and sorts the data according to the
    key.
    """
    data_points = get_f_points(features)
    for i in range(len(data_points)):
        data_points[i].append(i)
    return sorted(data_points, key=key)

def point_to_feature(point):
    """Translates a feature, represented by a list, into an OpenMS Feature object.

    Args:
        point (list<float>): Three floats, representing RT, MZ, and intensity, in that
            order.

    Returns:
        Feature: An equivalent OpenMS Feature object that only contains RT, MZ, and
            intensity information.
    """
    f = ms.Feature()
    f.setRT(point[0])
    f.setMZ(point[1])
    f.setIntensity(point[2])
    return f

def binary_search_leftmost(arr, idx, target):
    """Binary search a single column (idx) in a 2D array.

    <arr> must be sorted in <idx>.

    Args:
        arr (list<list<float>>): The 2D array to be searched.
        idx (int): The index of the column of the array to be searched.
        target (float): The target number to be searched for.

    Returns:
        int: The index of the leftmost element equal to the target, if it exists.
            Otherwise, the number of elements less than the target is returned.
    """
    l, r = 0, len(arr)
    while l < r:
        m = floor((l + r) / 2)
        if arr[m][idx] < target:
            l = m + 1
        else:
            r = m
    return l if l < len(arr) else l - 1

def compare_features(found_fm, openms_fm, baseline_fm, truth_fm, rt_threshold=5,
                     mz_threshold=0.01, brute_force=False):
    """Compares the sets of features found by the new feature finder (im_binning.py),
    OpenMS's FeatureFinderCentroided, baseline.py, and the souce mzML's corresponding
    tsv file. Performs various comparisons and prints their results to stdout.

    Args:
        found_fm (FeatureMap): An OpenMS feature map from the new feature finder.
        openms_fm (FeatureMap): An OpenMS feature map from FeatureFinderCentroided.
        baseline_fm (FeatureMap): An OpenMS feature map from the baseline.
        truth_fm (FeatureMap): An OpenMS feature map from the corresponding tsv file.
        rt_threshold (float): The retention time threshold to determine if features are
            similar.
        mz_threshold (float): The mass to charge threshold to determine if features are
            similar.
        brute_force (bool): Whether or not to hide log messages (hidden if True).

    Returns:
        FeatureMap: An OpenMS feature map containing features found in all four input
            feature maps.
    """
    if not brute_force:
        print('Extracting data points from feature maps', flush=True, end='\n\n')
    found = get_fsi_points(found_fm, itemgetter(1, 0))
    openms = get_fsi_points(openms_fm, itemgetter(1, 0))
    baseline = get_fsi_points(baseline_fm, itemgetter(1, 0))
    truth = get_fsi_points(truth_fm, itemgetter(1, 0))

    # Compare new feature finder to existing
    if not brute_force:
        print('Comparing found features to OpenMS features', flush=True)
    common_new = 0

    for i in range(len(found)):
        j = binary_search_leftmost(openms, 1, found[i][1] - mz_threshold)

        if openms[j][1] == found[i][1] - mz_threshold:
            for k in range(j, len(openms)):
                if openms[k][1] > found[i][1] + mz_threshold:
                    break
                if similar_features(found[i], openms[k]):
                    common_new += 1
        else:
            for k in range(j, -1, -1):
                if openms[k][1] < found[i][1] - mz_threshold:
                    break
                if similar_features(found[i], openms[k]):
                    common_new += 1

    # Compare baseline to existing feature finder
    if not brute_force:
        print('Comparing baseline features to OpenMS features', flush=True)
    common_base = 0

    for i in range(len(baseline)):
        j = binary_search_leftmost(openms, 1, baseline[i][1] - mz_threshold)

        if openms[j][1] == baseline[i][1] - mz_threshold:
            for k in range(j, len(openms)):
                if openms[k][1] > baseline[i][1] + mz_threshold:
                    break
                if similar_features(baseline[i], openms[k]):
                    common_base += 1
        else:
            for k in range(j, -1, -1):
                if openms[k][1] < baseline[i][1] - mz_threshold:
                    break
                if similar_features(baseline[i], openms[k]):
                    common_base += 1

    # Compare new feature finder to tsv real features
    if not brute_force:
        print('Comparing found features to tsv features', flush=True)
    common_newtruth = 0
    good_features = ms.FeatureMap()

    for i in range(len(found)):
        j = binary_search_leftmost(truth, 1, found[i][1] - mz_threshold)

        if truth[j][1] == found[i][1] - mz_threshold:
            for k in range(j, len(truth)):
                if truth[k][1] > found[i][1] + mz_threshold:
                    break
                if similar_features(found[i], truth[k]):
                    common_newtruth += 1
                    #good_features.push_back(point_to_feature(found[i]))
                    good_features.push_back(found_fm[found[i][3]])
        else:
            for k in range(j, -1, -1):
                if truth[k][1] < found[i][1] - mz_threshold:
                    break
                if similar_features(found[i], truth[k]):
                    common_newtruth += 1
                    #good_features.push_back(point_to_feature(found[i]))
                    good_features.push_back(found_fm[found[i][3]])

    # Compare baseline features to tsv real features
    if not brute_force:
        print('Comparing baseline features to tsv features', flush=True)
    common_basetruth = 0

    for i in range(len(baseline)):
        j = binary_search_leftmost(truth, 1, baseline[i][1] - mz_threshold)

        if truth[j][1] == baseline[i][1] - mz_threshold:
            for k in range(j, len(truth)):
                if truth[k][1] > baseline[i][1] + mz_threshold:
                    break
                if similar_features(baseline[i], truth[k]):
                    common_basetruth += 1
        else:
            for k in range(j, -1, -1):
                if truth[k][1] < baseline[i][1] - mz_threshold:
                    break
                if similar_features(baseline[i], truth[k]):
                    common_basetruth += 1

    # Control test: compare OpenMS features to tsv real features
    if not brute_force:
        print('Comparing OpenMS features to tsv features', flush=True)
    common_control = 0

    for i in range(len(openms)):
        j = binary_search_leftmost(truth, 1, openms[i][1] - mz_threshold)

        if truth[j][1] == openms[i][1] - mz_threshold:
            for k in range(j, len(truth)):
                if truth[k][1] > openms[i][1] + mz_threshold:
                    break
                if similar_features(openms[i], truth[k]):
                    common_control += 1
        else:
            for k in range(j, -1, -1):
                if truth[k][1] < openms[i][1] - mz_threshold:
                    break
                if similar_features(openms[i], truth[k]):
                    common_control += 1

    # Check for features found by all feature finders
    if not brute_force:
        print('Comparing all features', flush=True)
    common_features, common_triad = 0, 0
    common = ms.FeatureMap()

    for i in range(len(truth)):
        a = binary_search_leftmost(found, 1, truth[i][1] - mz_threshold)
        found_found = False

        if found[a][1] == truth[i][1] - mz_threshold:
            for j in range(a, len(found)):
                if found[j][1] > truth[i][1] + mz_threshold:
                    break
                if similar_features(truth[i], found[j]):
                    found_found = True
                    break

        else:
            for j in range(a, -1, -1):
                if found[j][1] < truth[i][1] - mz_threshold:
                    break
                if similar_features(truth[i], found[j]):
                    found_found = True
                    break

        if not found_found:
            continue

        b = binary_search_leftmost(openms, 1, truth[i][1] - mz_threshold)
        found_openms = False

        if openms[b][1] == truth[i][1] - mz_threshold:
            for j in range(b, len(openms)):
                if openms[j][1] > truth[i][1] + mz_threshold:
                    break
                if similar_features(truth[i], openms[j]):
                    found_openms = True
                    break
        else:
            for j in range(b, -1, -1):
                if openms[j][1] < truth[i][1] - mz_threshold:
                    break
                if similar_features(truth[i], openms[j]):
                    found_openms = True
                    break

        if not found_openms:
            continue

        # Feature found by the new feature finder, FeatureFinderCentroided, and the tsv
        common_triad += 1

        c = binary_search_leftmost(baseline, 1, truth[i][1] - mz_threshold)
        found_baseline = False

        if baseline[c][1] == truth[i][1] - mz_threshold:
            for j in range(c, len(baseline)):
                if baseline[j][1] > truth[i][1] + mz_threshold:
                    break
                if similar_features(truth[i], baseline[j]):
                    found_baseline = True
                    break
        else:
            for j in range(c, -1, -1):
                if baseline[j][1] < truth[i][1] - mz_threshold:
                    break
                if similar_features(truth[i], baseline[j]):
                    found_baseline = True
                    break

        if not found_baseline:
            continue

        common_features += 1
        common.push_back(truth[i])

    if not brute_force:
        print('\nFound features:    ', len(found))
        print('Baseline features: ', len(baseline))
        print('OpenMS features:   ', len(openms))
        print('tsv features:      ', len(truth))

        print('\nFound - OpenMS:    ', common_new)
        print('Baseline - OpenMS: ', common_base)
        print('Found - tsv:       ', common_newtruth)
        print('Baseline - tsv:    ', common_basetruth)
        print('OpenMS - tsv:      ', common_control)
        print('All common (-BL):  ', common_triad)
        print('All common:        ', common_features, end='\n\n')
    else:
        print(rt_threshold, mz_threshold)
        print(common_new, common_base, common_newtruth, common_basetruth, common_control,
              common_triad, common_features)

    return good_features

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Feature comparison.')
    parser.add_argument('--found', action='store', required=True, type=str)
    parser.add_argument('--openms', action='store', required=False, type=str)
    parser.add_argument('--source', action='store', required=False, type=str)
    parser.add_argument('--baseline', action='store', required=True, type=str)
    parser.add_argument('--truth', action='store', required=True, type=str)
    parser.add_argument('--outdir', action='store', required=True, type=str)
    parser.add_argument('--brute_force', action='store', required=False)

    args = parser.parse_args()
    found_features, openms_features, baseline_features, truth_features = \
        ms.FeatureMap(), ms.FeatureMap(), ms.FeatureMap(), ms.FeatureMap()

    brute_force = False
    if args.brute_force is not None:
        brute_force = True

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
        ms.FeatureXMLFile().store(args.outdir + '/openms.featureXML', openms_features)
    else:
        ms.FeatureXMLFile().load(args.openms + '.featureXML', openms_features)

    ms.FeatureXMLFile().load(args.found + '.featureXML', found_features)
    ms.FeatureXMLFile().load(args.baseline + '.featureXML', baseline_features)
    ms.FeatureXMLFile().load(args.truth + '.featureXML', truth_features)
    
    if not brute_force:
        print('Features loaded, beginning comparison', flush=True)
        common_features = compare_features(found_features, openms_features,
                                       baseline_features, truth_features)
        common_features.setUniqueIds()
        ms.FeatureXMLFile().store(args.outdir + '/common.featureXML', common_features)
    else:
        rt_start = 1
        mz_start = 0.005

        for rt_threshold in np.arange(rt_start, 12, 0.5):
            for mz_threshold in np.arange(mz_start, 0.5, 0.05):
                compare_features(found_features, openms_features, baseline_features,
                                 truth_features, rt_threshold, mz_threshold, True)
