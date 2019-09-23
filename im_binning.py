import argparse

import pyopenms as ms
import numpy as np

from operator import itemgetter
import peak_picker

# Globals and constants
bins, exps = [], []
first_im, last_im, delta_im = 0, 0, 0
num_bins, bin_size = 0, 0

# For the second pass (shift the bins by 50%)
bins2, exps2 = [], []
offset_im = 0

def get_points(spec):
    """Data preprocessing to extract the retention time, mass to charge, intensity,
    and ion mobility for each peak in a spectrum.

    Args:
        spec (MSSpectrum): The OpenMS spectrum to extract data from.

    Returns:
        list<list<float, float, float, float>>: A list of lists, where each interior
        list holds RT, MZ, intensity, and IM information (in that order) for a single
        peak in the spectrum. The exterior list is unsorted.
    """
    point_data = zip(*spec.get_peaks(), spec.getFloatDataArrays()[0])
    return [[spec.getRT(), mz, intensity, im] for mz, intensity, im in point_data]

def get_points_pp(spec_pp, spec):
    """Does the same thing as get_points(), but for a spectrum that has been peak-picked
    (and no longer has IM information). The two spectra must correspond to one another.
    """
    point_data = zip(*spec_pp.get_peaks(), spec.getFloatDataArrays()[0])
    return [[spec.getRT(), mz, intensity, im] for mz, intensity, im in point_data]

def get_extrema(spectra):
    """Finds the smallest and largest IM values in an experiment.

    Args:
        spectra (list<MSSpectrum>): A list of OpenMS spectra.

    Returns:
        tuple<float, float>: The smallest and largest IM values in the experiment.
    """
    smallest_im, largest_im = float('inf'), -1.0

    for i in range(len(spectra)):
        spec = spectra[i]
        new_points = get_points(spec)

        for point in new_points:
            if point[3] < smallest_im:
                smallest_im = point[3]
            if point[3] > largest_im:
                largest_im = point[3]

    return smallest_im, largest_im

def setup_bins(spectra):
    """Sets up the global bins using the smallest and largest ion mobility values.

    Args:
        spectra (list<MSSpectrum>): A list of OpenMS spectra.
    """
    global first_im, last_im, delta_im, bin_size, bins, exps
    global offset_im, bins2, exps2

    print('Getting binning bounds.....................', end='', flush=True)
    first_im, last_im = get_extrema(spectra)
    print('Done')

    delta_im = last_im - first_im
    bin_size = delta_im / num_bins
    print("  Smallest IM:", first_im)
    print("  Largest IM:", last_im, end='\n\n')

    for i in range(num_bins):
        bins.append([])
        exps.append(ms.MSExperiment())

    offset_im = bin_size / 2.0 + first_im
    
    # Need to take into account the first and last half-bins
    for i in range(num_bins + 1):
        bins2.append([])
        exps2.append(ms.MSExperiment())

def run_ff(exp, type):
    """Runs a feature finder on the given input map.

    Args:
        exp (MSExperiment): The OpenMS experiment to run the feature finder on.
        type (string): The name of the feature finder to run.

    Returns:
        FeatureMap: (OpenMS) Contains the found features from the given experiment.
    """
    ff = ms.FeatureFinder()
    ff.setLogType(ms.LogType.CMD)

    features = ms.FeatureMap()
    seeds = ms.FeatureMap()
    params = ms.FeatureFinder().getParameters(type)

    # Parameters for FeatureFinderCentroided
    params.__setitem__(b'mass_trace:min_spectra', 5)
    params.__setitem__(b'mass_trace:max_missing', 2)
    params.__setitem__(b'seed:min_score', 0.5)
    params.__setitem__(b'feature:min_score', 0.5)
    
    exp.updateRanges()
    ff.run(type, exp, features, params, seeds)

    features.setUniqueIds()
    return features

def within_range(target, x, epsilon):
    """Checks if <x> is within a distance of <epsilon> of <target>; returns True if this
    is the case, False otherwise.
    """
    return True if target - epsilon <= x <= target + epsilon else False

def bin_spectrum(spec, mz_epsilon=0.001):
    """Makes two passes at binning a single spectrum (with the second pass shifting the
    bins by 50%).

    Results are saved in the global arrays.

    Args:
        spec (MSSpectrum): The OpenMS spectrum to bin.
    """
    global bins, exps
    global bins2, exps2

    points = get_points(spec)
    # Sort points by IM ascending
    sorted_points = sorted(points, key=itemgetter(3))

    # To circumvent python's aliasing
    temp_bins = []
    for i in range(num_bins):
        temp_bins.append([])

    temp_bins2 = []
    for i in range(num_bins + 1):
        temp_bins2.append([])

    # Step 1: assign points to bins
    for i in range(len(sorted_points)):
        bin_idx = int((sorted_points[i][3] - first_im) / bin_size)
        if bin_idx >= num_bins:
            bin_idx = num_bins - 1
        # Need to cast the list to list to prevent aliasing
        temp_bins[bin_idx].append(list(sorted_points[i]))

        if sorted_points[i][3] < offset_im:
            temp_bins2[0].append(list(sorted_points[i]))
        else:
            bin_idx = int((sorted_points[i][3] - offset_im) / bin_size) + 1
            if bin_idx > num_bins:
                bin_idx = num_bins
            temp_bins2[bin_idx].append(list(sorted_points[i]))

    # Step 2: for each m/z, sum the intensities (pass 1)
    for i in range(num_bins):
        if len(temp_bins[i]) == 0:
            continue

        temp_bins[i] = sorted(temp_bins[i], key=itemgetter(1))
        mz_start, num_mz, curr_mz = 0, 0, temp_bins[i][0][1]
        running_intensity = 0

        for j in range(len(temp_bins[i])):
            if within_range(curr_mz, temp_bins[i][j][1], mz_epsilon):
            #if (temp_bins[i][j][1] == curr_mz):
                num_mz += 1
                running_intensity += temp_bins[i][j][2]
            else:
                # Reached a new m/z slice; update the previous intensities
                for k in range(mz_start, mz_start + num_mz):
                    temp_bins[i][k][2] = running_intensity

                # Update the current counters
                mz_start, num_mz, curr_mz = j, 1, temp_bins[i][j][1]
                running_intensity = temp_bins[i][j][2]

        # Take care of the last slice (if required)
        if num_mz > 0:
            for k in range(mz_start, mz_start + num_mz):
                temp_bins[i][k][2] = running_intensity

        bins[i].extend(temp_bins[i])

        # Build and add a new spectrum
        transpose = list(zip(*temp_bins[i]))

        new_spec = ms.MSSpectrum()
        new_spec.setRT(spec.getRT())
        new_spec.set_peaks((list(transpose[1]), list(transpose[2])))

        # Add IM data
        fda = ms.FloatDataArray()
        for j in list(transpose[3]):
            fda.push_back(j)
        new_spec.setFloatDataArrays([fda])

        exps[i].addSpectrum(new_spec)

    # Step 3: for each m/z, sum the intensities (pass 2)
    for i in range(num_bins + 1):
        if len(temp_bins2[i]) == 0:
            continue

        temp_bins2[i] = sorted(temp_bins2[i], key=itemgetter(1))
        mz_start, num_mz, curr_mz = 0, 0, temp_bins2[i][0][1]
        running_intensity = 0

        for j in range(len(temp_bins2[i])):
            if within_range(curr_mz, temp_bins2[i][j][1], mz_epsilon):
            #if (temp_bins2[i][j][1] == curr_mz):
                num_mz += 1
                running_intensity += temp_bins2[i][j][2]
            else:
                for k in range(mz_start, mz_start + num_mz):
                    temp_bins2[i][k][2] = running_intensity

                mz_start, num_mz, curr_mz = j, 1, temp_bins2[i][j][1]
                running_intensity = temp_bins2[i][j][2]

        if num_mz > 0:
            for k in range(mz_start, mz_start + num_mz):
                temp_bins2[i][k][2] = running_intensity

        bins2[i].extend(temp_bins2[i])

        # Build and add a new spectrum
        transpose = list(zip(*temp_bins2[i]))

        new_spec = ms.MSSpectrum()
        new_spec.setRT(spec.getRT())
        new_spec.set_peaks((list(transpose[1]), list(transpose[2])))

        fda = ms.FloatDataArray()
        for j in list(transpose[3]):
            fda.push_back(j)
        new_spec.setFloatDataArrays([fda])

        exps2[i].addSpectrum(new_spec)

def combine_spectra(exp1, exp2):
    """Adds all of the spectra in <exp2> to <exp1>.

    Args:
        exp1 (MSExperiment): The OpenMS experiment to be added to.
        exp2 (MSExperiment): The OpenMS experiment providing spectra.
    """
    spectra = exp2.getSpectra()
    for spec in spectra:
        exp1.addSpectrum(spec)

# Is there a better way to do this?
def similar_features(feature1, feature2, rt_threshold=5, mz_threshold=0.01):
    """Determines if two features are "similar"; i.e. both their RTs and M/Zs are within a
    certain threshold of each other.

    Args:
        feature1 (Feature or list<float>): An OpenMS feature or a list of floats
            representing a feature.
        feature2 (Feature or list<float>): An OpenMS feature or a list of floats
            representing a feature.
        rt_threshold (float): The maximum threshold in the RT dimension.
        mz_threshold (float): The maximum threshold in the m/z dimension.

    Returns:
        bool: True iff feature1 and feature2 are similar.
    """
    if isinstance(feature1, ms.Feature) and isinstance(feature2, ms.Feature):
        return (abs(feature1.getRT() - feature2.getRT()) < rt_threshold and
                abs(feature1.getMZ() - feature2.getMZ()) < mz_threshold)
    elif isinstance(feature1, list) and isinstance(feature2, list):
        return (abs(feature1[0] - feature2[0]) < rt_threshold and
                abs(feature1[1] - feature2[1]) < mz_threshold)
    else:
        return False

def hull_area(hull):
    """Calculates the area of a convex hull (as a polygon) using the shoelace formula.

    Args:
        hull (libcpp_vector[DPosition2]): A list of points of the convex hull.

    Returns:
        float: The area of the convex hull.
    """
    area = 0.0
    for i in range(len(hull)):
        area += hull[i][0] * hull[(i + 1) % len(hull)][1]
        area -= hull[i][1] * hull[(i + 1) % len(hull)][0]
    return abs(area) / 2.0

def bb_area(box):
    """Calculates the area of a convex hull's bounding box.

    Args:
        box (DBoundingBox2): An OpenMS bounding box for a convex hull.

    Returns:
        float: The area of the bounding box.
    """
    _min = box.minPosition()
    _max = box.maxPosition()
    return abs(_min[0] - _max[0]) * abs(_min[1] * _max[1])

def match_features_internal(features, rt_threshold=0.1, mz_threshold=0.01):
    """Matches overlapping features within a single bin. Intended to correct satellite
    features (a tight group of multiple features becomes a single feature: the one in the
    group with the largest convex hull).

    Args:
        features (FeatureMap): An OpenMS feature map representing a single bin.
        rt_threshold (float): The threshold in the RT dimension.
        mz_threshold (float): The threshold in the m/z dimension.

    Returns:
        FeatureMap: An OpenMS feature map containing all of the uniquely found features.
    """
    matched = ms.FeatureMap()

    for i in range(features.size()):
        f1 = features[i]
        similar = []

        max_area = hull_area(f1.getConvexHull().getHullPoints())
        max_feature = f1

        # Start above f1 to prevent double-matching
        for j in range(i + 1, features.size()):
            f2 = features[j]
            if similar_features(f1, f2, rt_threshold, mz_threshold):
                similar.append(f2)

        for f2 in similar:
            hp = hull_area(f2.getConvexHull().getHullPoints())
            if hp > max_area:
                max_area = hp
                max_feature = f2

        if max_feature not in matched:
            matched.push_back(max_feature)

    return matched

def match_features(features1, features2, rt_threshold=5, mz_threshold=0.01):
    """Matches overlapping features from adjacent bins to each other.
    
    For example, in the first pass of binning, a feature may not be contained entirely
    within its bin. Thus, it must be mapped to itself in an adjacent bin produced by the
    second pass.

    Args:
        features1 (list<FeatureMap>): A list of OpenMS feature maps from the first pass.
        features2 (list<FeatureMap>): A list of OpenMS feature maps from the second pass.
        rt_threshold (float): The threshold in the RT dimension.
        mz_threshold (float): The threshold in the m/z dimension.

    Returns:
        FeatureMap: An OpenMS feature map containing all of the uniquely found features.
    """
    # One bin was used, so there's no need to match anything
    if len(features1) == 1:
        return features1[0]
    
    features = ms.FeatureMap()

    for i in range(len(features1)):
        for f1 in features1[i]:
            similar = []
            max_area = hull_area(f1.getConvexHull().getHullPoints())
            max_feature = f1

            # Should test to see if this gets rid of satellite features
            for f2 in features2[i]:
                if similar_features(f1, f2, rt_threshold, mz_threshold):
                    similar.append(f2)

            for f2 in similar:
                hp = hull_area(f2.getConvexHull().getHullPoints())
                if hp > max_area:
                    max_area = hp
                    max_feature = f2

            if max_feature not in features:
                features.push_back(max_feature)
            # No need to map to the right bin if a match was found in the left
            if len(similar) > 0:
                continue

            similar = []
            max_area = hull_area(f1.getConvexHull().getHullPoints())
            max_feature = f1

            for f2 in features2[i + 1]:
                if similar_features(f1, f2, rt_threshold, mz_threshold):
                    similar.append(f2)

            for f2 in similar:
                hp = hull_area(f2.getConvexHull().getHullPoints())
                if hp > max_area:
                    max_area = hp
                    max_feature = f2

            if len(similar) > 0 and max_feature not in features:
                features.push_back(max_feature)

    return features

def has_peaks(exp):
    """Checks if any of the spectra in an experiment have peaks.

    Args:
        exp (MSExperiment): An OpenMS experiment containing the spectra to check.

    Returns:
        bool: True if any of the spectra have peaks; False otherwise.
    """
    spectra = exp.getSpectra()
    for spec in spectra:
        if spec.size() > 0:
            return True

    return False

def find_features(outdir, outfile, ff_type='centroided', pick=0, imatch=0, min_req=1,
                  window_size=float('Inf'), strict=True, sequential=True):
    """Runs an existing OpenMS feature finder on each of the experiment bins and writes
    the found features to files. Each bin (for each pass) gets its own featureXML and
    mzML files, each pass gets combined files, and the overall experiment gets a matched
    featureXML file from both passes.

    Args:
        outdir (string): The directory to write files to. It must already exist.
        outfile (string): The name of this experiment.
        ff_type (string): The name of the OpenMS feature finder to use.
        pick (boolean): Determines whether or not to pick peak the data before running
            the feature finder.

        min_req (int): see peak_picker.py
        window_size (float): see peak_picker.py
        strict (boolean): see peak_picker.py
        sequential (boolean): see peak_picker.py
    """
    pp = ms.PeakPickerHiRes()
    total_exp = [ms.MSExperiment(), ms.MSExperiment()]
    features = [[], []]
    total_features = [ms.FeatureMap(), ms.FeatureMap()]

    for i in range(num_bins):
        new_exp = ms.MSExperiment()

        if pick == 1:
            pp.pickExperiment(exps[i], new_exp)
        elif pick == 2:
            ms.MzMLFile().store(outdir + '/' + outfile + '-pass1-bin' + str(i) +
                                '-prepick.mzML', exps[i])
            new_exp = peak_picker.peak_pick(exps[i], min_req, window_size, 0.1, strict,
                                            sequential)
        else:
            new_exp = exps[i]

        ms.MzMLFile().store(outdir + '/' + outfile + '-pass1-bin' + str(i) + '.mzML',
                            new_exp)

        temp_features = ms.FeatureMap()
        if has_peaks(new_exp):
            temp_features = run_ff(new_exp, ff_type)
        # Added internal matching here
        if imatch == 1:
            temp_features = match_features_internal(temp_features)

        temp_features.setUniqueIds()
        ms.FeatureXMLFile().store(outdir + '/' + outfile + '-pass1-bin' + str(i) +
                                  '.featureXML', temp_features)

        combine_spectra(total_exp[0], new_exp)
        features[0].append(temp_features)
        total_features[0] += temp_features

    if num_bins == 1:
        return

    # Second pass
    for i in range(num_bins + 1):
        new_exp = ms.MSExperiment()

        if pick == 1:
            pp.pickExperiment(exps2[i], new_exp)
        elif pick == 2:
            ms.MzMLFile().store(outdir + '/' + outfile + '-pass2-bin' + str(i) +
                                '-prepick.mzML', exps[i])
            new_exp = peak_picker.peak_pick(exps2[i], min_req, window_size, 0.1, strict,
                                            sequential)
        else:
            new_exp = exps2[i]

        ms.MzMLFile().store(outdir + '/' + outfile + '-pass2-bin' + str(i) + '.mzML',
                            new_exp)

        temp_features = ms.FeatureMap()
        if has_peaks(new_exp):
            temp_features = run_ff(new_exp, ff_type)
        if imatch == 1:
            temp_features = match_features_internal(temp_features)

        temp_features.setUniqueIds()
        ms.FeatureXMLFile().store(outdir + '/' + outfile + '-pass2-bin' + str(i) +
                                  '.featureXML', temp_features)

        combine_spectra(total_exp[1], new_exp)
        features[1].append(temp_features)
        total_features[1] += temp_features

    # Combine spectra
    ms.MzMLFile().store(outdir + '/' + outfile + '-pass1.mzML', total_exp[0])
    ms.MzMLFile().store(outdir + '/' + outfile + '-pass2.mzML', total_exp[1])

    # Combine features
    total_features[0].setUniqueIds()
    total_features[1].setUniqueIds()
    ms.FeatureXMLFile().store(outdir + '/' + outfile + '-pass1.featureXML', total_features[0])
    ms.FeatureXMLFile().store(outdir + '/' + outfile + '-pass2.featureXML', total_features[1])

    print("Starting cross-bin feature matching. This may take a while.")
    matched_features = match_features(features[0], features[1])
    matched_features.setUniqueIds()
    ms.FeatureXMLFile().store(outdir + '/' + outfile + '.featureXML', matched_features)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='4D LC-IMS/MS Feature Finder.')
    parser.add_argument('--infile', action='store', required=True, type=str)
    parser.add_argument('--outfile', action='store', required=True, type=str)
    parser.add_argument('--outdir', action='store', required=True, type=str)

    parser.add_argument('--num_bins', action='store', required=True, type=int)
    parser.add_argument('--mz_eps', action='store', required=True, type=float)
    parser.add_argument('--int_match', action='store', required=True, type=int)
    parser.add_argument('--peak_pick', action='store', required=True, type=int)

    parser.add_argument('--min_req', action='store', required=True, type=int)
    parser.add_argument('--window_size', action='store', required=True, type=float)
    parser.add_argument('--strict', action='store', required=True, type=int)
    parser.add_argument('--sequential', action='store', required=True, type=int)
    parser.add_argument('--match_only', action='store', required=False)

    args = parser.parse_args()

    num_bins = args.num_bins
    do_matching = True if args.match_only is not None else False

    window_size = args.window_size
    if window_size == -1:
        window_size = float('Inf')

    strict = True if args.strict == 1 else False
    sequential = True if args.sequential == 1 else False
    
    if not do_matching:
        exp = ms.MSExperiment()
        print('Loading mzML input file....................', end='', flush=True)
        ms.MzMLFile().load(args.infile + '.mzML', exp)
        print('Done')

        spectra = exp.getSpectra()
        setup_bins(spectra)

        for i in range(len(spectra)):
            spec = spectra[i]
            # Currently only process MS1 spectra
            if (spec.getMSLevel() != 1):
                continue

            print('Binning MS', end='')
            print(spec.getMSLevel(), 'RT', spec.getRT())
            bin_spectrum(spec, args.mz_eps)

        find_features(args.outdir, args.outfile, 'centroided', args.peak_pick,
                      args.int_match, args.min_req, window_size, strict, sequential)
    else:
        print('Starting feature matcher')
        fm1, fm2 = [], []

        for i in range(args.num_bins):
            fm1.append(ms.FeatureMap())
            ms.FeatureXMLFile().load(args.outdir + '/' + args.outfile + '-pass1-bin' +
                                   str(i) + '.featureXML', fm1[i])

        for i in range(args.num_bins + 1):
            fm2.append(ms.FeatureMap())
            ms.FeatureXMLFile().load(args.outdir + '/' + args.outfile + '-pass2-bin' +
                                   str(i) + '.featureXML', fm2[i])

        f = match_features(fm1, fm2, 10, 0.5)
        f.setUniqueIds()
        ms.FeatureXMLFile().store(args.outdir + '/' + args.outfile + '-m.featureXML', f)
