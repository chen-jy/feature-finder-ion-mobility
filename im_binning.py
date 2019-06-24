from baseline import *

import time
from operator import itemgetter

# Globals and constants
bins, exps = [], []
first_im, last_im, delta_im = 0, 0, 0
num_bins, bin_size = 10, 0

# For the second pass (shift the bins by 50%)
bins2, exps2 = [], []
offset_im = 0

def get_points(spec):
    """Data preprocessing to extract the retention time, mass to charge, intensity,
    and ion mobility for each peak in a spectrum.

    Args:
        spec (MSSpectrum): An OpenMS MSSpectrum object.

    Returns:
        list<list<float, float, float, float>>: A list of lists, where each interior
        list holds RT, MZ, intensity, and IM information (in that order) for a single
        peak in the spectrum. The exterior list is unsorted.
    """
    point_data = zip(*spec.get_peaks(), spec.getFloatDataArrays()[0])
    return [[spec.getRT(), mz, intensity, im] for mz, intensity, im in point_data]

def get_points_pp(spec_pp, spec):
    """Does the same thing as get_points(), but for a spectrum that has been peak-
    picked (losing its IM information in the process).

    Args:
        spec_pp (MSSpectrum): An OpenMS MSSpectrum object that has been peak-picked.
        spec (MSSpectrum): An OpenMS MSSpectrum object, not peak-picked, corresponding
            to spec_pp.

    Returns:
        list<list<float, float, float, float>>: A list of lists, where each interior
        list holds RT, MZ, intensity, and IM information (in that order) for a single
        peak in the spectrum. The exterior list is unsorted.
    """
    point_data = zip(*spec_pp.get_peaks(), spec.getFloatDataArrays()[0])
    return [[spec.getRT(), mz, intensity, im] for mz, intensity, im in point_data]

def get_extrema(spectra):
    """Finds the smallest and largest IM values in an experiment.

    Args:
        spectra (list<MSSpectrum>): A list of OpenMS MSSpectrum objects.

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
    """Sets up the global bins using the smallest and greatest ion mobility values.

    Args:
        spectra (list<MSSpectrum>): A list of OpenMS MSSpectrum objects.
    """
    global first_im, last_im, delta_im, bin_size, bins, exps
    global offset_im, bins2, exps2

    print('Getting binning bounds.....................', end='', flush=True)
    first_im, last_im = get_extrema(spectra)
    # Approximation for debugging
    #first_im, last_im = 0.5, 1.7
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
        exp (MSExperiment): An OpenMS MSExperiment object.
        type (string): The name of the feature finder to run.

    Returns:
        FeatureMap: Contains the found features from the given experiment.
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

def bin_spectrum(spec, outdir, outfile):
    """Makes a single pass at binning a single spectrum. Needs to eventually support
    an overlapping series of bins (multiple passes).

    Results are saved in the global array <bins>.

    Args:
        spec (MSSpectrum): An OpenMS MSSpectrum object.
        outdir (string): The output directory for writing FeatureXML files.
        outfile (string): An identifier for this series of runs.
    """
    global bins, exps
    global bins2, exps2

    #pp = ms.PeakPickerHiRes()
    #spec_pp = ms.MSSpectrum()
    #pp.pick(spec, spec_pp)
    #points = get_points_pp(spec_pp, spec)

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
            if (temp_bins[i][j][1] == curr_mz):
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

        # Step 2.5: build and add a new spectrum
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
            if (temp_bins2[i][j][1] == curr_mz):
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

        # Step 3.5: build and add a new spectrum
        transpose = list(zip(*temp_bins2[i]))

        new_spec = ms.MSSpectrum()
        new_spec.setRT(spec.getRT())
        new_spec.set_peaks((list(transpose[1]), list(transpose[2])))

        fda = ms.FloatDataArray()
        for j in list(transpose[3]):
            fda.push_back(j)
        new_spec.setFloatDataArrays([fda])

        exps2[i].addSpectrum(new_spec)

def similar_features(feature1, feature2):
    rt_threshold = 0.1
    mz_threshold = 0.1
    intensity_threshold = 0.1

    return (abs(feature1.getRT() - feature2.getRT()) < rt_threshold and
            abs(feature1.getMZ() - feature2.getMZ()) < mz_threshold and
            abs(feature1.getIntensity() - feature2.getIntensity()) < intensity_threshold)

def match_features(features1, features2):
    if len(features1) == 1:
        return features1[0]
    
    features = []
    for i in range(len(features1)):
        for j in range(len(features1[i])):
            f1 = features1[i][j]

            for k in range(len(features2[i])):
                f2 = features2[i][k]

                if similar_features(f1, f2):
                    pass
                elif f1 not in features:
                    features.append(f1)

            for k in range(len(features2[i + 1])):
                f2 = features2[i + 1][k]

                if similar_features(f1, f2):
                    pass
                elif f1 not in features:
                    features.append(f1)

    return features

def find_features(outdir, outfile, ff_type='centroided'):
    """Runs an existing feature finder on the experiment bins and writes the features to
    files. Each bin (for each pass) gets its own mzML and featureXML files, each pass
    gets a combined featureXML file, and the overall experiment gets a matched featureXML
    file (from both passes).

    Args:
        outdir (string): The targeted output directory. It must already exist.
        outfile (string): The name of this experiment.
        ff_type (string): The name of the existing feature finder to use. Defaults to
            centroided for peak-picked data.
    """

    #pp = ms.PeakPickerHiRes()
    features = [[], []]
    total_features = [[], []]

    for i in range(num_bins):
        #new_exp = ms.MSExperiment()
        #pp.pickExperiment(exps[i], new_exp)
        ms.MzMLFile().store(outdir + '/' + outfile + '-pass1-bin' + str(i) + '.mzML',
                            exps[i])

        temp_features = run_ff(exps[i], ff_type)
        ms.FeatureXMLFile().store(outdir + '/' + outfile + '-pass1-bin' + str(i) +
                                  '.featureXML', temp_features)

        features[0].append(temp_features)
        total_features[0].extend(temp_features)

    # Second pass
    for i in range(num_bins + 1):
        #new_exp = ms.MSExperiment()
        #pp.pickExperiment(exps2[i], new_exp)
        ms.MzMLFile().store(outdir + '/' + outfile + '-pass2-bin' + str(i) + '.mzML',
                            exps2[i])

        temp_features = run_ff(exps2[i], ff_type)
        ms.FeatureXMLFile().store(outdir + '/' + outfile + '-pass2-bin' + str(i) +
                                  '.featureXML', temp_features)

        features[1].append(temp_features)
        total_features[1].extend(temp_features)

    # Combine features
    ms.FeatureXMLFile().store(outdir + '/' + outfile + '-pass1.mzML', total_features[0])
    ms.FeatureXMLFile().store(outdir + '/' + outfile + '-pass2.mzML', total_features[1])

    matched_features = match_features(features[0], features[1])
    ms.FeatureXMLFile().store(outdir + '/' + outfile + '.mzML', matched_features)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='4D LC-IMS/MS Feature Finder.')
    parser.add_argument('--infile', action='store', required=True, type=str)
    parser.add_argument('--outfile', action='store', required=True, type=str)
    parser.add_argument('--outdir', action='store', required=True, type=str)

    args = parser.parse_args()
    
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

        print("Processing", spec.getMSLevel(), "RT", spec.getRT())
        bin_spectrum(spec, args.outdir, args.outfile)

    find_features(args.outdir, args.outfile, 'centroided')
