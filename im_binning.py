from baseline import *

import time
from operator import itemgetter

# Globals and constants
bins, exps = [], []
first_im, last_im, delta_im = 0, 0, 0
num_bins, bin_size = 50, 0

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

def find_features_old(spec, outdir, outfile, spec_idx=0):
    """Make one pass at binning a spectrum and finding its features.

    Args:
        spec (MSSpectrum): An OpenMS MSSpectrum object.
        outdir (string): The output directory for FeatureXML files.
        outfile (string): A string to identify this series of runs.
        spec_idx (int): The index of this spectrum in a series of spectra.
    """
    points = get_points(spec)
    # Sort points by IM ascending
    sorted_points = sorted(points, key=itemgetter(3))

    # Position of bin i (0-indexed) = i * bin_size + first_im
    first_im, last_im = sorted_points[0][3], sorted_points[len(points) - 1][3]
    delta_im, offset_im = last_im - first_im, 0
    
    num_bins = 50
    bin_size = delta_im / num_bins
    # For successive passes, use offset_im to shift bins
    offset_delta, pass_num = 0.05, 1

    bins, new_exp = [], []
    for i in range(num_bins):
        bins.append([])
        new_exp.append(ms.MSExperiment())

    # Step 1: assign points to bins
    for i in range(len(points)):
        # Need to adapt this formula for offset_im
        bin_idx = int((sorted_points[i][3] - first_im) / bin_size)
        if bin_idx >= num_bins:
            bin_idx = num_bins - 1
        bins[bin_idx].append(sorted_points[i])

    # Step 2: for each m/z, average the intensities
    for i in range(num_bins):
        if len(bins[i]) == 0:
            continue

        bins[i] = sorted(bins[i], key=itemgetter(1))
        mz_start, num_mz, curr_mz = 0, 0, bins[i][0][1]
        run_intensity = 0

        for j in range(len(bins[i])):
            if (bins[i][j][1] == curr_mz):
                num_mz += 1
                run_intensity += bins[i][j][2]
            else:
                # Reached a new m/z slice; update the previous intensities
                run_intensity /= num_mz
                for k in range(mz_start, mz_start + num_mz):
                    bins[i][k][2] = run_intensity

                mz_start, num_mz, curr_mz = j, 1, bins[i][j][1]
                run_intensity = bins[i][j][2]

        # Takes care of the last slice (if required)
        if num_mz > 0:
            run_intensity /= num_mz
            for k in range(mz_start, mz_start + num_mz):
                bins[i][k][2] = run_intensity

        # Get the arrays of RT, MZ, Intensity, and IM
        transpose = list(zip(*bins[i]))

        new_spec = ms.MSSpectrum()
        new_spec.setRT(spec.getRT())
        new_spec.set_peaks((list(transpose[1]), list(transpose[2])))

        fda = ms.FloatDataArray()
        for j in list(transpose[3]):
            fda.push_back(j)
        new_spec.setFloatDataArrays([fda])

        new_exp[i].addSpectrum(new_spec)
        ms.MzMLFile().store(outdir + '/' + outfile + '-spec' + str(spec_idx) + '-pass' +
                            str(pass_num) + '-bin' + str(i) + '.mzML', new_exp[i])

    # Step 3: find the features for each bin
    ff_type = 'centroided'
    for i in range(num_bins):
        features = run_ff(new_exp[i], ff_type)
        ms.FeatureXMLFile().store(outdir + '/' + outfile + '-spec' + str(spec_idx) +
                                  '-pass' + str(pass_num) + '-bin' + str(i) +
                                  '.featureXML', features)

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

def find_features(outdir, outfile, ff_type='centroided'):
    #pp = ms.PeakPickerHiRes()

    for i in range(num_bins):
        #new_exp = ms.MSExperiment()
        #pp.pickExperiment(exps[i], new_exp)
        ms.MzMLFile().store(outdir + '/' + outfile + '-pass1-bin' + str(i) + '.mzML',
                            exps[i])

        features = run_ff(exps[i], ff_type)
        ms.FeatureXMLFile().store(outdir + '/' + outfile + '-pass1-bin' + str(i) +
                                  '.featureXML', features)

    # Second pass
    for i in range(num_bins + 1):
        #new_exp = ms.MSExperiment()
        #pp.pickExperiment(exps2[i], new_exp)
        ms.MzMLFile().store(outdir + '/' + outfile + '-pass2-bin' + str(i) + '.mzML',
                            exps2[i])

        features = run_ff(exps2[i], ff_type)
        ms.FeatureXMLFile().store(outdir + '/' + outfile + '-pass2-bin' + str(i) +
                                  '.featureXML', features)

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
