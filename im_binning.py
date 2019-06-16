from baseline import *

import time
from operator import itemgetter

def get_points(spec):
    """Data preprocessing to extract the retention time, mass to charge, intensity,
    and ion mobility for each peak in a spectrum.

    Args:
        spec (MSSpectrum): An OpenMS MSSpectrum object.

    Returns:
        list<list<double, double, double, double>>: A list of lists, where each
        interior list holds RT, MZ, intensity, and IM information (in that order)
        for a single peak in the spectrum. The exterior list is unsorted.
    """
    point_data = zip(*spec.get_peaks(), spec.getFloatDataArrays()[0])
    return [[spec.getRT(), mz, intensity, im] for mz, intensity, im in point_data]

@deprecated
def run_old(args):
    """Collects all points from all spectra.
    """
    exp = ms.MSExperiment()
    ms.MzMLFile().load(args.infile + '.mzML', exp)
    print("Raw data file loaded; beginning execution")

    # Store the RT, MZ, intensity, and IM data for every peak in every spectrum
    point_cloud = []

    spectra = exp.getSpectra()
    for i in range(args.num_frames):
        spec = spectra[i]

        new_points = get_points(spec)
        point_cloud.extend(new_points)

    print("Sorting data points")
    # Sort points by IM ascending (using lambda significantly slower)
    start = time.time()
    #sorted_cloud = sorted(point_cloud, key=lambda x: x[3])
    sorted_cloud = sorted(point_cloud, key=itemgetter(3))
    end = time.time()
    print("Number of data points:", len(sorted_cloud))
    print("Time to sort:", end - start)

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
    
    exp.update_ranges()
    ff.run(name, input_map, features, params, seeds)

    features.setUniqueIds()
    return features

def find_features(spec):
    """Make one pass at binning a spectrum and finding its features.

    Args:
        spec (MSSpectrum): An OpenMS MSSpectrum object.
    """
    points = get_points(spec)
    # Sort points by IM ascending
    sorted_points = sorted(points, key=itemgetter(3))

    # Position of bin i (0-indexed) = i * bin_size + first_im
    first_im, last_im = sorted_points[0][3], sorted_points[len(points) - 1][3]
    delta_im = last_im - first_im
    
    num_bins = 50
    bin_size = delta_im / num_bins
    bins = [[]] * num_bins

    new_exp = []
    for i in range(num_bins):
        new_exp.append(ms.MSExperiment())

    # Step 1: assign points to bins
    for i in range(len(points)):
        bin_idx = int((sorted_points[i][3] - first_im) / delta_im)
        bins[bin_idx].append(sorted_points[i])

    # Step 2: for each m/z, average the intensities
    for i in range(num_bins):
        pass

    # Step 3: find the features for each bin
    for i in range(num_bins):
        pass

if __name__ == "__main__":
    # Includes legacy arguments from baseline.py
    parser = argparse.ArgumentParser(description='FragToPre Clustering Baseline')
    parser.add_argument('--infile', action='store', required=True, type=str)
    parser.add_argument('--outfile', action='store', required=True, type=str)
    parser.add_argument('--outdir', action='store', required=True, type=str)
    parser.add_argument('--mz_epsilon', action='store', required=False, type=float)
    parser.add_argument('--im_epsilon', action='store', required=False, type=float)
    parser.add_argument('--num_frames', action='store', required=False, type=int)
    parser.add_argument('--window_size', action='store', required=False, type=int)
    parser.add_argument('--rt_length', action='store', required=False, type=int)

    args = parser.parse_args()

    exp = ms.MSExperiment()
    ms.MzMLFile().load(args.infile + '.mzML', exp)

    spectra = exp.getSpectra()
    for i in range(args.num_frames):
        spec = spectra[i]
        find_features(spec)
