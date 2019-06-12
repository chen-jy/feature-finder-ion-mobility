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
        interior lists holds RT, MZ, intensity, and IM information (in that order)
        for a single peak in the spectrum. The exterior list is unsorted.
    """
    point_data = zip(*spec.get_peaks(), spec.getFloatDataArrays()[0])
    return [[spec.getRT(), mz, intensity, im] for mz, intensity, im in point_data]

def find_features(args):
    """The primary point of execution for the experimental feature finder.
    """
    exp = ms.MSExperiment()
    ms.MzMLFile().load(args.infile + '.mzML', exp)

    # Store the RT, MZ, intensity, and IM data for every peak in every spectrum
    point_cloud = []

    spectra = exp.getSpectra()
    for i in range(args.num_frames):
        spec = spectra[i]

        new_points = get_points(spec)
        point_cloud.extend(new_points)

    # Sort points by IM ascending (using lambda significantly slower)
    start = time.time()
    #sorted_cloud = sorted(point_cloud, key=lambda x: x[3])
    sorted_cloud = sorted(point_cloud, key=itemgetter(3))
    end = time.time()
    print("Length of data:", len(sorted_cloud))
    print("Time to sort:", end - start)

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
    find_features(args)
