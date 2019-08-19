import argparse

import pyopenms as ms
import numpy as np

def peak_pick(exp):
    """A custom peak picker for use with im_binning, since PeakPickerHiRes always
    destroys the data. The idea is to get rid of satellite peaks so that matching
    features within a bin is not required.

    1. Peaks are sorted.
    2. A boolean array is created with a False for every existing peak.
    3. Iterate forward and find a peak with the local highest intensity.
        - Need to check further forward (more than just the next peak)
    4. Go left and right (within a window) until the current peak is less than 5-10% of
        the most intense local peak.
    5. Mark all of the peaks as True.
    6. Create a new peak with its intensity being the sum of the intensities of the newly
        marked peaks and its m/z being a weighted average.

    Args:
        exp (MSExperiment): The OpenMS experiment to be peak picked.

    Returns:
        MSExperiment: An OpenMS experiment corresponding to the input, but peak-picked.
    """
    new_exp = ms.MSExperiment()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Custom peak picker.')
    parser.add_argument('--in', action='store', required=True, type=str)
    parser.add_argument('--out', action='store', required=True, type=str)

    args = parser.parse_args()

    exp = ms.MSExperiment()
    print('Loading mzML input file....................', end='', flush=True)
    ms.MzMLFile().load(args.infile + '.mzML', exp)
    print('Done')
