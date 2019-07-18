import argparse

import pyopenms as ms
import numpy as np

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='PeakPickerHiRes brute-force parameter \
        finder.')
    parser.add_argument('--input', action='store', required=True, type=str)
    parser.add_argument('--output', action='store', required=True, type=str)

    args = parser.parse_args()

    exp = ms.MSExperiment()
    print('Loading mzML input file....................', end='', flush=True)
    ms.MzMLFile().load(args.input + '.mzML', exp)
    print('Done')

    # Note: the peaks must be sorted to ascending m/z!
    spectra = exp.getSpectra()
    for spec in spectra:
        spec.sortByPosition()
        if not spec.isSorted():
            print("A spectrum is not sorted")
            quit()

    # Parameters, default
    signal_to_noise = 0 # float, 0+, 0
    ms_levels = [ 1 ]
    win_len = 200 # float, 1+, 200
    bin_count = 30 # int, 3+, 30
    min_required_elements = 10 # int, 1+, 10

    # Parameters, advanced
    spacing_difference_gap = 0 # float, 0+, 4
    spacing_difference = 0 # float, 0+, 1.5
    missing = 1 # int, 0+, 1

    pp = ms.PeakPickerHiRes()
    params = ms.PeakPickerHiRes().getParameters()

    params.__setitem__(b'signal_to_noise', signal_to_noise)
    params.__setitem__(b'ms_levels', ms_levels)
    params.__setitem__(b'SignalToNoise:win_len', win_len)
    params.__setitem__(b'SignalToNoise:bin_count', bin_count)
    params.__setitem__(b'SignalToNoise:min_required_elements', min_required_elements)

    params.__setitem__(b'spacing_difference_gap', spacing_difference_gap)
    params.__setitem__(b'spacing_difference', spacing_difference)
    params.__setitem__(b'missing', missing)

    # Iterate through the params and create a new mzML file for each
    new_exp = ms.MSExperiment()
    pp.pickExperiment(exp, new_exp)
    ms.MzMLFile().store(args.output + '.mzML', new_exp)
