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
    ms.MzMLFile().load(args.input + '.mzML', exp)

    # Note: the peaks must be sorted to ascending m/z!
    spectra = exp.getSpectra()
    for spec in spectra:
        spec.sortByPosition()
        if not spec.isSorted():
            print("A spectrum is not sorted")
            quit()

    # Parameters, default [type, constraints, default value]
    signal_to_noise = 0.0 # float, 0+, 0
    ms_levels = [ 1 ]
    win_len = 150.0 # float, 1+, 200
    bin_count = 25 # int, 3+, 30
    min_required_elements = 8 # int, 1+, 10

    # Parameters, advanced (not currently iterating through these)
    spacing_difference_gap = 4.0 # float, 0+, 4
    spacing_difference = 1.5 # float, 0+, 1.5
    missing = 1 # int, 0+, 1

    params = ms.PeakPickerHiRes().getParameters()

    params.__setitem__(b'spacing_difference_gap', spacing_difference_gap)
    params.__setitem__(b'spacing_difference', spacing_difference)
    params.__setitem__(b'missing', missing)

    # Iterate through the params and create a new mzML file for each
    # Idea: run FeatureFinderCentroided after each peak-picking and delete the file if
    # no features/a small number of features is found
    params.__setitem__(b'ms_levels', ms_levels)
    for stn in range(signal_to_noise, 1.5, 0.5): # 3 loops
        params.__setitem__(b'signal_to_noise', stn)

        for wl in range(win_len, 300.0, 50.0): # 3 loops
            params.__setitem__(b'SignalToNoise:win_len', wl)

            for bc in range(bin_count, 40, 5): # 3 loops
                params.__setitem__(b'SignalToNoise:bin_count', bc)

                for mre in range(min_required_elements, 13, 1): # 5 loops
                    params.__setitem__(b'SignalToNoise:min_required_elements', mre)

                    pp = ms.PeakPickerHiRes(params)
                    new_exp = ms.MSExperiment()
                    pp.pickExperiment(exp, new_exp)
                    ms.MzMLFile().store(args.output + '/' + 'stn' + str(stn) + 'wl' +
                                        str(wl) + 'bc' + str(bc) + 'mre' + str(mre) +
                                        '.mzML', new_exp)
