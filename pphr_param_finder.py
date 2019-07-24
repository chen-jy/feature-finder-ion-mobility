import argparse

import pyopenms as ms
import numpy as np

from im_binning import run_ff

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='PeakPickerHiRes brute-force parameter \
        finder.')
    parser.add_argument('--input', action='store', required=True, type=str)
    parser.add_argument('--output', action='store', required=True, type=str)
    parser.add_argument('--target', action='store', required=False, type=int, default=250)

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
    win_len = 100.0 # float, 1+, 200
    bin_count = 20 # int, 3+, 30
    min_required_elements = 5 # int, 1+, 10

    # Parameters, advanced (not currently iterating through these)
    spacing_difference_gap = 4.0 # float, 0+, 4
    spacing_difference = 1.5 # float, 0+, 1.5
    missing = 1 # int, 0+, 1

    params = ms.PeakPickerHiRes().getParameters()
    
    params.__setitem__(b'ms_levels', ms_levels)

    params.__setitem__(b'spacing_difference_gap', spacing_difference_gap)
    params.__setitem__(b'spacing_difference', spacing_difference)
    params.__setitem__(b'missing', missing)

    # Iterate through the params, and for each set, run PeakPickerHiRes and
    # FeatureFinderCentroided, and save the param sets for those that find at least the
    # target number of features
    file = open("suitable_params.txt", "w+")

    for stn in range(signal_to_noise, 1.5, 0.25):
        params.__setitem__(b'signal_to_noise', stn)

        for wl in range(win_len, 301.0, 25.0):
            params.__setitem__(b'SignalToNoise:win_len', wl)

            for bc in range(bin_count, 41, 1):
                params.__setitem__(b'SignalToNoise:bin_count', bc)

                for mre in range(min_required_elements, 15, 1):
                    params.__setitem__(b'SignalToNoise:min_required_elements', mre)

                    pp = ms.PeakPickerHiRes(params)
                    new_exp = ms.MSExperiment()

                    pp.pickExperiment(exp, new_exp)
                    features = run_ff(new_exp, 'centroided')

                    if features.size() >= args.target:
                        file.write("Number of features: " + str(features.size()) + "\n")
                        file.write("signal_to_noise:    " + str(stn) + "\n")
                        file.write("win_len:            " + str(wl) + "\n")
                        file.write("bin_count:          " + str(bc) + "\n")
                        file.write("min_req_elements:   " + str(mre) + "\n\n")
                        
    file.close()
