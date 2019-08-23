import argparse
import pyopenms as ms
import numpy as np

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Extract MS1 spectra from an mzML file.')
    parser.add_argument('--input', action='store', required=True, type=str)
    parser.add_argument('--output', action='store', required=True, type=str)

    args = parser.parse_args()

    exp = ms.MSExperiment()
    ms.MzMLFile().load(args.input + '.mzML', exp)

    spectra = exp.getSpectra()
    new_spectra = []

    for spec in spectra:
        if spec.getMSLevel() == 1:
             new_spectra.append(spec)

    exp.setSpectra(new_spectra)
    ms.MzMLFile().store(args.output + '.mzML', exp)
