import argparse

import pyopenms as ms
import numpy as np

# Try window_size=0.02
def peak_pick(exp, maxima_threshold=3, window_size=float(inf), strict=True):
    """A custom peak picker for use with im_binning, since PeakPickerHiRes always
    destroys the data. The idea is to get rid of satellite peaks so that matching
    features within a bin is not required.

    1. Peaks are sorted by m/z.
    2. A boolean array is created with a False for every existing peak.
    3. Iterate forward and find a peak with the local highest intensity.
        - Need to check further forward (more than just the next peak)
    4. Go left and right (within a window) until the current peak is less than 10% of the
        most intense local peak.
    5. Mark all of the peaks as True.
    6. Create a new peak with its intensity being the sum of the intensities of the newly
        marked peaks and its m/z being a weighted average.

    Args:
        exp (MSExperiment): The OpenMS experiment to be peak picked.
        maxima_threshold (int): The required number of decreasing/non-increasing peaks to
            either side of a peak in order to be considered a local maximum.
        window_size (float): The maximum m/z distance left/right from the initial peak to
            consider.
        strict (bool): If true, peaks must be non-increasing from the initial peak.
            Otherwise, a single peak is allowed to break this rule.

    Returns:
        MSExperiment: An OpenMS experiment corresponding to the input, but peak-picked.
    """
    exp.sortSpectra()
    spectra = exp.getSpectra()
    new_exp = ms.MSExperiment()

    for spec in spectra:
        if spec.getMSLevel() != 1:
            continue
        print("Peak picking RT", spec.getRT())

        spec.sortByPosition()
        num_peaks = spec.size()

        new_spec = ms.MSSpectrum()
        new_spec.setMSLevel(1)
        new_spec.setRT(spec.getRT())

        picked = [False] * num_peaks

        for i in range(0, num_peaks):
            if picked[i]:
                continue

            init_intensity = spec[i].getIntensity()
            total_intensity = spec[i].getIntensity()
            init_position = spec[i].getPos()
            total_position = spec[i].getPos()
            left_picked, right_picked = 0, 0
            low_bound, high_bound = i, i

            # Walk left
            # Flag for a single increase in intensity (when strict=False)
            sflag = False
            for j in range(i - 1, -1, -1):
                if abs(spec[j].getPos() - init_position) > window_size:
                    break

                if spec[j].getIntensity() > spec[j + 1].getIntensity():
                    if strict or sflag:
                        break
                    sflag = True

                total_intensity += spec[j].getIntensity()
                total_position += spec[j].getPos()
                left_picked += 1
                low_bound -= 1

                if left_picked >= maxima_threshold and \
                    spec[j].getIntensity() <= init_intensity * 0.1:
                    break

            if left_picked < maxima_threshold:
                continue

            # Walk right
            sflag = False
            for j in range(i + 1, num_peaks):
                if abs(spec[j].getPos() - init_position) > window_size:
                    break

                if spec[j].getIntensity() > spec[j - 1].getIntensity():
                    if strict or sflag:
                        break
                    sflag = True

                total_intensity += spec[j].getIntensity()
                total_position += spec[j].getPos()
                right_picked += 1
                high_bound += 1

                if right_picked >= maxima_threshold and \
                    spec[j].getIntensity() <= init_intensity * 0.1:
                    break

            if right_picked < maxima_threshold:
                continue

            for j in range(low_bound, high_bound + 1):
                picked[j] = True

            p = ms.Peak1D()
            p.setIntensity(total_intensity)
            p.setPos(total_position / (left_picked + right_picked + 1))
            new_spec.push_back(p)

        new_exp.addSpectrum(new_spec)

    return new_exp

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Custom peak picker.')
    parser.add_argument('--input', action='store', required=True, type=str)
    parser.add_argument('--output', action='store', required=True, type=str)

    args = parser.parse_args()

    exp = ms.MSExperiment()
    print('Loading mzML input file....................', end='', flush=True)
    ms.MzMLFile().load(args.input + '.mzML', exp)
    print('Done')

    new_exp = peak_pick(exp)
    ms.MzMLFile().store(args.output + '.mzML', new_exp)
