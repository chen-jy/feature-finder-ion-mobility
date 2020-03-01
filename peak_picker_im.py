"""A custom peak picker for (optional) use with the LC-IMS-MS/MS feature finder.
"""

import argparse
import os

import pyopenms as ms


class PeakPickerIonMobility:
    """A custom peak picker for use with feature_finder_im. This tries to get rid of satellite
    peaks so that erroneous/duplicate features are minimized.

    The algorithm:
    1. Sort peaks by increasing m/z (or decreasing intensity) and mark all peaks as "not picked".
    2. Walk from left to right and check if the current peak is a local maximum.
        - Not necessary in decreasing intensity mode.
    3. For each local maximum, walk left (and right) while peak intensities decrease.
        - If the params are satisfied, construct a "peak set" and mark all peaks as "picked".
    4. Replace the peak set with a single peak with:
        - intensity = the sum of all the intensities of the peaks in the peak set.
        - position = the intensity-weighted average of all of the positions in the peak set.

    There are no attributes and only two (public) methods.
    """

    def __init__(self) -> None:
        pass

    def pick_spectra(self, spec: ms.MSSpectrum, peak_radius: int = 1, window_radius: float = 0.015,
                     pp_mode: str = 'int', min_int_mult: float = 0.10, strict: bool = True) -> ms.MSSpectrum():
        """Peak picks a single spectrum.

        Keyword arguments:
        spec: the spectrum to peak pick
        peak_radius: the minimum peak radius of a peak set
        window_radius: the maximum m/z window radius of a peak set
        pp_mode: the mode to use ('ltr' or 'int')
        min_int_mult: a multiplier to the maximum peak intensity in a set (for differentiating
            between signal and noise)
        strict: if False, allow a single increase in intensity in either direction

        Returns: the peak picked spectrum.
        """
        num_peaks = spec.size()
        spec.sortByPosition()

        peak_idx = []  # Intensity lookup table
        picked = [False] * num_peaks

        if pp_mode == 'int':
            for i in range(num_peaks):
                peak_idx.append([spec[i].getIntensity(), i])
            peak_idx = sorted(peak_idx, reverse=True)

        picked_spec = ms.MSSpectrum()
        picked_spec.setMSLevel(1)
        picked_spec.setRT(spec.getRT())

        for idx in range(num_peaks):  # Begin peak picking
            i = idx if pp_mode == 'ltr' else peak_idx[idx][1]
            if picked[i]:
                continue

            init_intensity = spec[i].getIntensity()
            total_intensity = spec[i].getIntensity()
            init_position = spec[i].getPos()
            left_picked, right_picked = 0, 0
            low_bound, high_bound = i, i

            sFlag = False  # Flag for when strict is False
            threshold = peak_radius

            for j in range(i - 1, -1, -1):  # Walk left
                if picked[j] or abs(spec[j].getPos() - init_position) > window_radius:
                    break

                if spec[j].getIntensity() > spec[j + 1].getIntensity():
                    if strict or sFlag or j + 1 == i:  # Don't start with an abnormal peak
                        break
                    sFlag = True
                    threshold += 1  # End the peak set with a lower peak ("increase peak_radius")

                total_intensity += spec[j].getIntensity()
                left_picked += 1
                low_bound -= 1

                if left_picked >= threshold and spec[j].getIntensity() <= init_intensity * min_int_mult:
                    break

            if left_picked < threshold:
                continue
            sFlag = False
            threshold = peak_radius

            for j in range(i + 1, num_peaks):  # Walk right
                if picked[j] or abs(spec[j].getPos() - init_position) > window_radius:
                    break

                if spec[j].getIntensity() > spec[j - 1].getIntensity():
                    if strict or sFlag or j - 1 == i:
                        break
                    sFlag = True
                    threshold += 1

                total_intensity += spec[j].getIntensity()
                right_picked += 1
                high_bound += 1

                if right_picked >= threshold and spec[j].getIntensity() <= init_intensity * min_int_mult:
                    break

            if right_picked < threshold:
                continue

            total_position = 0
            for j in range(low_bound, high_bound + 1):
                picked[j] = True
                if total_intensity != 0:
                    total_position += spec[j].getPos() * (spec[j].getIntensity() / total_intensity)

            p = ms.Peak1D()
            p.setIntensity(total_intensity)
            p.setPos(total_position)
            picked_spec.push_back(p)

        return picked_spec

    def pick_experiment(self, exp: ms.MSExperiment, peak_radius: int = 1, window_radius: float = 0.015,
            pp_mode: str = 'int', min_int_mult: float = 0.10, strict: bool = True) -> ms.MSExperiment():
        """Peak picks an experiment.

        Keyword arguments:
        exp: the experiment to peak pick
        peak_radius: the minimum peak radius of a peak set
        window_radius: the maximum m/z window radius of a peak set
        pp_mode: the mode to use ('ltr' or 'int')
        min_int_mult: a multiplier to the maximum peak intensity in a set (for differentiating
            between signal and noise)
        strict: if False, allow a single increase in intensity in either direction

        Returns: the peak picked experiment.
        """
        exp.sortSpectra()
        spectra = exp.getSpectra()

        picked_exp = ms.MSExperiment()
        for spec in spectra:
            if spec.getMSLevel() != 1:
                continue
            picked_exp.addSpectrum(self.pick_spectra(spec, peak_radius, window_radius, pp_mode, min_int_mult, strict))

        return picked_exp


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Custom LC-IMS-MS/MS peak picker.')
    parser.add_argument('-i', '--in', action='store', required=True, type=str, dest='in_',
                        help='the input mzML file')
    parser.add_argument('-o', '--out', action='store', required=True, type=str,
                        help='the output mzML file')
    parser.add_argument('-r', '--peak_radius', action='store', required=False, type=int, default=1,
                        help='the minimum peak radius of a peak set')
    parser.add_argument('-w', '--window_radius', action='store', required=False, type=float, default=0.015,
                        help='the maximum window radius of a peak set')
    parser.add_argument('-m', '--pp_mode', action='store', required=False, type=str, default='int',
                        choices=['ltr', 'int'], help='the peak picking mode to use')
    args = parser.parse_args()

    if not os.path.isfile(args.in_):
        print('Error:', args.in_, 'is not a file')
        exit(1)
    if not args.out.endswith('.mzML'):
        print('Error:', args.out, 'must be an mzML file')
        exit(1)

    exp = ms.MSExperiment()
    print('Loading mzML input file', flush=True)
    ms.MzMLFile().load(args.in_, exp)

    pp = PeakPickerIonMobility()
    print('Running peak picker', flush=True)
    picked_exp = pp.run()

    print('Writing mzML output file', flush=True)
    ms.MzMLFile().store(args.output, picked_exp)
