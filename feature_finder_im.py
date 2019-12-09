"""LC-IMS-MS/MS feature finder.

Extends OpenMS (via pyOpenMS) feature finding capabilities to work on 4D LC-IMS-MS/MS data.
"""

import argparse
from operator import itemgetter
import os
from typing import Any, List, Optional, Tuple

import pyopenms as ms

import peak_picker


class FeatureFinderIonMobility:
    """The LC-IMS-MS/MS feature finder.

    No attributes are public, and the only methods that should be called are run() and reset().
    """

    MZ_EPSILON = 0.001
    MIN_INTENSITY = 0.1
    RT_THRESHOLD = 5.0
    MZ_THRESHOLD = 0.01

    def __init__(self) -> None:
        self.reset()

    def reset(self) -> None:
        """Resets the feature finder to its default state."""
        self.bins = [[], []]  # bins[0] is the first pass and bins[1] is the second
        self.exps = [[], []]
        self.num_bins, self.bin_size = 0, 0
        self.im_start, self.im_end = 0, 0
        self.im_delta, self.im_offset = 0, 0

    def get_spectrum_points(self, spec: ms.MSSpectrum) -> List[List[float]]:
        """Extracts the retention times, mass to charges, intensities, and ion mobility values of
        all peaks in a spectrum.

        Keyword arguments:
        spec: the spectrum to extract data points from

        Returns: a list of lists, where each interior list holds RT, m/z, intensity, and IM data
        (in that order) for a single peak in the spectrum.
        """
        point_data = zip(*spec.get_peaks(), spec.getFloatDataArrays()[0])
        return [[spec.getRT(), mz, intensity, im] for mz, intensity, im in point_data]

    def get_im_extrema(self, spectra: List[ms.MSSpectrum]) -> Tuple[float, float]:
        """Finds the smallest and largest IM values in a list of spectra.

        Keyword arguments:
        spectra: the list of spectra with IM data to scan through

        Returns: a tuple of the smallest and largest IM values, in that order.
        """
        smallest_im, largest_im = float('inf'), float('-inf')

        for spec in spectra:
            points = self.get_spectrum_points(spec)
            for point in points:
                if point[3] < smallest_im:
                    smallest_im = point[3]
                if point[3] > largest_im:
                    largest_im = point[3]

        return smallest_im, largest_im

    def setup_bins(self, spectra: List[ms.MSSpectrum]) -> None:
        """Sets up the IM bins for feature finding.

        Keyword arguments:
        spectra: the list of spectra to bin
        """
        print('Getting IM bounds', flush=True)
        self.im_start, self.im_end = self.get_im_extrema(spectra)

        self.im_delta = self.im_end - self.im_start
        self.bin_size = self.im_delta / self.num_bins
        self.im_offset = self.bin_size / 2.0 + self.im_start

        for i in range(self.num_bins):
            for j in range(2):
                self.bins[j].append([])
                self.exps[j].append(ms.MSExperiment())
        self.bins[1].append([])  # The last half-bin in the second pass
        self.exps[1].append(ms.MSExperiment())

    def within_epsilon(self, target: float, var: float) -> bool:
        """Checks if var is within the m/z epsilon of target."""
        return target - self.MZ_EPSILON <= var <= target + self.MZ_EPSILON

    def bin_spectrum(self, spec: ms.MSSpectrum) -> None:
        """Bins a single spectrum in two passes.

        Results are stored in the private bin lists.

        Keyword arguments:
        spec: the spectrum to bin
        """
        points = self.get_spectrum_points(spec)
        points = sorted(points, key=itemgetter(3))  # Ascending IM

        temp_bins = [[], [[]]]  # New bins to prevent aliasing
        new_bins = [[], [[]]]
        for i in range(self.num_bins):
            for j in range(2):
                temp_bins[j].append([])
                new_bins[j].append([])

        for i in range(len(points)):  # Assign points to bins
            bin_idx = int((points[i][3] - self.im_start) / self.bin_size)
            if bin_idx >= self.num_bins:
                bin_idx = self.num_bins - 1
            temp_bins[0][bin_idx].append(list(points[i]))  # Cast list to list to prevent aliasing

            bin_idx = int((points[i][3] - self.im_offset) / self.bin_size) + 1
            if points[i][3] < self.im_offset:
                bin_idx = 0
            elif bin_idx > self.num_bins:
                bin_idx = self.num_bins
            temp_bins[1][bin_idx].append(list(points[i]))

        for i in range(self.num_bins):  # First pass
            if len(temp_bins[0][i]) == 0:
                continue

            temp_bins[0][i] = sorted(temp_bins[0][i], key=itemgetter(1))  # Ascending m/z
            mz_start, num_mz, curr_mz = 0, 0, temp_bins[0][i][0][1]
            running_intensity = 0

            for j in range(len(temp_bins[0][i])):
                if self.within_epsilon(curr_mz, temp_bins[0][i][j][1]):
                    num_mz += 1
                    running_intensity += temp_bins[0][i][j][2]
                else:  # Reached a new m/z slice
                    temp_bins[0][i][mz_start][2] = running_intensity
                    new_bins[0][i].append(list(temp_bins[0][i][mz_start]))
                    mz_start, num_mz, curr_mz = j, 1, temp_bins[0][i][j][1]
                    running_intensity = temp_bins[0][i][j][2]

            temp_bins[0][i][mz_start][2] = running_intensity  # Take care of the last slice
            new_bins[0][i].append(list(temp_bins[0][i][mz_start]))
            self.bins[0][i].extend(new_bins[0][i])

            transpose = list(zip(*new_bins[0][i]))
            new_spec = ms.MSSpectrum()  # The final binned spectrum
            new_spec.setRT(spec.getRT())
            new_spec.set_peaks((list(transpose[1]), list(transpose[2])))
            self.exps[0][i].addSpectrum(new_spec)

        for i in range(self.num_bins + 1):  # Second pass
            if len(temp_bins[1][i]) == 0:
                continue

            temp_bins[1][i] = sorted(temp_bins[1][i], key=itemgetter(1))
            mz_start, num_mz, curr_mz = 0, 0, temp_bins[1][i][0][1]
            running_intensity = 0

            for j in range(len(temp_bins[1][i])):
                if self.within_epsilon(curr_mz, temp_bins[1][i][j][1]):
                    num_mz += 1
                    running_intensity += temp_bins[1][i][j][2]
                else:
                    temp_bins[1][i][mz_start][2] = running_intensity
                    new_bins[1][i].append(list(temp_bins[1][i][mz_start]))
                    mz_start, num_mz, curr_mz = j, 1, temp_bins[1][i][j][1]
                    running_intensity = temp_bins[1][i][j][2]

            temp_bins[1][i][mz_start][2] = running_intensity
            new_bins[1][i].append(list(temp_bins[1][i][mz_start]))
            self.bins[1][i].extend(new_bins[1][i])

            transpose = list(zip(*new_bins[1][i]))
            new_spec = ms.MSSpectrum()
            new_spec.setRT(spec.getRT())
            new_spec.set_peaks((list(transpose[1]), list(transpose[2])))
            self.exps[1][i].addSpectrum(new_spec)

    def similar_features(self, feature1: Any, feature2: Any) -> bool:
        """Checks if the RTs and m/zs of two features are within fixed thresholds of each other."""
        if isinstance(feature1, ms.Feature) and isinstance(feature2, ms.Feature):
            return (abs(feature1.getRT() - feature2.getRT()) < self.RT_THRESHOLD and
                    abs(feature1.getMZ() - feature2.getMZ()) < self.MZ_THRESHOLD)
        elif isinstance(feature1, list) and isinstance(feature2, list):
            return (abs(feature1[0] - feature2[0]) < self.RT_THRESHOLD and
                    abs(feature1[1] - feature2[1]) < self.MZ_THRESHOLD)
        elif isinstance(feature1, ms.Feature) and isinstance(feature2, list):
            return (abs(feature1.getRT() - feature2[0]) < self.RT_THRESHOLD and
                    abs(feature1.getMZ() - feature2[1]) < self.MZ_THRESHOLD)
        elif isinstance(feature1, list) and isinstance(feature2, ms.Feature):
            return self.similar_features(feature2, feature1)
        else:
            return False

    def polygon_area(self, polygon: List[Tuple[float, float]]) -> float:
        """Computes the area of a convex polygon using the shoelace formula."""
        area = 0.0
        for i in range(len(polygon)):
            area += polygon[i][0] * polygon[(i + 1) % len(polygon)][1]
            area -= polygon[i][1] * polygon[(i + 1) % len(polygon)][0]
        return abs(area) / 2.0

    def match_features_internal(self, features: ms.FeatureMap) -> ms.FeatureMap:
        """Matches features in a single bin; intended to correct satellite features.

        The feature in each feature set with the largest convex hull becomes the 'representative'
        feature of that set.

        Keyword arguments:
        features: the features of a single bin for intra-bin matching

        Returns: a matched set of features.
        """
        matched = ms.FeatureMap()
        for i in range(features.size()):
            feature1 = features[i]
            similar = []
            max_area = self.polygon_area(feature1.getConvexHull().getHullPoints())
            max_feature = feature1

            for j in range(i + 1, features.size()):  # Start above to prevent double matching
                feature2 = features[j]
                if self.similar_features(feature1, feature2):
                    similar.append(feature2)
            for feature2 in similar:
                area = self.polygon_area(feature2.getConvexHull().getHullPoints())
                if area > max_area:
                    max_area = area
                    max_feature = feature2

            if max_feature not in matched:
                matched.push_back(max_feature)

        return matched

    def match_features(self, features1: List[ms.FeatureMap], features2: List[ms.FeatureMap]) -> ms.FeatureMap:
        """Matches features from adjacent bins across two passes.

        The feature in each feature set with the largest convex hull becomes the 'representative'
        feature of that set.

        Keyword arguments:
        features1: the features of the first pass
        features2: the features of the second pass

        Returns: a matched set of features.
        """
        if len(features1) == 1:
            return features1[0]
    
        matched = ms.FeatureMap()

        for i in range(len(features1)):  # Match the first pass against the second
            for feature1 in features1[i]:
                similar = []
                max_area = self.polygon_area(feature1.getConvexHull().getHullPoints())
                max_feature = feature1

                for feature2 in features2[i]:  # The half-overlapping bin to the left
                    if self.similar_features(feature1, feature2):
                        similar.append(feature2)
                for feature2 in similar:
                    area = self.polygon_area(feature2.getConvexHull().getHullPoints())
                    if area > max_area:
                        max_area = area
                        max_feature = feature2

                if max_feature not in matched:
                    matched.push_back(max_feature)
                if len(similar) > 0:
                    continue

                similar = []
                max_area = self.polygon_area(feature1.getConvexHull().getHullPoints())
                max_feature = feature1

                for feature2 in features2[i + 1]:  # The half-overlapping bin to the right
                    if self.similar_features(feature1, feature2):
                        similar.append(feature2)
                for feature2 in similar:
                    area = self.polygon_area(feature2.getConvexHull().getHullPoints())
                    if area > max_area:
                        max_area = area
                        max_feature = feature2

                # If len(similar) == 0, then max_feature has already been added
                if len(similar) > 0 and max_feature not in matched:
                    matched.push_back(max_feature)

        for i in range(len(features2)):  # Match the second pass against the first
            for feature1 in features2[i]:
                similar = []
                max_area = self.polygon_area(feature1.getConvexHull().getHullPoints())
                max_feature = feature1

                if i > 0:
                    for feature2 in features1[i - 1]:
                        if self.similar_features(feature1, feature2):
                            similar.append(feature2)
                    for feature2 in similar:
                        area = self.polygon_area(feature2.getConvexHull().getHullPoints())
                        if area > max_area:
                            max_area = area
                            max_feature = feature2

                    is_new = True  # Avoid duplicate features from the first matching pass
                    for feature3 in matched:
                        if self.similar_features(max_feature, feature3):
                            is_new = False
                            break

                    if is_new:
                        matched.push_back(max_feature)
                    if len(similar) > 0:
                        continue

                similar = []
                max_area = self.polygon_area(feature1.getConvexHull().getHullPoints())
                max_feature = feature1

                if i < len(features1):
                    for feature2 in features1[i]:
                        if self.similar_features(feature1, feature2):
                            similar.append(feature2)
                    for feature2 in similar:
                        area = self.polygon_area(feature2.getConvexHull().getHullPoints())
                        if area > max_area:
                            max_area = area
                            max_feature = feature2

                    is_new = True
                    for feature3 in matched:
                        if self.similar_features(max_feature, feature3):
                            is_new = False
                            break

                    if len(similar) > 0 and is_new:
                        matched.push_back(max_feature)

        return matched

    def run_ff(self, exp: ms.MSExperiment, type: str) -> ms.FeatureMap:
        """Runs an existing OpenMS feature finder on an experiment.

        Keyword arguments:
        exp: the experiment to run the existing feature finder on
        type: the name of the existing feature finder to run ('centroided')

        Returns: the features in the experiment.
        """
        ff = ms.FeatureFinder()
        ff.setLogType(ms.LogType.NONE)
        features, seeds = ms.FeatureMap(), ms.FeatureMap()

        # TODO: tighten up these parameters
        params = ms.FeatureFinder().getParameters(type)
        params.__setitem__(b'mass_trace:min_spectra', 5)
        params.__setitem__(b'mass_trace:max_missing', 2)
        params.__setitem__(b'seed:min_score', 0.5)
        params.__setitem__(b'feature:min_score', 0.5)
    
        exp.updateRanges()
        ff.run(type, exp, features, params, seeds)

        features.setUniqueIds()
        return features

    def has_peaks(self, exp: ms.MSExperiment) -> bool:
        """Checks if any spectrum in an experiment has peaks."""
        spectra = exp.getSpectra()
        for spec in spectra:
            if spec.size() > 0:
                return True
        return False

    def find_features(self, dir: str, peak_pick: str, peak_radius: int, window_size: float, cus_mode: int,
                      ff_type: str, suppress: bool) -> List[List[ms.FeatureMap]]:
        """Runs optional peak picking and an existing feature finder on each IM bin.

        Keyword arguments:
        dif: the directory to write the intermediate output files to
        peak_pick: the peak picker to use ('none', 'pphr', or 'cus')
        peak_radius: for the custom peak picker, the minimum peak radius of a peak set
        window_size: for the custom peak picker, the maximum m/z window size to consider
        cus_mode: for the custom peak picker, the mode to use (0 for m/z or 1 for intensity)
        ff_type: the existing feature finder to use ('centroided')
        suppress: determine if intermediate output files should be written

        Returns: a list of two lists (for the passes), each containing the features for all of
        their bins.
        """
        features = [[], []]
        total_features = [ms.FeatureMap(), ms.FeatureMap()]
        pp = ms.PeakPickerHiRes()

        for i in range(self.num_bins):  # First pass
            new_exp = ms.MSExperiment()
            if peak_pick != 'none' and not suppress:
                ms.MzMLFile().store(dir + '/pass1-bin' + str(i) + '-prepick.mzML', self.exps[0][i])

            if peak_pick == 'pphr':
                pp.pickExperiment(self.exps[0][i], new_exp)
            elif peak_pick == 'cus':
                # TODO: replace this line after refactoring peak_picker.py
                new_exp = peak_picker.peak_pick(self.exps[0][i], peak_radius, window_size, self.MIN_INTENSITY, True,
                                                False)
            else:
                new_exp = self.exps[0][i]

            # TODO: maybe call the prepicked file '/pass1-bin${i}' and add '-picked' to this file.
            if not suppress:
                ms.MzMLFile().store(dir + '/pass1-bin' + str(i) + '.mzML', new_exp)

            temp_features = ms.FeatureMap()
            if self.has_peaks(new_exp):
                temp_features = self.run_ff(new_exp, ff_type)

            temp_features = self.match_features_internal(temp_features)
            temp_features.setUniqueIds()
            if not suppress:
                ms.FeatureXMLFile().store(dir + '/pass1-bin' + str(i) + '.featureXML', temp_features)

            features[0].append(temp_features)
            total_features[0] += temp_features

        if self.num_bins == 1:
            return

        # TODO: this for loop is essentially the same as above; maybe merge them?
        for i in range(self.num_bins + 1):  # Second pass
            new_exp = ms.MSExperiment()
            if peak_pick != 'none' and not suppress:
                ms.MzMLFile().store(dir + '/pass2-bin' + str(i) + '-prepick.mzML', self.exps[1][i])

            if peak_pick == 'pphr':
                pp.pickExperiment(self.exps[1][i], new_exp)
            elif peak_pick == 'cus':
                # TODO: replace this line after refactoring peak_picker.py
                new_exp = peak_picker.peak_pick(self.exps[1][i], peak_radius, window_size, self.MIN_INTENSITY, True,
                                                False)
            else:
                new_exp = self.exps[1][i]

            if not suppress:
                ms.MzMLFile().store(dir + '/pass2-bin' + str(i) + '.mzML', new_exp)

            temp_features = ms.FeatureMap()
            if self.has_peaks(new_exp):
                temp_features = self.run_ff(new_exp, ff_type)

            temp_features = self.match_features_internal(temp_features)
            temp_features.setUniqueIds()
            if not suppress:
                ms.FeatureXMLFile().store(dir + '/pass2-bin' + str(i) + '.featureXML', temp_features)

            features[1].append(temp_features)
            total_features[1] += temp_features

        if not suppress:
            total_features[0].setUniqueIds()
            ms.FeatureXMLFile().store(dir + '/pass1.featureXML', total_features[0])
            total_features[1].setUniqueIds()
            ms.FeatureXMLFile().store(dir + '/pass2.featureXML', total_features[1])

        return features[0], features[1]

    def run(self, exp: ms.MSExperiment(), dir: str, num_bins: int, peak_pick: str, peak_radius: int,
            window_size: float, cus_mode: int, ff_type: str, suppress: bool) -> ms.FeatureMap:
        """Runs the feature finder on an experiment.

        Keyword arguments:
        exp: the experiment to run the feature finder on
        dir: the directory to write the intermediate output files to
        num_bins: the number of IM bins to use
        peak_pick: the peak picker to use ('none', 'pphr', or 'cus')
        peak_radius: for the custom peak picker, the minimum peak radius of a peak set
        window_size: for the custom peak picker, the maximum m/z window size to consider
        cus_mode: for the custom peak picker, the mode to use (0 for m/z or 1 for intensity)
        ff_type: the existing feature finder to use ('centroided')
        suppress: determine if intermediate output files should be written

        Returns: the features found by the feature finder.
        """
        self.reset()
        self.num_bins = num_bins

        spectra = exp.getSpectra()
        self.setup_bins(spectra)

        print('Starting binning')
        for spec in spectra:
            if spec.getMSLevel() != 1:  # Currently only works on MS1 scans
                continue
            self.bin_spectrum(spec)

        print('Starting feature finding')
        features1, features2 = self.find_features(dir, peak_pick, peak_radius, window_size, cus_mode, ff_type,
                                                  suppress)

        print('Starting feature matching')
        matched = self.match_features(features1, features2)
        matched.setUniqueIds()

        return matched


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='LC-IMS-MS/MS Feature Finder.')

    parser.add_argument('-i', '--in', action='store', required=True, type=str, dest='in_',
                        help='the input mzML file')
    parser.add_argument('-o', '--out', action='store', required=True, type=str,
                        help='the output featureXML file')
    parser.add_argument('-d', '--dir', action='store', required=False, type=str, default='.',
                        help='the output directory')
    parser.add_argument('-n', '--num_bins', action='store', required=True, type=int,
                        help='the number of IM bins to use')
    parser.add_argument('-f', '--finder', action='store', required=False, type=str, default='centroided',
                        choices=['centroided'], help='the existing feature finder to use')

    parser.add_argument('-p', '--peak_pick', action='store', required=False, type=str, default='none',
                        choices=['none', 'pphr', 'cus'], help='the peak picker to use')
    parser.add_argument('-r', '--peak_radius', action='store', required=False, type=int, default=1,
                        help='the peak radius for the custom peak picker')
    parser.add_argument('-w', '--window_size', action='store', required=False, type=float, default=0.015,
                        help='the window size for the custom peak picker')
    parser.add_argument('-m', '--cus_mode', action='store', required=False, type=int, default=0,
                        choices=[0, 1], help='the mode of the custom peak picker')
    parser.add_argument('-s', '--suppress', action='store_true', required=False, default=False,
                        help='suppress the writing of intermediate mzML and featureXML files')

    args = parser.parse_args()

    if not os.path.isfile(args.in_):
        print('Error:', args.in_, 'is not a file')
        exit(1)
    if not os.path.isdir(args.dir):
        print('Error:', args.dir, 'is not an existing directory')
        exit(1)
    if not args.out.endswith('.featureXML'):
        print('Error:', args.out, 'must be a featureXML file')
        exit(1)
    
    exp = ms.MSExperiment()
    print('Loading mzML input file', flush=True)
    ms.MzMLFile().load(args.in_, exp)

    ff = FeatureFinderIonMobility()
    # TODO: maybe make a parameter class so that run_ff doesn't need so many arguments
    features = ff.run(exp, args.dir, args.num_bins, args.peak_pick, args.peak_radius, args.window_size,
                         args.cus_mode, args.finder, args.suppress)
    ms.FeatureXMLFile().store(args.dir + '/' + args.out, features)
    print('Found', features.size(), 'features')
