"""LC-IMS-MS/MS feature finder.

Extends OpenMS (via pyOpenMS) feature finding capabilities to work on 4D LC-IMS-MS/MS data.
"""

import argparse
import csv
from operator import itemgetter
import os
import psutil
import time
from typing import List, Optional, Tuple

import pyopenms as ms

import common_utils_im as util
import peak_picker_im as ppim


class FeatureFinderIonMobility:
    """The LC-IMS-MS/MS feature finder.

    There are no public attributes, and the only public method is run().

    TODO: don't load everything into memory at once; maybe cache to disk
    TODO: make a separate parameter class so that run() doesn't require so many arguments.
    """

    MZ_EPSILON = 0.001  # For binning
    MIN_INTENSITY = 0.1  # For the custom peak picker
    RT_THRESHOLD = 5.0  # For feature matching
    MZ_THRESHOLD = 0.01  # For feature matching

    def __init__(self) -> None:
        self.reset()

    def reset(self) -> None:
        """Resets the feature finder to its default state."""
        self.bins = [[], []]  # bins[0] is the first pass and bins[1] is the second
        self.exps = [[], []]  # Same as above
        self.num_bins, self.bin_size = 0, 0
        self.im_start, self.im_end = 0, 0
        self.im_delta, self.im_offset = 0, 0
        self.im_scan_nums = [[], []]  # Keep the midpoint IM value for each bin

    def setup_bins(self, spectra: List[ms.MSSpectrum]) -> None:
        """Sets up the IM bins for feature finding.

        Keyword arguments:
        spectra: the list of spectra to bin
        """
        print('Getting IM bounds.', end=' ', flush=True)
        self.im_start, self.im_end = util.get_im_extrema(spectra)

        self.im_delta = self.im_end - self.im_start
        self.bin_size = self.im_delta / self.num_bins
        self.im_offset = self.im_start + self.bin_size / 2.0

        for i in range(self.num_bins):
            for j in range(2):
                self.bins[j].append([])
                self.exps[j].append(ms.MSExperiment())
            bin_start, bin_stop = i * self.bin_size + self.im_start, (i + 1) * self.bin_size + self.im_start
            self.im_scan_nums[0].append((bin_start + bin_stop) / 2.0)

        self.im_scan_nums[1].append((self.im_start + self.im_offset) / 2.0)
        for i in range(1, self.num_bins):
            bin_start, bin_stop = i * self.bin_size + self.im_offset, (i + 1) * self.bin_size + self.im_offset
            self.im_scan_nums[1].append((bin_start + bin_stop) / 2.0)
        self.bins[1].append([])  # The last half-bin in the second pass
        self.exps[1].append(ms.MSExperiment())

        print('Done', flush=True)
        print('\tSmallest IM:', self.im_start)
        print('\tLargest IM:', self.im_end, end='\n\n')

    def within_epsilon(self, target: float, var: float) -> bool:
        """Checks if var is within the m/z epsilon of target."""
        return target - self.MZ_EPSILON <= var <= target + self.MZ_EPSILON

    def bin_spectrum(self, spec: ms.MSSpectrum) -> None:
        """Bins a single spectrum in two passes.

        Results are stored in the private bin lists.

        Keyword arguments:
        spec: the spectrum to bin
        """
        points = util.get_spectrum_points(spec)
        points = sorted(points, key=itemgetter(3))  # Ascending IM

        temp_bins = [[], [[]]]  # New bins for each pass to prevent aliasing
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
            mz_start, curr_mz = 0, temp_bins[0][i][0][1]
            running_intensity = 0

            for j in range(len(temp_bins[0][i])):
                if self.within_epsilon(curr_mz, temp_bins[0][i][j][1]):
                    running_intensity += temp_bins[0][i][j][2]
                else:  # Reached a new m/z slice
                    temp_bins[0][i][mz_start][2] = running_intensity
                    new_bins[0][i].append(list(temp_bins[0][i][mz_start]))
                    mz_start, curr_mz = j, temp_bins[0][i][j][1]
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
            mz_start, curr_mz = 0, temp_bins[1][i][0][1]
            running_intensity = 0

            for j in range(len(temp_bins[1][i])):
                if self.within_epsilon(curr_mz, temp_bins[1][i][j][1]):
                    running_intensity += temp_bins[1][i][j][2]
                else:
                    temp_bins[1][i][mz_start][2] = running_intensity
                    new_bins[1][i].append(list(temp_bins[1][i][mz_start]))
                    mz_start, curr_mz = j, temp_bins[1][i][j][1]
                    running_intensity = temp_bins[1][i][j][2]

            temp_bins[1][i][mz_start][2] = running_intensity
            new_bins[1][i].append(list(temp_bins[1][i][mz_start]))
            self.bins[1][i].extend(new_bins[1][i])

            transpose = list(zip(*new_bins[1][i]))
            new_spec = ms.MSSpectrum()
            new_spec.setRT(spec.getRT())
            new_spec.set_peaks((list(transpose[1]), list(transpose[2])))
            self.exps[1][i].addSpectrum(new_spec)

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
            max_area = util.polygon_area(feature1.getConvexHull().getHullPoints())
            max_feature = feature1

            for j in range(features.size()):
                if i == j:
                    continue
                feature2 = features[j]
                if util.similar_features(feature1, feature2, self.RT_THRESHOLD, self.MZ_THRESHOLD):
                    similar.append(feature2)

            for feature2 in similar:
                area = util.polygon_area(feature2.getConvexHull().getHullPoints())
                if area > max_area:
                    max_area = area
                    max_feature = feature2

            if max_feature not in matched:
                matched.push_back(max_feature)

        return matched

    def match_features_pass(self, features: List[ms.FeatureMap]) -> List[Tuple[ms.Feature, int]]:
        """Matches features in a contiguous sequence of adjacent bins in a single pass. This should
        reduce the amount of redundant features.

        Keyword arguments:
        features: the list of feature maps (one per bin) to match.

        Returns: a list of features and their bin indices for which their intensities are highest.
        """
        used = []
        for bin_idx in range(len(features)):
            features[bin_idx].sortByPosition()  # Ascending m/z
            used.append([False] * features[bin_idx].size())

        matched = []  # All features
        not_unique = []  # Features found in at least two contiguous bins

        for bin_idx in range(len(features)):
            for i in range(features[bin_idx].size()):
                if used[bin_idx][i]:
                    continue

                feature1 = features[bin_idx][i]
                feature_indices = [(bin_idx, i)]
                next_idx = bin_idx + 1

                while next_idx < len(features):
                    similar = []
                    for j in range(features[next_idx].size()):
                        if used[next_idx][j]:
                            continue

                        feature2 = features[next_idx][j]
                        if util.similar_features(feature1, feature2, self.RT_THRESHOLD, self.MZ_THRESHOLD):
                            similar.append((feature2, j))
                            used[next_idx][j] = True

                    if len(similar) == 0:  # Cannot extend the chain any further
                        break

                    max_area = util.polygon_area(similar[0][0].getConvexHull().getHullPoints())
                    max_feature = similar[0]
                    for j in range(1, len(similar)):
                        area = util.polygon_area(similar[j][0].getConvexHull().getHullPoints())
                        if area > max_area:
                            max_area = area
                            max_feature = similar[j]  # (feature, index)

                    feature_indices.append((next_idx, max_feature[1]))
                    next_idx += 1

                max_intensity, max_feature, max_idx = feature1.getIntensity(), feature1, bin_idx
                for j in range(len(feature_indices)):
                    b_idx, f_idx = feature_indices[j]
                    intensity = features[b_idx][f_idx].getIntensity()
                    if intensity > max_intensity:
                        max_intensity = intensity
                        max_feature = features[b_idx][f_idx]
                        max_idx = b_idx

                matched.append((max_feature, max_idx))
                if len(feature_indices) > 1:
                    not_unique.append((max_feature, max_idx))

        return matched

    def match_features(self, features1: List[ms.FeatureMap], features2: List[ms.FeatureMap]) -> \
            Tuple[ms.FeatureMap, List[Tuple[ms.Feature, float]]]:
        """Matches found features across passes to reduce the amount of redundant features.

        Keyword arguments:
        features1: the list of feature maps (one per bin) for the first pass
        features2: the list of feature maps (one per bin) for the second pass

        Returns: the matched feature map, as well as a list of features and their IM values
            (corresponding to the midpoint IM value of the bin that each feature is located in; we
            can't use the exact IM value because they are not computed by FeatureFinderCentroided).
        """
        pass1 = self.match_features_pass(features1)  # List[(ms.Feature, int)]
        pass2 = self.match_features_pass(features2)  # Interior tuples are (feature, bin index)

        used = [False] * len(pass2)
        feature_bins = []  # Holds (feature, IM)

        for (feature1, bin1) in pass1:
            similar = []
            for j in range(len(pass2)):
                if used[j]:
                    continue
                feature2, bin2 = pass2[j][0], pass2[j][1]
                if util.similar_features(feature1, feature2, self.RT_THRESHOLD, self.MZ_THRESHOLD) and \
                    (bin1 == bin2 or bin1 + 1 == bin2):
                    similar.append((feature2, bin2, j))
                    used[j] = True

            max_area = util.polygon_area(feature1.getConvexHull().getHullPoints())
            max_feature = (feature1, self.im_scan_nums[0][bin1])
            max_idx = -1

            for (feature2, bin2, j) in similar:
                area = util.polygon_area(feature2.getConvexHull().getHullPoints())
                if area > max_area:
                    max_area = area
                    max_feature = (feature2, self.im_scan_nums[1][bin2])
                    max_idx = j

            feature_bins.append(max_feature)

        for j in range(len(pass2)):  # Features unique to the second pass
            if not used[j]:
                feature_bins.append((pass2[j][0], self.im_scan_nums[1][pass2[j][1]]))

        cleaned, clean_bins = ms.FeatureMap(), []  # Clean up potential duplicates
        used = [False] * len(feature_bins)

        for i in range(len(feature_bins)):
            if used[i]:
                continue
            used[i] = True

            similar = []
            for j in range(len(feature_bins)):
                if used[j]:
                    continue
                if util.similar_features(feature_bins[i][0], feature_bins[j][0], self.RT_THRESHOLD,
                                          self.MZ_THRESHOLD) and feature_bins[i][1] == feature_bins[j][1]:
                    similar.append(feature_bins[j])
                    used[j] = True

            max_feature, max_area = feature_bins[i], \
                                    util.polygon_area(feature_bins[i][0].getConvexHull().getHullPoints())
            for feature in similar:
                area = util.polygon_area(feature[0].getConvexHull().getHullPoints())
                if area > max_area:
                    max_feature, max_area = feature, area

            cleaned.push_back(max_feature[0])
            clean_bins.append(max_feature)

        return cleaned, clean_bins

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

        params = ms.FeatureFinder().getParameters(type)  # default (Leon's) (modified)
        params.__setitem__(b'mass_trace:min_spectra', 7)  # 10 (5) (7)
        params.__setitem__(b'mass_trace:max_missing', 1)  # 1 (2) (1)
        params.__setitem__(b'seed:min_score', 0.65)  # 0.8 (0.5) (0.65)
        params.__setitem__(b'feature:min_score', 0.6)  # 0.7 (0.5) (0.6)
    
        exp.updateRanges()
        ff.run(type, exp, features, params, seeds)

        features.setUniqueIds()
        return features

    def run_ffm(self, exp: ms.MSExperiment) -> ms.FeatureMap:
        """Runs FeatureFinderMultiplex on an experiment.

        Keyword arguments:
        exp: the experiment to run the feature finder on

        Returns: the features in the experiment.
        """
        ffm = ms.FeatureFinderMultiplexAlgorithm()
        params = ffm.getDefaults()
        params.setValue(b'algorithm:labels', b'[]')
        params.setValue(b'algorithm:spectrum_type', b'centroid')
        ffm.setParameters(params)

        _ = False
        ffm.run(exp, _)

        return ffm.getFeatureMap()

    def find_features(self, pp_type: str, peak_radius: int, window_radius: float, pp_mode: str, ff_type: str,
                      dir: str, filter: str, debug: bool) -> List[List[ms.FeatureMap]]:
        """Runs optional peak picking and then an existing feature finder on each IM bin.

        Keyword arguments:
        pp_type: the peak picker to use ('none', 'pphr', or 'custom')
        peak_radius: for the custom peak picker, the minimum peak radius of a peak set
        window_radius: for the custom peak picker, the maximum m/z window radius to consider
        pp_mode: for the custom peak picker, the mode to use ('ltr' or 'int')
        ff_type: the existing feature finder to use ('centroided')
        dir: the directory to write the intermediate output files to
        filter: the noise filter to use
        debug: determine if intermediate output files should be written

        Returns: a list of two lists (for the passes), each containing the features for all of
        their bins.
        """
        features = [[], []]
        total_features = [ms.FeatureMap(), ms.FeatureMap()]

        if filter == 'gauss':
            filter_g = ms.GaussFilter()
            params_g = filter_g.getDefaults()
            params_g.setValue(b'ppm_tolerance', 20.0)
            params_g.setValue(b'use_ppm_tolerance', b'true')
            filter_g.setParameters(params_g)
        
        if filter == 'sgolay':
            filter_s = ms.SavitzkyGolayFilter()
            params_s = filter_s.getDefaults()
            params_s.setValue(b'frame_length', 7)
            params_s.setValue(b'polynomial_order', 3)
            filter_s.setParameters(params_s)

        pick_hr = ms.PeakPickerHiRes()
        pick_im = ppim.PeakPickerIonMobility()  # Maybe use a parameter class?

        nb = [self.num_bins, 0 if self.num_bins == 1 else self.num_bins + 1]

        for j in range(2):  # Pass index
            for i in range(nb[j]):  # Bin index
                new_exp = ms.MSExperiment()
                if debug:
                    ms.MzMLFile().store(dir + '/pass' + str(j + 1) + '-bin' + str(i) + '.mzML', self.exps[j][i])

                # Optional noise filtering
                if filter == 'gauss':
                    filter_g.filterExperiment(self.exps[j][i])
                elif filter == 'sgolay':
                    filter_s.filterExperiment(self.exps[j][i])

                if filter != 'none' and debug:
                    ms.MzMLFile().store(dir + '/pass' + str(j + 1) + '-bin' + str(i) + '-filtered.mzML',
                                        self.exps[j][i])

                # Optional peak picking
                if pp_type == 'pphr':
                    pick_hr.pickExperiment(self.exps[j][i], new_exp)
                elif pp_type == 'custom':
                    new_exp = pick_im.pick_experiment(self.exps[j][i], peak_radius, window_radius, pp_mode,
                                                      self.MIN_INTENSITY, strict=True)
                else:
                    new_exp = self.exps[j][i]

                if pp_type != 'none' and debug:
                    ms.MzMLFile().store(dir + '/pass' + str(j + 1) + '-bin' + str(i) + '-picked.mzML', new_exp)

                # Feature finding
                temp_features = ms.FeatureMap()
                if util.has_peaks(new_exp):
                    if ff_type == 'centroided':
                        temp_features = self.run_ff(new_exp, ff_type)
                    elif ff_type == 'multiplex':
                        temp_features = self.run_ffm(new_exp)

                temp_features = self.match_features_internal(temp_features)
                temp_features.setUniqueIds()
                if debug:
                    ms.FeatureXMLFile().store(dir + '/pass' + str(j + 1) + '-bin' + str(i) + '.featureXML',
                                              temp_features)

                features[j].append(temp_features)
                total_features[j] += temp_features

        if debug:
            for j in range(2):
                total_features[j].setUniqueIds()
                ms.FeatureXMLFile().store(dir + '/pass' + str(j + 1) + '.featureXML', total_features[j])

        return features[0], features[1]

    def run(self, exp: ms.MSExperiment(), num_bins: int = 50, pp_type: str = 'pphr', peak_radius: int = 1,
            window_radius: float = 0.015, pp_mode: str = 'int', ff_type: str = 'centroided', dir: str = '.',
            filter: str = 'none', debug: bool = False) -> ms.FeatureMap:
        """Runs the feature finder on an experiment.

        Keyword arguments:
        exp: the experiment to run the feature finder on
        num_bins: the number of IM bins to use
        pp_type: the peak picker to use ('none', 'pphr', or 'custom')
        peak_radius: for the custom peak picker, the minimum peak radius of a peak set
        window_radius: for the custom peak picker, the maximum m/z window radius to consider
        pp_mode: for the custom peak picker, the mode to use ('ltr' or 'int')
        ff_type: the existing feature finder to use ('centroided')
        dir: the directory to write the intermediate output files to
        filter: the noise filter to use
        debug: determine if intermediate output files should be written

        Returns: the features found by the feature finder.
        """
        time_out = open(dir + '/timing.log', 'a')
        pymem = psutil.Process(os.getpid())
        mem_use = pymem.memory_info()[0] / 2.0 ** 30
        time_out.write(f'mzml load: {mem_use} GiB\n')
        start_tt = time.time()

        self.reset()
        self.num_bins = num_bins

        start_t = time.time()
        spectra = exp.getSpectra()
        self.setup_bins(spectra)

        total_t = time.time() - start_t
        time_out.write(f'setup bins: {total_t}s\n')
        mem_use = pymem.memory_info()[0] / 2.0 ** 30
        time_out.write(f'setup bins: {mem_use} GiB\n')

        print('Starting binning.', flush=True)
        start_t = time.time()
        for spec in spectra:
            if spec.getMSLevel() != 1:  # Currently only works on MS1 scans
                continue
            print('Binning RT', spec.getRT(), flush=True)
            self.bin_spectrum(spec)
        print('Done')

        total_t = time.time() - start_t
        time_out.write(f'binning: {total_t}s\n')
        mem_use = pymem.memory_info()[0] / 2.0 ** 30
        time_out.write(f'binning: {mem_use} GiB\n')

        print('Starting feature finding.', end=' ', flush=True)
        start_t = time.time()
        features1, features2 = self.find_features(pp_type, peak_radius, window_radius, pp_mode, ff_type, dir, filter,
                                                  debug)
        print('Done')

        total_t = time.time() - start_t
        time_out.write(f'feature finding: {total_t}s\n')
        mem_use = pymem.memory_info()[0] / 2.0 ** 30
        time_out.write(f'feature finding: {mem_use} GiB\n')

        if self.num_bins == 1:  # Matching between passes for one bin results in no features
            features1[0].setUniqueIds()
            return features1[0]

        print('Starting feature matching.', end=' ', flush=True)
        start_t = time.time()
        all_features, feature_bins = self.match_features(features1, features2)
        all_features.setUniqueIds()
        print('Done')

        total_t = time.time() - start_t
        time_out.write(f'feature matching: {total_t}s\n')
        mem_use = pymem.memory_info()[0] / 2.0 ** 30
        time_out.write(f'feature matching: {mem_use} GiB\n')

        indexed_bins = [[f.getRT(), f.getMZ(), bin] for f, bin in feature_bins]
        with open(dir + '/feature-bins.csv', 'w', newline='') as file:
            writer = csv.writer(file)
            writer.writerow(['RT', 'm/z', 'im'])
            writer.writerows(indexed_bins)
            
        total_t = time.time() - start_tt
        time_out.write(f'total: {total_t}s\n')
        time_out.close()

        return all_features


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='LC-IMS-MS/MS Feature Finder.')

    parser.add_argument('-i', '--in', action='store', required=True, type=str, dest='in_',
                        help='the input mzML file')
    parser.add_argument('-o', '--out', action='store', required=True, type=str,
                        help='the output featureXML file')
    parser.add_argument('-d', '--dir', action='store', required=False, type=str, default='.',
                        help='the output directory')
    parser.add_argument('-n', '--num_bins', action='store', required=False, type=int, default=50,
                        help='the number of IM bins to use')

    parser.add_argument('-p', '--pp_type', action='store', required=False, type=str, default='pphr',
                        choices=['none', 'pphr', 'custom'], help='the peak picker to use')
    parser.add_argument('-r', '--peak_radius', action='store', required=False, type=int, default=1,
                        help='the peak radius for the custom peak picker')
    parser.add_argument('-w', '--window_radius', action='store', required=False, type=float, default=0.015,
                        help='the window radius for the custom peak picker')
    parser.add_argument('-m', '--pp_mode', action='store', required=False, type=str, default='int',
                        choices=['ltr', 'int'], help='the mode of the custom peak picker')

    parser.add_argument('-f', '--ff_type', action='store', required=False, type=str, default='centroided',
                        choices=['centroided', 'multiplex'], help='the existing feature finder to use')
    parser.add_argument('-e', '--filter', action='store', required=False, type=str, default='none',
                        choices=['none', 'gauss', 'sgolay'], help='the noise filter to use')

    parser.add_argument('--debug', action='store_true', required=False, default=False,
                        help='write intermediate mzML and featureXML files')

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

    time_out = open(args.dir + '/timing.log', 'w')
    start_t = time.time()
    
    exp = ms.MSExperiment()
    print('Loading mzML input file.', end=' ', flush=True)
    ms.MzMLFile().load(args.in_, exp)
    print('Done', flush=True)

    total_t = time.time() - start_t
    time_out.write(f'mzml load: {total_t}s\n')
    time_out.close()

    ff = FeatureFinderIonMobility()
    features = ff.run(exp, args.num_bins, args.pp_type, args.peak_radius, args.window_radius, args.pp_mode,
                      args.ff_type, args.dir, args.filter, args.debug)

    ms.FeatureXMLFile().store(args.dir + '/' + args.out, features)
    print('Found', features.size(), 'features')
