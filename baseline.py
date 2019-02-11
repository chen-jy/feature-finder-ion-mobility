import argparse
import pickle
from collections import defaultdict

import matplotlib.pyplot as plt
import numpy as np
import pyopenms as ms
from mpl_toolkits.mplot3d import Axes3D
from scipy.signal import argrelextrema

ISOLATION_WINDOWS = \
    {0: 0, 1: 437.5, 2: 487.5, 3: 537.5, 4: 587.5, 5: 637.5, 6: 687.5, \
     7: 737.5, 8: 787.5, 9: 837.5, 10: 887.5, 11: 937.5, 12: 987.5, \
     13: 1037.5, 14: 1087.5, 15: 1137.5, 16: 1187.5}

def four_d_spectrum_to_experiment(spec):
    """Function that converts a 4D spectrum object which contains retention
    time, ion mobility, mass to charge, and intensity data, into a new
    experiment where the ion mobility becomes the retention time of
    each of its new spectra and vice versa.

    Args:
        spec (MSSpectrum): An OpenMS MSSpectrum object.

    Returns:
        MSExperiment: A new MSExperiment where each MSSpectrum object is a 3D
        spectrum from the original where the retention time and ion mobility
        has been swapped.
    """
    ion_mobility_to_peaks = zip(spec.getFloatDataArrays()[0], *spec.get_peaks())

    new_exp = ms.MSExperiment()
    new_mz, new_int = [], []
    curr_im = None

    for im, mz, intensity in sorted(
        ion_mobility_to_peaks, key = lambda x: (x[0], x[1], x[2])):
        if im != curr_im:
            if curr_im:
                new_spec = ms.MSSpectrum()
                new_spec.setRT(curr_im)
                new_spec.set_peaks((new_mz, new_int))

                rt_fda = ms.FloatDataArray()
                for i in new_mz:
                    rt_fda.push_back(spec.getRT())

                new_spec.setFloatDataArrays([rt_fda])
                new_exp.addSpectrum(new_spec)
                new_mz, new_int = [], []
            
            curr_im = im
            
        new_mz.append(mz)
        new_int.append(intensity)

    return new_exp

def run_feature_finder_centroided_on_experiment(input_map):
    """Function that runs FeatureFinderCentroided on the given input map.

    Args:
        input_map (MSExperiment): An OpenMS MSExperiment object.

    Returns:
        FeatureMap: A FeatureMap containing the found features from the given
        experiment.
    """
    # Load data
    input_map.updateRanges()

    ff = ms.FeatureFinder()
    ff.setLogType(ms.LogType.CMD)

    # Run the feature finder
    name = 'centroided'
    features = ms.FeatureMap()
    seeds = ms.FeatureMap()
    params = ms.FeatureFinder().getParameters(name)
    params.__setitem__(b'mass_trace:min_spectra', 5)
    params.__setitem__(b'mass_trace:max_missing', 2)
    params.__setitem__(b'seed:min_score', 0.5)
    params.__setitem__(b'feature:min_score', 0.5)
    ff.run(name, input_map, features, params, seeds)

    features.setUniqueIds()

    return features

def find_local_maxima_indices(feature_arrays):
    """Function that finds local intensity maxima.

    Args:
        feature_arrays: Array of arrays of points.

    Returns:
        Array of tuples, where the first element is
        the index of the maxima, and the second element
        is the point index.
    """
    max_ints, max_int_pt_idxs = [], []

    for i in range(len(feature_arrays)):
        if len(feature_arrays[i]) > 0:
            pt_idx, max_val = max(
                enumerate(feature_arrays[i]), key=lambda p:p[1][2])

            max_ints.append(max_val[2])
            max_int_pt_idxs.append(pt_idx)
        else:
            max_ints.append(0)
            max_int_pt_idxs.append(None)

    max_ints, max_int_pt_idxs = np.array(max_ints), np.array(max_int_pt_idxs)

    # determine the indices of the local maxima
    max_int_idxs = argrelextrema(max_ints, np.greater)

    # get the actual values using these indices
    local_max_idxs = \
        zip(*max_int_idxs, max_int_pt_idxs[max_int_idxs])

    return local_max_idxs

def link_to_peak(
    rt_idx, pt_idx, \
    rt_idx_to_points, \
    rt_idx_to_num_points_left, \
    possible_species, \
    species_counter, \
    rt_idx_to_rt, \
    mz_epsilon, \
    im_epsilon):
    """Function links points adjacent to a local maximum.
       Tries to link everything in a peak shape (e.g. decreasing on both sides)
    """
    if rt_idx_to_num_points_left[rt_idx] > 0:
        peak = rt_idx_to_points[rt_idx].pop(pt_idx)
        possible_species[species_counter] = [[rt_idx_to_rt[rt_idx]], [peak]]

        rt_idx_to_num_points_left[rt_idx]-= 1

        walk_idx = rt_idx - 1
        skipped = 0
        while walk_idx >= 0:
            if len(rt_idx_to_points[walk_idx]) == 0:
                break

            prev_val = peak[2]
            point_found = False

            for i in range(len(rt_idx_to_points[walk_idx])):
                point = rt_idx_to_points[walk_idx][i]

                if peak[0] - mz_epsilon <= point[0] <= peak[0] + mz_epsilon \
                and peak[1] - im_epsilon <= point[1] <= peak[1] + im_epsilon \
                and point[2] <= peak[2] and point[3] == peak[3] \
                and point[2] < prev_val:
                    point = rt_idx_to_points[walk_idx].pop(i)
                    possible_species[species_counter][0].append(
                        rt_idx_to_rt[walk_idx])
                    possible_species[species_counter][1].append(point)
                    rt_idx_to_num_points_left[walk_idx]-= 1
                    prev_val = point[2]
                    point_found = True
                    walk_idx-= 1
                    break
            
            if not point_found:
                skipped+= 1
                if skipped > 2:
                    break
                else:
                    walk_idx-= 1


        walk_idx = rt_idx + 1
        skipped = 0
        while walk_idx < len(rt_idx_to_points):
            if len(rt_idx_to_points[walk_idx]) == 0:
                break

            prev_val = peak[2]
            point_found = False

            for i in range(len(rt_idx_to_points[walk_idx])):
                point = rt_idx_to_points[walk_idx][i]

                if peak[0] - mz_epsilon <= point[0] <= peak[0] + mz_epsilon \
                and peak[1] - im_epsilon <= point[1] <= peak[1] + im_epsilon \
                and point[2] <= peak[2] and point[3] == peak[3] \
                and point[2] < prev_val:
                    point = rt_idx_to_points[walk_idx].pop(i)
                    possible_species[species_counter][0].append(
                        rt_idx_to_rt[walk_idx])
                    possible_species[species_counter][1].append(point)
                    rt_idx_to_num_points_left[walk_idx]-= 1
                    prev_val = point[2]
                    point_found = True
                    walk_idx+= 1
                    break
            
            if not point_found:
                skipped+= 1
                if skipped > 2:
                    break
                else:
                    walk_idx+= 1

def link_between_frames(feature_maps, rt_idx_to_rt, mz_epsilon, im_epsilon):
    """Function that extracts possible species of targets from feature maps.
    """
    rt_idx_to_points = []
    rt_idx_to_num_points_left, possible_species = {}, {}
    num_unclassified, species_counter = None, 0

    for i in range(len(feature_maps)):
        rt = []

        for feature in feature_maps[i]:
            point = []

            point.append(feature.getMZ())
            point.append(feature.getRT())
            point.append(feature.getIntensity())
            point.append(feature.getCharge())

            rt.append(point)

        rt_idx_to_points.append(rt)
        rt_idx_to_num_points_left[i] = len(rt)

    num_unclassified = len([key for key in rt_idx_to_num_points_left \
        if rt_idx_to_num_points_left[key] > 0])
    while num_unclassified > 0:
        local_maxima = list(find_local_maxima_indices(rt_idx_to_points))

        if not local_maxima:
            break

        for rt_idx, pt_idx in local_maxima:
            link_to_peak(
                rt_idx, \
                pt_idx, \
                rt_idx_to_points, \
                rt_idx_to_num_points_left, \
                possible_species, \
                species_counter, \
                rt_idx_to_rt, \
                mz_epsilon, \
                im_epsilon)

            species_counter+= 1

    for key in [key for key in rt_idx_to_num_points_left \
        if rt_idx_to_num_points_left[key] > 0]:
        points = rt_idx_to_points[key]

        for point in points:
            possible_species[species_counter] = [[rt_idx_to_rt[key]], [point]]
            species_counter+= 1

    return possible_species

def split_precursors_and_fragments(
    possible_species, window_size, rt_length, counter_to_og_rt_ms):
    """Function that splits the list of possible species into precursors and
    fragments.
    """
    precursors, fragments = {}, {}
    precursor_counter, fragment_counter = 0, 0

    for species in possible_species:
        rts = possible_species[species][0]

        if len(rts) > rt_length:
            for rt in rts:
                if counter_to_og_rt_ms[rt][1] == 1:
                    precursors[precursor_counter] = possible_species[species]
                    precursor_counter+= 1
                    break
                else:
                    fragments[fragment_counter] = possible_species[species]
                    fragment_counter+= 1
                    break

    return precursors, fragments

def link_frag_to_prec(dir, fragments, precursors, window_size, im_epsilon, threshold):
    precursor_to_fragments = defaultdict(list)

    for precursor in precursors:
        precursor_rt_range, precursor_points = precursors[precursor]

        for fragment in fragments:
            fragment_rt_range, fragment_points = fragments[fragment]
            isolation_mz = ISOLATION_WINDOWS[fragment_rt_range[0] % window_size]
            lower, upper = 25, 25

            if min(precursor_rt_range) < fragment_rt_range[0] < max(precursor_rt_range) \
            and np.abs(
                np.mean(
                    [x[1] for x in precursor_points]) - np.mean(
                    [x[1] for x in fragment_points])) <= im_epsilon \
            and isolation_mz - lower <= np.mean(
                [x[0] for x in precursor_points]) <= isolation_mz + upper:
                precursor_to_fragments[precursor].append(fragment)

    with open(dir + '/results.txt', 'w') as outfile:
        for precursor in precursor_to_fragments:
            outfile.write(str(precursor) + ': ' + str(precursor_to_fragments[precursor]) + '\r\n')

    with open(dir + '/results.pkl', 'wb') as handle:
        pickle.dump(precursor_to_fragments, handle, protocol=pickle.HIGHEST_PROTOCOL)

    return precursor_to_fragments

def plot_3d_intensity_map(feature_maps, rt_idx_to_rt, mode='3d'):
    rt, mz, im, intensity = [], [], [], []

    for i in range(len(feature_maps)):
        for feature in feature_maps[i]:
            rt.append(rt_idx_to_rt[i])
            mz.append(feature.getMZ())
            im.append(feature.getRT())
            intensity.append(feature.getIntensity())

    fig = plt.figure()

    if mode == '3d':
        ax = fig.add_subplot(111, projection='3d')
        line = ax.scatter(mz, rt , im , c=intensity, s=25, marker='.', edgecolors='none', depthshade=0)
        ax.set_xlabel('mz')
        ax.set_ylabel('rt')
        ax.set_zlabel('im')
    else:
        ax = fig.add_subplot(111)
        line = ax.scatter(mz, rt, c=intensity, s=25, marker='.', edgecolors='none')
        ax.set_xlabel('mz')
        ax.set_ylabel('rt')

    cb = plt.colorbar(line)

    plt.show()

def parse_openms(openms_filename):
    prec_id_to_mz, prec_id_to_frag_anno = {}, {}
    with open(openms_filename, 'r') as infile:
        next(infile)
        for line in infile:
            line = line.rstrip('\n')
            line = line.split(',')
            prec_id_to_mz[line[1]] = float(line[9])
            prec_id_to_frag_anno[line[1]] = line[60]
                    
    return prec_id_to_mz, prec_id_to_frag_anno

def process_openms_frag_anno(anno):
    anno = [float(x.split('_')[2]) for x in anno.split(';')]

    return anno

def compare_baseline_to_openms(
    args, precursors, fragments, results, openms_filename):
    num_openms = 0
    shared_pairs = []

    openms_id_to_mz, openms_id_to_frag_anno = parse_openms(openms_filename)

    for openms_id in openms_id_to_mz:
        openms_precursor_mz = openms_id_to_mz[openms_id]
        match_found = False

        for precursor in precursors:
            _, precursor_points = precursors[precursor]
            precursor_mz = np.mean([x[0] for x in precursor_points])

            if np.abs(precursor_mz - openms_precursor_mz) < args.mz_epsilon:
                shared_pairs.append((openms_id, precursor))
                match_found = True

        if match_found:
            num_openms+= 1

    print('Number of precursors detected in both OpenMS and Baseline: ' + \
        str(num_openms) + ' out of ' + str(len(openms_id_to_mz)) + '\n\r')
    print(len(precursors))

    num_openms = 0
    openms_frag_mzs = []

    for openms_id in openms_id_to_mz:
        for openms_frag_mz in process_openms_frag_anno(
            openms_id_to_frag_anno[openms_id]):
            openms_frag_mzs.append(openms_frag_mz)

    openms_frag_mzs = list(set(openms_frag_mzs))
    openms_frag_mzs.sort()

    baseline_frag_mzs = []

    for fragment in fragments:
        _, fragment_points = fragments[fragment]
        fragment_mz = np.mean([x[0] for x in fragment_points])
        baseline_frag_mzs.append(fragment_mz)

    baseline_frag_mzs = list(set(baseline_frag_mzs))
    baseline_frag_mzs.sort()

    i, j = 0, 0

    while i < len(openms_frag_mzs) and j < len(baseline_frag_mzs):
        openms_frag_mz, baseline_frag_mz = (
            openms_frag_mzs[i], baseline_frag_mzs[j])

        if np.abs(baseline_frag_mz - openms_frag_mz) < args.mz_epsilon:
            num_openms+= 1
            i+= 1
            j+= 1
        elif min(openms_frag_mz, baseline_frag_mz) == openms_frag_mz:
            i+= 1
        else:
            j+= 1

    # for openms_frag_mz in openms_frag_mzs:
    #     if counter % 100 == 0:
    #         print(counter)
    #     counter+= 1
    #     match_found = False

    #     for fragment in fragments:
    #         _, fragment_points = fragments[fragment]
    #         fragment_mz = np.mean([x[0] for x in fragment_points])

    #         if np.abs(fragment_mz - openms_frag_mz) < args.mz_epsilon:
    #             match_found = True

    #     if match_found:
    #         num_openms+= 1

    print('Number of fragments detected in both OpenMS and Baseline: ' + \
        str(num_openms) + ' out of ' + str(len(openms_frag_mzs)) + '\n\r')
    
    baseline_frag_mzs = list(set([float('%.6g' % x) for x in baseline_frag_mzs]))
    print(len(baseline_frag_mzs))

    # shared_pairs_to_counts = {}

    # for pair in shared_pairs:
    #     openms_id, prec_id = pair
    #     openms_frag_mzs = process_openms_frag_anno(
    #         openms_id_to_frag_anno[openms_id])
    #     baseline_frags = results[prec_id]

    #     total = len(openms_frag_mzs)
    #     matching, extra = 0, 0

    #     for openms_frag_mz in openms_frag_mzs:
    #         match_found = False

    #         for baseline_frag in baseline_frags:
    #             _, frag_points = fragments[baseline_frag]
    #             frag_mz = np.mean([x[0] for x in frag_points])

    #             if np.abs(frag_mz - openms_frag_mz) < args.mz_epsilon:
    #                 match_found = True

    #         if match_found:
    #             matching+= 1
    #         else:
    #             extra+= 1

    #     shared_pairs_to_counts[pair] = [total, matching, extra]

    # print(len(set([x[0] for x in shared_pairs if shared_pairs_to_counts[x][1] > 5])))
    # print('\n\r')
    # print(len(set([x[0] for x in shared_pairs if shared_pairs_to_counts[x][1] > 0])))
    # print('\n\r')

    # for pair in shared_pairs:
    #     if shared_pairs_to_counts[pair][1] > 5:
    #         print(str(pair) + ':' + str(shared_pairs_to_counts[pair]) + '\n\r')
    #         openms_id, prec_id = pair
    #         openms_frag_mzs = process_openms_frag_anno(
    #             openms_id_to_frag_anno[openms_id])
    #         print(str(openms_frag_mzs) + '\n\r')
    #         baseline_frags = results[prec_id]
    #         baseline_frag_mzs = []
    #         for baseline_frag in baseline_frags:
    #             _, frag_points = fragments[baseline_frag]
    #             frag_mz = np.mean([x[0] for x in frag_points])
    #             baseline_frag_mzs.append(frag_mz)
    #         print(str(sorted(baseline_frag_mzs)) + '\n\r')

    # return shared_pairs_to_counts
    return

def driver(args):
    # exp = ms.MSExperiment()
    # ms.MzMLFile().load(args.infile + '.mzML', exp)

    # counter_to_og_rt_ms = {}

    # start_idx = 0
    # spectra = exp.getSpectra()
    # for i in range(start_idx, start_idx + args.num_frames):
    #     spec = spectra[i]
    #     new_exp = four_d_spectrum_to_experiment(spec)
    #     ms.MzMLFile().store(args.outdir + '/' + str(i) + '_' + args.outfile + '.mzML', new_exp)

    #     new_features = run_feature_finder_centroided_on_experiment(new_exp)
    #     ms.FeatureXMLFile().store(args.outdir + '/' + str(i) + '_' + args.outfile + '.featureXML', new_features)

    #     counter_to_og_rt_ms[i] = [spec.getRT(), spec.getMSLevel()]

    # with open(args.outdir + '/counter_to_og_rt_ms.pkl', 'wb') as handle:
    #     pickle.dump(counter_to_og_rt_ms, handle, protocol=pickle.HIGHEST_PROTOCOL)

    # #####################################################################################

    # feature_maps = []
    # rt_idx_to_rt, counter_to_og_rt_ms = {}, {}
    # counter = 0

    # with open(args.outdir + '/counter_to_og_rt_ms.pkl', 'rb') as handle:
    #     counter_to_og_rt_ms = pickle.load(handle)

    # with open(args.outdir + '/frames.txt', 'w') as infile:
    #     for i in range(0, args.window_size):
    #         for j in range(i, args.num_frames, args.window_size):
    #             infile.write(str(j) + "\r\n")
    #             features = ms.FeatureMap()
    #             ms.FeatureXMLFile().load(
    #                 args.outdir + '/' + str(j) + '_' + args.outfile + '.featureXML',
    #                 features)

    #             for feature in features:
    #                 infile.write(str(feature.getMZ()) \
    #                 + ',' + str(feature.getRT()) + ',' \
    #                 + str(feature.getIntensity()) + "\r\n")

    #             feature_maps.append(features)
    #             rt_idx_to_rt[counter] = j
    #             counter+= 1

    # possible_species = \
    #     link_between_frames(feature_maps, rt_idx_to_rt, args.mz_epsilon, args.im_epsilon)
    # # plot_3d_intensity_map(feature_maps, rt_idx_to_rt)

    # precursors, fragments = split_precursors_and_fragments(
    #     possible_species, args.window_size, args.rt_length, counter_to_og_rt_ms)

    # with open(args.outdir + '/precursors.txt', 'w') as outfile:
    #     outfile.write(str(len(precursors)) + "\r\n")
    #     for precursor in precursors:
    #         outfile.write(str(precursor) + ": " + str(precursors[precursor]) + "\r\n")

    # with open(args.outdir + '/precursors.pkl', 'wb') as handle:
    #     pickle.dump(precursors, handle, protocol=pickle.HIGHEST_PROTOCOL)

    # with open(args.outdir + '/fragments.txt', 'w') as outfile:
    #     outfile.write(str(len(fragments)) + "\r\n")
    #     for fragment in fragments:
    #         outfile.write(str(fragment) + ": " + str(fragments[fragment]) + "\r\n")

    # with open(args.outdir + '/fragments.pkl', 'wb') as handle:
    #     pickle.dump(fragments, handle, protocol=pickle.HIGHEST_PROTOCOL)

    # ####################################################################################

    # link_frag_to_prec(args.outdir, fragments, precursors, args.window_size, args.im_epsilon, 0)

    with open(args.outdir + '/precursors.pkl', 'rb') as handle:
        precursors = pickle.load(handle)

    with open(args.outdir + '/fragments.pkl', 'rb') as handle:
        fragments = pickle.load(handle)

    with open(args.outdir + '/results.pkl', 'rb') as handle:
        results = pickle.load(handle)

    compare_baseline_to_openms(args, precursors, fragments, results, '../../data/3001_to_3098_openms.csv')


if __name__ == "__main__":
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

    driver(args)
