import argparse

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyopenms as ms
import seaborn as sns
from mpl_toolkits.mplot3d import Axes3D
from scipy.signal import argrelextrema


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
    params.__setitem__(b'intensity:bins', 1)
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

def split_precursors_and_fragments(possible_species, window_size, rt_length):
    """Function that splits the list of possible species into precursors and
    fragments.
    """
    precursors, fragments = {}, {}
    precursor_counter, fragment_counter = 0, 0

    for species in possible_species:
        rts = possible_species[species][0]

        if len(rts) > rt_length:
            for rt in rts:
                if rt % window_size == 0:
                    precursors[precursor_counter] = possible_species[species]
                    precursor_counter+= 1
                    break
                else:
                    fragments[fragment_counter] = possible_species[species]
                    fragment_counter+= 1
                    break

    return precursors, fragments

def link_frag_to_prec(fragments, precursors):
    pass

def plot_3d_intensity_map(feature_maps, rt_idx_to_rt):
    rt, mz, im, intensity = [], [], [], []

    for i in range(len(feature_maps)):
        for feature in feature_maps[i]:
            rt.append(rt_idx_to_rt[i])
            mz.append(feature.getMZ())
            im.append(feature.getRT())
            intensity.append(feature.getIntensity())

    fig = plt.figure()
    ax = fig.add_subplot(111,projection='3d')
    line = ax.scatter(mz, rt , im , c=intensity, s=25, marker='.', edgecolors='none', depthshade=0)
    ax.set_xlabel('mz')
    ax.set_ylabel('rt')
    ax.set_zlabel('im')
    cb = plt.colorbar(line)

    plt.show()
        
def driver(args):
    # exp = ms.MSExperiment()
    # ms.MzMLFile().load(args.infile + '.mzML', exp)

    # counter_to_og_rt_and_mslevel = {}

    # counter = 0
    # for spec in exp:
    #     new_exp = four_d_spectrum_to_experiment(spec)
    #     ms.MzMLFile().store(args.outdir + '/' + str(counter) + '_' + args.outfile + '.mzML', new_exp)

    #     new_features = run_feature_finder_centroided_on_experiment(new_exp)
    #     ms.FeatureXMLFile().store(args.outdir + '/' + str(counter) + '_' + args.outfile + '.featureXML', new_features)

    #     counter_to_og_rt_and_mslevel[counter] = [spec.getRT(), spec.getMSLevel()]
        
    #     counter+= 1
        
    # with open(args.outdir + '/counter_to_og_rt_and_mslevel.txt', 'w') as infile:
    #     infile.write(str(counter_to_og_rt_and_mslevel))

    feature_maps = []
    rt_idx_to_rt = {}
    counter = 0

    with open('frames.txt', 'w') as infile:
        for i in range(0, 17):
            for j in range(i, args.num_frames, 17):
                infile.write(str(j) + "\r\n")
                features = ms.FeatureMap()
                ms.FeatureXMLFile().load(str(j) + '_' + args.outfile + '.FeatureXML', features)

                for feature in features:
                    infile.write(str(feature.getMZ()) + ',' + str(feature.getRT()) + ',' + str(feature.getIntensity()) + "\r\n")

                feature_maps.append(features)
                rt_idx_to_rt[counter] = j
                counter+= 1

    possible_species = link_between_frames(feature_maps, rt_idx_to_rt, args.mz_epsilon, args.im_epsilon)
    plot_3d_intensity_map(feature_maps, rt_idx_to_rt)

    precursors, fragments = split_precursors_and_fragments(
        possible_species, args.window_size, args.rt_length)

    print(str(len(precursors)) + "\r\n")
    for precursor in precursors:
        print(str(precursor) + ": " + str(precursors[precursor]) + "\r\n")

    print(str(len(fragments)) + "\r\n")
    for fragment in fragments:
        print(str(fragment) + ": " + str(fragments[fragment]) + "\r\n")


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