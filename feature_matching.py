from im_binning import *

run_num, num_bins = -1, -1
fm1, fm2 = [], []

def load_features(base_fp, run_name):
    global fm1, fm2

    for i in range(num_bins):
        fm = ms.FeatureMap()
        ms.FeatureXMLFile().load(base_fp + str(run_num) + '/' + run_name + str(run_num) +
                                 '-pass1-bin' + str(i) + '.featureXML', fm)
        fm1.append(fm)

    for i in range(num_bins + 1):
        fm = ms.FeatureMap()
        ms.FeatureXMLFile().load(base_fp + str(run_num) + '/' + run_name + str(run_num) +
                                 '-pass2-bin' + str(i) + '.featureXML', fm)
        fm2.append(fm)

if __name__ == '__main__':
    parser = argparse.ArgumentParser('Feature matching threshold finder.')
    parser.add_argument('--base_fp', action='store', required=True, type=str)
    parser.add_argument('--run_name', action='store', required=True, type=str)
    parser.add_argument('--run_num', action='store', required=True, type=str)
    parser.add_argument('--num_bins', action='store', required=True, type=int)
    parser.add_argument('--output', action='store', required=True, type=str)

    args = parser.parse_args()
    run_num, num_bins = args.run_num, args.num_bins

    load_features(args.base_fp, args.run_name)

    rt_start, mz_start = 1, 0.005
    min_features, min_rt, min_mz = float('inf'), -1, -1

    f = open(args.output + '/fm_output.txt', 'w+')

    for rt_threshold in np.arange(rt_start, 12, 0.5):
        for mz_threshold in np.arange(mz_start, 0.5, 0.05):
            features = match_features(fm1, fm2, rt_threshold, mz_threshold)
            if features.size() < min_features:
                min_features = features.size()
                min_rt, min_mz = rt_threshold, mz_threshold

            f.write(str(rt_threshold) + ' ' + str(mz_threshold) + '\n')
            f.write(str(features.size()) + '\n')

    f.write('\n' + str(min_rt) + ' ' + str(min_mz) + '\n')
    f.write(str(min_features) + '\n')

    f.close()

    features = match_features(fm1, fm2, min_rt, min_mz)
    ms.FeatureXMLFile().store(args.output + '/common.featureXML', features)
