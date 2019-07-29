from im_binning import *

import multiprocessing as mp

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

def match_work(rt_start, rt_stop, rt_inc, q):
    for rt_threshold in np.arange(rt_start, rt_stop, rt_inc):
        for mz_threshold in np.arange(0.005, 0.5, 0.05):
            features = match_features(fm1, fm2, rt_threshold, mz_threshold)
            q.put([rt_threshold, mz_threshold, features.size()])

    q.put('DONE')

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

    q = mp.Queue()
    procs = [mp.Process(target=match_work, args=(0, 4, 0.5, q,)),
             mp.Process(target=match_work, args=(4, 8, 0.5, q,)),
             mp.Process(target=match_work, args=(8, 12, 0.5, q,))]

    for i in range(1):
        procs[i].start()

    procs_done = 0
    while True:
        msg = q.get()
        if (msg == 'DONE'):
            procs_done += 1
            if (procs_done == 1):
                break
            continue

        print(msg[0], msg[1], msg[2])

    for i in range(len(procs)):
        procs[i].join()

    #rt_start, mz_start = 1, 0.005
    #min_features, min_rt, min_mz = float('inf'), -1, -1

    #f = open(args.output + '/fm_output.txt', 'w+')

    #for rt_threshold in np.arange(rt_start, 12, 0.5):
    #    for mz_threshold in np.arange(mz_start, 0.5, 0.05):
    #        features = match_features(fm1, fm2, rt_threshold, mz_threshold)
    #        if features.size() < min_features:
    #            min_features = features.size()
    #            min_rt, min_mz = rt_threshold, mz_threshold

    #        f.write(str(rt_threshold) + ' ' + str(mz_threshold) + '\n')
    #        f.write(str(features.size()) + '\n')

    #f.write('\n' + str(min_rt) + ' ' + str(min_mz) + '\n')
    #f.write(str(min_features) + '\n')

    #f.close()

    #features = match_features(fm1, fm2, min_rt, min_mz)
    #ms.FeatureXMLFile().store(args.output + '/common.featureXML', features)
