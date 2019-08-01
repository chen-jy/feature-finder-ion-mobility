from im_binning import *

import multiprocessing as mp
import platform

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

# Linux-based systems spawn processes through fork(), so copying globals is unnecessary
if platform.system() == 'Windows':
    def init_globals(m1, m2):
        global fm1, fm2
        fm1, fm2 = m1, m2

def match_work(rt_start, rt_stop, rt_inc, q):
    for rt_threshold in np.arange(rt_start, rt_stop, rt_inc):
        for mz_threshold in np.arange(0.005, 0.5, 0.05):
            features = match_features(fm1, fm2, rt_threshold, mz_threshold)
            q.put([rt_threshold, mz_threshold, features.size()])

    # A sentinel to signal end of writing
    q.put(None)

def match_work_q(fm1, fm2, rt, mz, q):
    features = match_features(fm1, fm2, rt, mz)
    q.put([rt, mz, features.size()])

if __name__ == '__main__':
    parser = argparse.ArgumentParser('Feature matching threshold finder.')
    # TODO: reduce the number of args
    parser.add_argument('--base_fp', action='store', required=True, type=str)
    parser.add_argument('--run_name', action='store', required=True, type=str)
    parser.add_argument('--run_num', action='store', required=True, type=str)
    parser.add_argument('--num_bins', action='store', required=True, type=int)
    parser.add_argument('--output', action='store', required=True, type=str)
    parser.add_argument('--nprocs', action='store', required=True, type=int)
    parser.add_argument('--mp_mode', action='store', required=True, type=int)

    args = parser.parse_args()
    run_num, num_bins = args.run_num, args.num_bins
    load_features(args.base_fp, args.run_name)

    # ========== POOL IMPLEMENTATION ==========

    if args.mp_mode == 1:
        man = mp.Manager()
        q = man.Queue()
        rt_start, mz_start = 1, 0.005

        if platform.system() == 'Windows':
            pool = mp.Pool(args.nprocs, init_globals(fm1, fm2))
        else:
            pool = mp.Pool(args.nprocs)

        for rt_threshold in np.arange(rt_start, 12, 0.25):
            for mz_threshold in np.arange(mz_start, 0.5, 0.05):
                pool.apply_async(match_work_q, (fm1, fm2, rt_threshold, mz_threshold,
                                                q,))

        pool.close()
        pool.join()

        # What if the queue fills up and a child blocks while writing?
        while not q.empty():
            output = q.get()
            print(*output)

    # ========== PROCESS IMPLEMENTATION ==========

    elif args.mp_mode == 2:
        q = mp.Queue()
        # Maybe leave one core to handle output so the queue never fills up
        procs = [mp.Process(target=match_work, args=(0, 3, 0.25, q,)),
                 mp.Process(target=match_work, args=(3, 6, 0.25, q,)),
                 mp.Process(target=match_work, args=(6, 9, 0.25, q,)),
                 mp.Process(target=match_work, args=(9, 12, 0.25, q,))]

        for p in procs:
            p.start()
        for p in procs:
            p.join()

        done = 0
        while True:
            output = q.get()

            if output is None:
                done += 1
                if done == args.nprocs:
                    break
            else:
                print(*output)

    # ========== SERIAL IMPLEMENTATION ==========

    elif args.mp_mode == 0:
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
