from im_binning import *

import multiprocessing as mp
import platform

run_num, num_bins = -1, -1
fm1, fm2 = [], []

def load_features(base_fp, run_name):
    """Loads both of the OpenMS FeatureMap objects created by im_binning into the two
    global lists.

    Args:
        base_fp (str): The base filepath for the selected experiment. Subject to change.
        run_name (str): The name of the experiment (--out flag). Subject to change.
    """
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

# Cannot use multithreading effectively due to the GIL, so shmem becomes a problem
# Linux-based systems spawn processes through fork(), so copying globals is unnecessary
if platform.system() == 'Windows':
    def init_globals(m1, m2):
        """Initializes the two global lists in all newly spawned child processes.

        Args:
            m1 (list<FeatureMap>): Pass in fm1.
            m2 (list<FeatureMap>): Pass in fm2.
        """
        global fm1, fm2
        fm1, fm2 = m1, m2

def match_work(rt_start, rt_stop, rt_inc, q):
    """Tries to match/map features in both binning runs together. Only tests a specific
    RT range; for a single MP Process.

    Args:
        rt_start (float): The starting RT value to test.
        rt_stop (float): The last (exclusive) value to test.
        rt_inc (float): How much to increase RT by each step.
        q (Queue): A multiprocessing.Manager.Queue object for communication.
    """
    for rt_threshold in np.arange(rt_start, rt_stop, rt_inc):
        for mz_threshold in np.arange(0.005, 0.5, 0.05):
            features = match_features(fm1, fm2, rt_threshold, mz_threshold)
            q.put([rt_threshold, mz_threshold, features.size()])

    # A sentinel to signal end of writing
    q.put(None)

def match_work_pool(fm1, fm2, rt_threshold, mz_threshold, q):
    """Tries to match/map featurse in both binning runs together. Only tests a specific
    set of parameters; this represents a single task in a work pool.

    Args:
        fm1 (list<FeatureMap>): Pass in fm1.
        fm2 (list<FeatureMap>): Pass in fm2.
        rt_threshold (float): The RT value to test.
        mz_threshold (float): The m/z value to test.
        q (Queue): A multiprocessing.Manager.Queue object for communication
    """
    features = match_features(fm1, fm2, rt_threshold, mz_threshold)
    q.put([rt_threshold, mz_threshold, features.size()])

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

    # ========== WORK POOL IMPLEMENTATION ==========

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
                pool.apply_async(match_work_pool, (fm1, fm2, rt_threshold, mz_threshold,
                                                   q,))

        pool.close()
        pool.join()

        # What if the queue fills up and a child blocks while writing?
        while not q.empty():
            output = q.get()
            print(*output)

    # ========== PROCESS IMPLEMENTATION ==========

    elif args.mp_mode == 2:
        man = mp.Manager()
        q = man.Queue()
        # Maybe leave one core to handle output so the queue never fills up
        procs = [mp.Process(target=match_work, args=(0, 3, 0.25, q,)),
                 mp.Process(target=match_work, args=(3, 6, 0.25, q,)),
                 mp.Process(target=match_work, args=(6, 9, 0.25, q,)),
                 mp.Process(target=match_work, args=(9, 12, 0.25, q,))]

        for p in procs:
            p.start()

        #done = 0
        #while True:
        #    output = q.get()

        #    if output is None:
        #        done += 1
        #        if done == args.nprocs:
        #            break
        #    else:
        #        print(*output)

        for p in procs:
            p.join()

        while not q.empty():
            output = q.get()
            print(*output)

    # ========== SERIAL IMPLEMENTATION ==========

    elif args.mp_mode == 0:
        rt_start, mz_start = 1, 0.005
        # Also keep track of the param set that produces the "best" mapping
        min_features, min_rt, min_mz = float('inf'), -1, -1

        for rt_threshold in np.arange(rt_start, 12, 0.5):
            for mz_threshold in np.arange(mz_start, 0.5, 0.05):
                features = match_features(fm1, fm2, rt_threshold, mz_threshold)
                if features.size() < min_features:
                    min_features = features.size()
                    min_rt, min_mz = rt_threshold, mz_threshold

                print(rt_threshold, mz_threshold, features.size())

        print(min_rt, min_mz, min_features)

        features = match_features(fm1, fm2, min_rt, min_mz)
        ms.FeatureXMLFile().store(args.output + '/common.featureXML', features)
