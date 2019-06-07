from baseline import *

import ransac
from plane_fitting import estimate, is_inlier

import random
import math
from sklearn.cluster import DBSCAN

# Data preprocessing
def get_points(spec):
    ion_mobility_to_peaks = zip(spec.getFloatDataArrays()[0], *spec.get_peaks())
    return [[im, mz, intensity] for im, mz, intensity in ion_mobility_to_peaks]

# Data postprocessing
def rm_outliers(points, coeffs, threshold):
    for i in reversed(range(len(points))):
        if not is_inlier(coeffs, points[i], threshold):
            del points[i]
    return points

def extract_outliers(points, coeffs, threshold):
    outl = []
    for i in reversed(range(len(points))):
        if not is_inlier(coeffs, points[i], threshold):
            outl.append(points[i])
    return outl

# Custom plane fitting (needs testing)
def cus_ransac(points, epsilon, num_iters):
    a, b, c, d = 0, 0, 0, 0
    max_inliers = []
    max_inlier_count = 0

    for i in range(num_iters):
        sample = random.choices(points, k=3)
        p1, p2, p3 = sample[0], sample[1], sample[2]
        inliers = []
        inlier_count = 0

        v1 = (p2[0] - p1[0], p2[1] - p1[1], p2[2] - p1[2])
        v2 = (p3[0] - p1[0], p3[1] - p1[1], p3[2] - p1[2])
        # Normal vector
        v3 = (v1[1] * v2[2] - v1[2] * v2[1],
              v1[2] * v2[0] - v1[0] * v2[2],
              v1[0] * v2[1] - v1[1] * v2[0])
        # Compute D
        cd = v3[0] * p1[0] + v3[1] * p1[1] + v3[2] * p1[2]
        # Coefficients A, B, C, D
        m = (v3[0], v3[1], v3[2], cd)

        for j in range(len(points)):
            s = points[j]
            dist = (m[0] * s[0] + m[1] * s[1] + m[2] * s[2] + m[3]) / (
                   math.sqrt(m[0] * m[0] + m[1] * m[1] + m[2] * m[2]))
            if dist < epsilon:
                inlier_count += 1
                inliers.append(s)

        if (inlier_count > max_inlier_count):
            max_inlier_count = inlier_count
            max_inliers = inliers
            a, b, c, d = m

    return (a, b, c, d), max_inliers

def fit_plane(coords, spec):
    lxyzs = [list(x) for x in coords]
    cxyzs = list(zip(*coords))
    points = [list(cxyzs[0]), list(cxyzs[1]), list(cxyzs[2])]
    
    # RANSAC not working; just use points 10% greater than avg
    avg = 0
    for i in range(len(lxyzs)):
        avg += lxyzs[i][2]
    avg /= len(lxyzs)

    pf = []
    for i in range(len(lxyzs)):
        if lxyzs[i][2] > avg * 1.1:
            pf.append(lxyzs[i])

    print(len(pf))
    lpf = list(zip(*pf))
    llpf = [list(x) for x in lpf]

    #fig = plt.figure()
    #ax = Axes3D(fig)
    #ax.scatter3D(llpf[0], llpf[1], llpf[2])

    # DBSCAN clustering
    db = DBSCAN(eps=1, min_samples=5).fit(lxyzs)
    core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
    core_samples_mask[db.core_sample_indices_] = True
    labels = db.labels_

    n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
    n_noise_ = list(labels).count(-1)

    print("Estimated number of clusters: %d" % n_clusters_)
    print("Estimated number of noise points: %d" % n_noise_)

    #n = len(coords)
    #max_iterations = 100
    #goal_inliers = n * 0.9

    #m, b = ransac.run_ransac(lxyzs, estimate, lambda x, y: is_inlier(x, y, 0.01),
    #                         3, goal_inliers, max_iterations)
    #a, b, c, d = m

    #m, ni = cus_ransac(lxyzs)
    #a, b, c, d = m

    #fig = plt.figure()
    #ax = Axes3D(fig)
    #ax.scatter3D(points[0], points[1], points[2])

    #def plot_plane(a, b, c, d):
    #    xx, yy = np.mgrid[:3, :1750]
    #    return xx, yy, (-d - a * xx - b * yy) / c

    #xx, yy, zz = plot_plane(a, b, c, d)
    #ax.plot_surface(xx, yy, zz, color=(0, 1, 0, 0.5))

    #xyzs_in = rm_outliers(lxyzs, m, 0.01)
    #coords = [tuple(x) for x in xyzs_in]
    #cxyzs = list(zip(*coords))
    #points = [list(cxyzs[0]), list(cxyzs[1]), list(cxyzs[2])]
    #ax.scatter3D(points[0], points[1], points[2])

    #m, best_inliers = ransac.run_ransac(xyzs_in, estimate, lambda x, y: is_inlier(x, y, 0.01),
    #                                    3, goal_inliers, max_iterations)
    #a, b, c, d = m
    #xyzs_in = rm_outliers(lxyzs, m, 0.01)
    #xx, yy, zz = plot_plane(a, b, c, d)
    #ax.plot_surface(xx, yy, zz, color=(0, 1, 0, 0.5))

    #plt.show()

    # Group by cluster
    f_groups = []
    for i in range(len(coords)):
        f_groups.append([])

    for i in range(len(coords)):
        if labels[i] != -1:
            f_groups[labels[i]].append(coords[i])

    features = ms.FeatureMap()

    # Attempt 1: only include apex points
    for i in range(len(coords)):
        if len(f_groups[i]) == 0:
            continue

        apex = f_groups[i][0]
        for j in range(len(f_groups[i])):
            if f_groups[i][j][2] > apex[2]:
                apex = f_groups[i][j]

        f = ms.Feature()
        f.setMZ(apex[1])
        f.setCharge(1)
        f.setRT(spec.getRT())
        f.setIntensity(apex[2])
        f.setOverallQuality(10)

        features.push_back(f)

    features.setUniqueIds()
    return features

# Begin execution
def init(args):
    exp = ms.MSExperiment()
    ms.MzMLFile().load(args.infile + '.mzML', exp)

    counter_to_og_rt_ms = {}
    start_idx = 0
    spectra = exp.getSpectra()

    for i in range(start_idx, start_idx + args.num_frames):
        spec = spectra[i]
        points = get_points(spec)

        # Quick graph of data
        #fig = plt.figure()
        #ax = fig.add_subplot(111, projection='3d')
        #cxyz = list(zip(*coords))
        #x, y, z = list(cxyz[0]), list(cxyz[1]), list(cxyz[2])
        #ax.scatter(x, y, z, c='r', marker='o')
        #ax.set_xlabel('IM')
        #ax.set_ylabel('MZ')
        #ax.set_zlabel('Intensity')
        #plt.show()

        new_features = fit_plane(points, spec)
        ms.FeatureXMLFile().store(args.outdir + '/' + str(i) + '_' + args.outfile +
                                  '.featureXML', new_features)

        counter_to_og_rt_ms[i] = [spec.getRT(), spec.getMSLevel()]

    with open(args.outdir + '/counter_to_og_rt_ms.pkl', 'wb') as handle:
        pickle.dump(counter_to_og_rt_ms, handle, protocol=pickle.HIGHEST_PROTOCOL)

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

    init(args)
    driver(args)
