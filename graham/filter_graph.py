import csv
import pyopenms as ms
import feature_finder_im as ffim

exp = ms.MSExperiment()
print('Loading mzML input file.', end=' ', flush=True)
ms.MzMLFile().load('mzML/2768-800-860.mzML', exp)
print('Done', flush=True)

results = []  # List of lists (ppm, number of features)
ff = ffim.FeatureFinderIonMobility()

for i in range(2, 30, 2):
    features = ff.run(exp, 10, 'pphr', 0, 0, 0, 'centroided', 'runs', 'gauss', False, i)
    results.append([i, features.size()])
    
with open('runs/gauss_ppm_pphr.csv', 'w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(['ppm', 'num features'])
    writer.writerows(results)

results = []

for i in range(3, 15, 1):
    features = ff.run(exp, 10, 'pphr', 0, 0, 0, 'centroided', 'runs', 'sgolay', False, 0, i)
    results.append([i, features.size()])
    
with open('runs/sgolay_pts_pphr.csv', 'w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(['num points', 'num features'])
    writer.writerows(results)
