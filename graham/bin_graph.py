import csv
import pyopenms as ms
import feature_finder_im as ffim

exp = ms.MSExperiment()
print('Loading mzML input file.', end=' ', flush=True)
ms.MzMLFile().load('mzML/2768-800-860.mzML', exp)
print('Done', flush=True)

results = []  # List of lists (number of bins, number of features)
ff = ffim.FeatureFinderIonMobility()

features = ff.run(exp, 1, 'pphr', 0, 0, 0, 'centroided', 'runs', 'none', False)
results.append([1, features.size()])
for i in range(5, 76, 5):
    features = ff.run(exp, i, 'pphr', 0, 0, 0, 'centroided', 'runs', 'none', False)
    results.append([i, features.size()])

with open('runs/bin_counts_pphr.csv', 'w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(['num bins', 'num_features'])
    writer.writerows(results)

results = []

features = ff.run(exp, 1, 'custom', 1, 0.015, 'int', 'centroided', 'runs', 'none', False)
results.append([1, features.size()])
for i in range(5, 76, 5):
    features = ff.run(exp, i, 'custom', 1, 0.015, 'int', 'centroided', 'runs', 'none', False)
    results.append([i, features.size()])

with open('runs/bin_counts_custom.csv', 'w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(['num bins', 'num features'])
    writer.writerows(results)
