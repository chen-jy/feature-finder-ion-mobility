# FragToPre
Approaches to linking fragments to their precursors in mass spec data using four-dimensional feature finding (with retention time, mass to charge, intensity, and ion mobility).

Tools:

**baseline**: splits data into frames by RT, then for each frame, swaps RT with ion mobility. FeatureFinderCentroided is then run and the results are linked together across frames.
```
export OMP_NUM_THREADS=1
python baseline.py --infile sample_3001_to_3098 --outfile rt --outdir run --mz_epsilon 2.5 --im_epsilon 0.025 --num_frames 833 --window_size 17 --rt_length 2
```

**im_binning**: splits data into frames by RT, then for each frame, bins data points in the IM dimension. All of the bins are then linked together, and FeatureFinderCentroided (with optional hi-res peak picking) is run.
```
python im_binning.py --infline sample --outfile rt --outdir run
```

**cmp_features**: compares features in two featureXML files. Thresholds in each dimension are used to determine similarity, and the number of intersecting features is found. Used to compare features found by im_binning with those found by OpenMS (via pyOpenMS).
```
python cmp_features --found my_features.featureXML --openms oms_features.featureXML --out results

# If an OpenMS-run featureXML cannot be provided, the tool can run FeatureFinderCentroided
# on the MS1 data of an mzML file
python cmp_features --found my_features.featureXML --source ms_data.mzML --out results
```
