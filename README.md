# FragToPre

**Warning: this README was last updated 2019-07-25 and is probably outdated; it is not updated regularly.**

Approaches to four-dimensional mass spectrometry feature finding (with retention time, mass to charge, intensity, and ion mobility data).

Utilities:

**baseline**: splits raw mzML data into frames by RT, then for each frame, swaps RT with IM. Next, FeatureFinderCentroided is run, and the results are linked together across frames.

Some arguments (obvious arguments will be skipped from hereon in):
- mode: 0 to find features and link fragments to precursors, 1 to only find features, 2 to only link fragments to precursors
- *Required*: infile, outfile, outdir, mode
- *Note*: outdir must already exist
```
export OMP_NUM_THREADS=1
python baseline.py --infile sample --outfile rt --outdir run1 --mode 1 --mz_epsilon 2.5 --im_epsilon 0.025 --num_frames 833 --window_size 17 --rt_length 2
```

**im_binning**: splits raw mzML data into frames by RT, then for each frame, bins data points in the IM dimension. All of the bins are then linked together and FeatureFinderCentroided (with optional hi-res peak picking) is run.

Some arguments:
- num_bins: the number of bins to use (int)
- peak_pick: run PeakPickerHiRes after binning, if this flag is present
- *Required*: infile, outfile, outdir, num_bins
- *Note*: outdir must already exist
```
python im_binning.py --infile sample --outfile im --outdir run2 --num_bins 10 --peak_pick
```

**cmp_features**: compares features in four featureXML files (intended for features found by im_binning, OpenMS's FeatureFinderCentroided, baseline, and the source file's corresponding .tsv file). Thresholds in the RT and m/z dimensions are used to determine "similarity", and the number of intersecting features is found.

Some arguments:
- found: the featureXML file produced by im_binning
- openms: the featureXML file produced by OpenMS
- source: the raw mzML file; only required if FeatureFinderCentroided has not already been run (the OpenMS featureXML file does not exist)
- baseline: the featureXML file produced by baseline
- truth: the featureXML file (filtered and transformed by csvFilter and translate_features) of the raw data's corresponding .tsv file
- brute_force: brute-forces various threshold values, if this flag is present
- *Required*: found, openms OR source, baseline, truth, outdir
- *Note*: outdir must already exist
```
python cmp_features.py --found run1/im --source sample --baseline run2/rt --truth filtered --outdir run3 --brute_force
```

**pphr_param_finder**: finds suitable parameters for running PeakPickerHiRes on IM-binned data by brute force.

Some arguments:
- target: the threshold number of features to be detected in order to be considered a suitable parameter set
- *Required*: input, output
- *Note*: output (dir) must already exist
```
python pphr_param_finder.py --input run1/im-pass1-bin1 --output run4 --target 250
```

**translate_features**: converts features in .csv format to .featureXML format.

Some arguments:
- *Required*: input, output
```
python translate_features.py --input sample_csv --output filtered
```

**csvFilter**: filters data from a .csv file (which has been converted from a .tsv file). Either filters by initial_peak_quality (q) or extracts RT, m/z, and intensity data; i.e. filter performs a "SELECT * FROM sample_csv WHERE q > *threshold*" and extract performs a "SELECT RT, m/z, intensity FROM sample_csv".

Some arguments:
- mode: 0 to filter by q, 1 to extract triplet data
- *Required*: mode
- *Note*: the threshold value can be changed within the source code
```
g++ csvFilter.cpp -o csvFilter
./csvFilter <mode>
```

**Legacy tools**: cluster_finder, plane_fitting, and ransac. This is an old approach and is no longer supported or documented.
