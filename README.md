# FragToPre

**Warning: this README was last updated 2019-08-21 and is probably outdated; it is not updated regularly.**

Approaches to four-dimensional mass spectrometry feature finding (with retention time, mass to charge, intensity, and ion mobility data).

Descriptions:

**baseline**: splits raw mzML data into frames by RT, then for each frame, swaps RT data with IM data. Next, FeatureFinderCentroided is run, and the results are linked together across frames.

- infile*: the name of the source mzML file without the file extension (no ".mzML")
- outfile*: the name to group the output by
- outdir*: the directory to output to (it must already exist)
- mode*: `0` to find features and link fragments to precursors, `1` to only find features, or `2` to only link fragments to precursors

```
export OMP_NUM_THREADS=1
python baseline.py --infile sample --outfile rt --outdir run1 --mode 1
```

**im_binning**: splits raw mzML data into frames by RT, then for each frame, bins data points in the IM dimension. All of the bins are then linked together and FeatureFinderCentroided (with optional peak picking) is run.

- infile*: the name of the source mzML file without the file extension
- outfile*: the name to group the output by
- outdir*: the directory to output to (it must already exist)
- num_bins*: the number of bins to use
- peak_pick: `0` to skip peak picking, `1` to use PeakPickerHiRes after binning, or `2` to use the custom peak picker (peak_picker.py) after binning
- match_only: if an argument is present for this flag, only the feature matcher will be run (it assumes that binning has already completed and all required files are present)

```
python im_binning.py --infile sample --outfile im --outdir run2 --num_bins 10 --peak_pick 0
```

**cmp_features**: compares features in four featureXML files (intended for features found by im_binning, FeatureFinderCentroided, baseline, and the source file's corresponding tsv file). Thresholds in the RT and m/z dimensions are used to determine "similarity", and the number of intersecting features is found.

- found*: the name of the featureXML file produced by im_binning, without the file extension
- openms: the featureXML file produced by FeatureFinderCentroided
- source: the raw mzML file (either --openms or --source must be used)
- baseline*: the featureXML file produced by baseline
- truth*: the featureXML file (filtered by csvFilter and translated by translate_features) of the raw data's corresponding tsv file
- outdir*: the directory to output to
- brute_force: if an argument is present for this flag, various threshold values will be tested

```
python cmp_features.py --found run1/im --source sample --baseline run2/rt --truth filtered --outdir cmp
```

**pphr_param_finder**: tries to find suitable parameters for running PeakPickerHiRes on IM-binned data using brute force.

- input*: the name of the raw mzML file, without the file extension
- output*: the directory to output to
- target: the minimum number of features required to be detected in order for a parameter set to be considered suitable

```
python pphr_param_finder.py --input run1/im-pass1-bin1 --output run3 --target 250
```

**feature_match**: tries to find suitable parameters for running im_binning's feature matcher.

- input*: the input featureXML files (see below for the proper format)
- output*: the name of the output directory
- nbins*: the number of bins used (the number of featureXML files per pass)
- nprocs: the number of processors to use in the work pool implementation
- mp_mode: `0` to use the serial implementation, `1` to use the work pool implementation, or `2` to use the process implementation

```
python feature_match.py --input run2/im --output run4 --nbins 10 --nprocs 16 --mp_mode 1
```

**peak_picker**: a custom peak picker for use on MS data with IM information, intended to replace PeakPickerHiRes.

- input*: the name of the raw mzML file without the file extension
- output*: the name of the peak-picked mzML file to output to

```
python peak_picker.py --input sample --output pp_sample
```

**utils/translate_features**: converts features in .csv format to .featureXML format.

- input*: the name of the input csv file without the file extension
- output*: the name of the output featureXML file, without the file extension

```
python translate_features.py --input sample_csv --output sample_fxml
```

**utils/csvFilter**: filters data from a .csv file (which has been converted from a .tsv file). Either filters by initial_peak_quality (q) or extracts RT, m/z, and intensity data; i.e. filter performs a "SELECT * FROM sample_csv WHERE q > 0.01" and extract performs a "SELECT RT, m/z, intensity FROM sample_csv".

Usage: `csvFilter <csv file> <mode>`
- mode*: 0 to filter by q or 1 to extract triplet data
Note: the threshold value can be changed within the source code

**utils/plot_data**: plots 3D data in a scatter plot.

- input*: the name of the input file (each line must have three numbers)
- output: the name of the output file to save the graph to
- xlabel: the label of the x-axis
- ylabel: the label of the y-axis
- zlabel: the label of the z-axis

```
python plot_data.py --input scatter --output graph.png --xlabel RT --ylabel mz --zlabel intensity
```

**utils/tsv_to_csv**: converts a tsv file to a csv file (if a machine does not have Excel to do this).

Usage: `python tsv_to_csv.py < input_tsv > output_csv`

**utils/extract_ms1**: discards all non-MS1 spectra from a raw mzML file.

- input*: the name of the raw mzML file without the file extension
- output*: the name of the new mzML file to output to

```
python extract_ms1.py --input sample --output small_sample
```

**legacy/***: cluster_finder, plane_fitter, and ransac. This is an old approach and is no longer supported or documented.
