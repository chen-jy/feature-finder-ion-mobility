# FragToPre

**Warning: this README was last updated 2019-09-24 and is probably outdated; it is not updated regularly.**

Approaches to four-dimensional mass spectrometry feature finding (with retention time, mass to charge, intensity, and ion mobility data).

Descriptions of python tools:

**im_binning**: splits raw mzML data into frames by RT, then for each frame, bins data points in the IM dimension. All of the respective bins from each frame are then linked together, and FeatureFinderCentroided (with optional peak picking) is run on each.

- infile*: the name of the source mzML file with the file extension
- outfile*: the name to group the output by (output prefix)
- outdir*: the directory to output to (it must already exist)
- num_bins*: the number of bins to use
- mz_eps*: the m/z epsilon to use when binning
- int_match*: `1` to perform feature matching within individual bins or `0` to skip it
- peak_pick*: `0` to skip peak picking, `1` to use PeakPickerHiRes, or `2` to use the custom peak picker (peak_picker.py)
- min_req*: (custom PP) the minimum number of peaks to either side of a start peak to be consider a peak set
- window_size*: (custom PP) the maximum m/z window to consider for a peak set (`-1` for infinity)
- strict*: (custom PP) `1` for all neighbouring peaks to be strictly non-increasing or `0` to allow a single increase in intensity on either side
- sequential*: (custom PP) `1` to use left-to-right m/z mode or `0` to use decreasing peak intensity mode
- match_only: if `1`, only the feature matcher will be run (it assumes that binning has already completed and all required files are present)

```
python im_binning.py --infile sample.mzML --outfile run --outdir im1 --num_bins 10 --mz_eps 0.001 --int_match 1 --peak_pick 2 --min_req 3 --window_size 0.015 --strict 1 --sequential 1
```

**peak_picker**: a simple custom peak picker for use on MS data with IM information, intended to replace PeakPickerHiRes in the case that it destroys the data.

- input*: the name of the raw mzML file
- output*: the name of the peak-picked mzML file to output to

```
python peak_picker.py --input sample.mzML --output picked_sample.mzML
```

**csv_cmp**: compares a set of features (as a csv or featureXML file) with another set of reference features (as a csv file). In csv-csv mode, the set of reference features is sorted by decreasing intensity, and a `[FOUND]` is appended to each line that represents a common feature. In featureXML-csv mode, features are separated into categories (all features, missing features, unique features, and common features), and are saved as featureXML files.
- input*: the name of the input file (can be csv or featureXML)
- ref*: the name of the csv file containing the reference feature set
- output*: if \<input\> is a csv file, the name of the output csv file. Otherwise, the prefix to use for the output featureXML files

```
python csv_cmp.py --input im1/features.csv --ref ref_features.csv --output im1/common_features.csv
python csv_cmp.py --input im1/run.featureXML --ref ref_features.csv --output im1/features_cmp
```

**baseline**: splits raw mzML data into frames by RT, then for each frame, swaps RT data with IM data. Next, FeatureFinderCentroided is run, and the results are linked together across frames.

- infile*: the name of the source mzML file without the file extension
- outfile*: the name to group the output by (output prefix)
- outdir*: the directory to output to (it must already exist)
- mode*: `0` to find features and link fragments to precursors, `1` to only find features, or `2` to only link fragments to precursors

```
export OMP_NUM_THREADS=1
python baseline.py --infile sample --outfile run --outdir baseline --mode 1
```

**utils/translate_features**: converts features between csv and featureXML formats.

- input*: the name of the input file
- output*: the name of the output file

```
python translate_features.py --input features.csv --output features.featureXML
python translate_features.py --input features.featureXML --output features.csv
```

**utils/extract_ms1**: removes all non-MS1 spectra from a raw mzML file.

- input*: the name of the raw mzML file
- output*: the name of the new mzML file to output to

```
python extract_ms1.py --input sample.mzML --output small_sample.mzML
```

**utils/csvFilter**: filters data from a csv file.
- i*: input csv
- o*: output csv
- m*: mode (see below)
- a: min RT
- b max RT

If $m$ is `0`, a selection of $a \leq$ RT $\leq b \,\wedge$ q_score $\leq$ 0.01 $\wedge$ peak_group_rank = 1 is done.
If $m$ is `1`, a projection of RT, m/z, and (precursor) Intensity of all the lines is done.
If $m$ is `2`, a sort by decreasing precusor intensity is done.

Note: compile with at least `-std=c++11`
```
./csvFilter -i full_sample.csv -o q_filtered.csv -m 0 -a 800 -b 900
./csvFilter -i q_filtered.csv -o q_filtered_small.csv -m 1
./csvFilter -i full_sample.csv -o sorted_sample.csv -m 2
```

**utils/legacy/***: tsv_to_csv and plot_data. These are old tools and are no longer documented (go through the README's history to see old documentation), although they may still be useful in some scenarios.

**legacy/***: binning/\* and clustering/\*. These are old tools and are no longer documented (go through the README's history to see old documentation), although the tools in the former subdirectory may still be useful in some scenarios.
