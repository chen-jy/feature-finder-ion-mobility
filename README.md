# LC-IMS-MS/MS Feature Finder

An approach to four-dimensional mass spectrometry feature finding (retention time, mass-to-charge, intensity, and ion mobility).

*Warning: the modules here have only been tested with pyOpenMS 2.4.0 (and OpenMS 2.4.0). In addition, analyzing full-length mass spectrometry runs (120 min/~12 GiB files) may take several days to complete, depending on the target machine.*

Descriptions of python modules (run each with --help for individual arguments):

**feature_finder_im**: the LC-IMS-MS/MS feature finder. Works by binning each RT spectrum in the raw data by IM, linking all the respective binned spectra together, running an OpenMS feature finder (with optional noise filtering and peak picking) on each bin, then matching all the resulting features together.
```
python feature_finder_im.py --in sample.mzML --out features.featureXML --dir run --num_bins 50 --pp_type pphr --ff_type centroided
```

**peak_picker_im**: a simple custom peak picker for use on MS data containing IM information. For comparison purposes with PeakPickerHiRes.
```
python peak_picker_im.py --in sample.mzML --out sample_picked.mzML
```

**compare_features**: an evaulation utility for the feature finder. Compares a set of found features with a set of reference features. Useful for comparing against converted MaxQuant output.
```
python compare_features.py --in run/features.featureXML --ref evidence.csv --out cmp/evidence
```

**baseline**: a different approach to feature finding (with development currently on hold). Works by splitting raw mzML data into frames by RT, swapping RT and IM data, running FeatureFinderCentroided, and linking the results together across frames. For comparison purposes with feature_finder_im.
```
export OMP_NUM_THREADS=1
python baseline.py --infile sample --outfile baserun --outdir baseline --mode 1
```
