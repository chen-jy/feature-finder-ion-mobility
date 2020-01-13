# FragToPre

**Warning: this README was last updated 2020-01-06 and is probably outdated; it is not updated regularly.**

Approaches to four-dimensional mass spectrometry feature finding (with retention time, mass to charge, intensity, and ion mobility data).

Descriptions of python modules (run each with --help for individual arguments):

**feature_finder_im**: the LC-IMS-MS/MS feature finder. Works by binning each RT spectrum in the raw data by IM, linking all the respective binned spectra together, running an OpenMS feature finder (with optional noise filtering and peak picking) on each bin, then matching all the resulting features together.
```
python feature_finder_im.py --in sample.mzML --out features.featureXML --dir runs/1 --num_bins 10 --pp_type pphr --ff_type centroided
```

**peak_picker_im**: a simple custom peak picker for use on MS data containing IM information. For comparison purposes with PeakPickerHiRes.
```
python peak_picker_im.py --in sample.mzML --out picked_sample.mzML
```

**compare_features**: an evaulation utility for the feature finder. Compares a set of found features with another set of reference features. Useful for comparing against MaxQuant output.
```
python compare_features.py --in found.featureXML --ref evidence.csv --out cmp/1/evidence
```

**baseline**: a different approach to feature finding (development currently on-hold). Works by splitting raw mzML data into frames by RT, swapping RT and IM data, running FeatureFinderCentroided, and linking the results together across frames. For comparison purposes with feature_finder_im.
```
export OMP_NUM_THREADS=1
python baseline.py --infile sample --outfile run --outdir baseline --mode 1
```

