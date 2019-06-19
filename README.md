# FragToPre
Approaches to linking fragments to their precursors in mass spec data.

Running the algorithms:

This splits frames by retention time, then for each, swaps RT with ion mobility to run FeatureFinderCentroided.
```export OMP_NUM_THREADS=1
python baseline.py --infile sample_3001_to_3098 --outfile rt --outdir run --mz_epsilon 2.5 --im_epsilon 0.025 --num_frames 833 --window_size 17 --rt_length 2
```

This splits frames by retention time, then bins each in the ion mobility dimension to run FeatureFinderCentroided.
```python im_binning.py --infline sample --outfile rt --outdir run
```
