#!/usr/bin/env bash

# Combined scripts (for three cores):

for ((i=0; i<5; i+=1)); do python im_binning.py --infile /mnt/d/sample.mzML --outfile run --outdir /mnt/d/pp$((i * 3 + 1)) --num_bins 10 --mz_eps 0.001 --int_match 1 --peak_pick 2 --min_req $((i + 1)) --window_size -1 --strict 1 --sequential 1; done; for ((i=1; i<5; i+=1)); do python im_binning.py --infile /mnt/d/sample.mzML --outfile run --outdir /mnt/d/bc$((i * 3 + 1)) --num_bins $((5 * i)) --mz_eps 0.001 --int_match 1 --peak_pick 2 --min_req 1 --window_size 0.015 --strict 1 --sequential 1; done; for ((i=0; i<3; i+=1)); do python im_binning.py --infile /mnt/d/sample.mzML --outfile run --outdir /mnt/d/pp$((i * 3 + 2))b --num_bins 10 --mz_eps 0.001 --int_match 1 --peak_pick 2 --min_req $((i + 1)) --window_size 0.015 --strict 1 --sequential 0; done

for ((i=0; i<5; i+=1)); do python im_binning.py --infile /mnt/d/sample.mzML --outfile run --outdir /mnt/d/pp$((i * 3 + 2)) --num_bins 10 --mz_eps 0.001 --int_match 1 --peak_pick 2 --min_req $((i + 1)) --window_size 0.015 --strict 1 --sequential 1; done; for ((i=1; i<5; i+=1)); do python im_binning.py --infile /mnt/d/sample.mzML --outfile run --outdir /mnt/d/bc$((i * 3 + 2)) --num_bins $((5 * i)) --mz_eps 0.001 --int_match 1 --peak_pick 2 --min_req 2 --window_size 0.015 --strict 1 --sequential 1; done; for ((i=1; i<3; i+=1)); do python im_binning.py --infile /mnt/d/sample.mzML --outfile run --outdir /mnt/d/pp$((i * 3 + 2))a --num_bins 10 --mz_eps 0.001 --int_match 1 --peak_pick 2 --min_req $((i + 1)) --window_size 0.015 --strict 0 --sequential 1; done

for ((i=0; i<5; i+=1)); do python im_binning.py --infile /mnt/d/sample.mzML --outfile run --outdir /mnt/d/pp$((i * 3 + 3)) --num_bins 10 --mz_eps 0.001 --int_match 1 --peak_pick 2 --min_req $((i + 1)) --window_size 0.02 --strict 1 --sequential 1; done; for ((i=1; i<5; i+=1)); do python im_binning.py --infile /mnt/d/sample.mzML --outfile run --outdir /mnt/d/bc$((i * 3 + 3)) --num_bins $((5 * i)) --mz_eps 0.001 --int_match 1 --peak_pick 2 --min_req 3 --window_size 0.015 --strict 1 --sequential 1; done; for ((i=1; i<4; i+=1)); do python im_binning.py --infile /mnt/d/sample.mzML --outfile run --outdir /mnt/d/bc$i --num_bins 1 --mz_eps 0.001 --int_match 1 --peak_pick 2 --min_req $i --window_size 0.015 --strict 1 --sequential 1; done; for ((i=1; i<3; i+=1)); do python im_binning.py --infile /mnt/d/sample.mzML --outfile run --outdir /mnt/d/pp$((i * 3 + 2))c --num_bins 10 --mz_eps 0.001 --int_match 1 --peak_pick 2 --min_req $((i + 1)) --window_size 0.015 --strict 0 --sequential 0; done

# Individual scripts

for ((i=0; i<5; i+=1)); do python im_binning.py --infile /mnt/d/sample.mzML --outfile run --outdir /mnt/d/pp$((i * 3 + 1)) --num_bins 10 --mz_eps 0.001 --int_match 1 --peak_pick 2 --min_req $((i + 1)) --window_size -1 --strict 1 --sequential 1; done

for ((i=0; i<5; i+=1)); do python im_binning.py --infile /mnt/d/sample.mzML --outfile run --outdir /mnt/d/pp$((i * 3 + 2)) --num_bins 10 --mz_eps 0.001 --int_match 1 --peak_pick 2 --min_req $((i + 1)) --window_size 0.015 --strict 1 --sequential 1; done

for ((i=0; i<5; i+=1)); do python im_binning.py --infile /mnt/d/sample.mzML --outfile run --outdir /mnt/d/pp$((i * 3 + 3)) --num_bins 10 --mz_eps 0.001 --int_match 1 --peak_pick 2 --min_req $((i + 1)) --window_size 0.02 --strict 1 --sequential 1; done

#

for ((i=1; i<5; i+=1)); do python im_binning.py --infile /mnt/d/sample.mzML --outfile run --outdir /mnt/d/bc$((i * 3 + 1)) --num_bins $((5 * i)) --mz_eps 0.001 --int_match 1 --peak_pick 2 --min_req 1 --window_size 0.015 --strict 1 --sequential 1; done

for ((i=1; i<5; i+=1)); do python im_binning.py --infile /mnt/d/sample.mzML --outfile run --outdir /mnt/d/bc$((i * 3 + 2)) --num_bins $((5 * i)) --mz_eps 0.001 --int_match 1 --peak_pick 2 --min_req 2 --window_size 0.015 --strict 1 --sequential 1; done

for ((i=1; i<5; i+=1)); do python im_binning.py --infile /mnt/d/sample.mzML --outfile run --outdir /mnt/d/bc$((i * 3 + 3)) --num_bins $((5 * i)) --mz_eps 0.001 --int_match 1 --peak_pick 2 --min_req 3 --window_size 0.015 --strict 1 --sequential 1; done

for ((i=1; i<4; i+=1)); do python im_binning.py --infile /mnt/d/sample.mzML --outfile run --outdir /mnt/d/bc$i --num_bins 1 --mz_eps 0.001 --int_match 1 --peak_pick 2 --min_req $i --window_size 0.015 --strict 1 --sequential 1; done

#

for ((i=0; i<3; i+=1)); do python im_binning.py --infile /mnt/d/sample.mzML --outfile run --outdir /mnt/d/pp$((i * 3 + 2))b --num_bins 10 --mz_eps 0.001 --int_match 1 --peak_pick 2 --min_req $((i + 1)) --window_size 0.015 --strict 1 --sequential 0; done

for ((i=1; i<3; i+=1)); do python im_binning.py --infile /mnt/d/sample.mzML --outfile run --outdir /mnt/d/pp$((i * 3 + 2))a --num_bins 10 --mz_eps 0.001 --int_match 1 --peak_pick 2 --min_req $((i + 1)) --window_size 0.015 --strict 0 --sequential 1; done

for ((i=1; i<3; i+=1)); do python im_binning.py --infile /mnt/d/sample.mzML --outfile run --outdir /mnt/d/pp$((i * 3 + 2))c --num_bins 10 --mz_eps 0.001 --int_match 1 --peak_pick 2 --min_req $((i + 1)) --window_size 0.015 --strict 0 --sequential 0; done

# Comparison scripts

for ((i=1; i<16; i+=1)); do python csv_cmp.py --input /mnt/d/pp$i/run.featureXML --ref /mnt/d/sample.csv --output /mnt/d/pp$i/features; done

for i in 2b 5a 5b 5c 8a 8b 8c; do python csv_cmp.py --input /mnt/d/pp$i/run.featureXML --ref /mnt/d/sample.csv --output /mnt/d/pp$i/features; done

for ((i=1; i<4; i+=1)); do python csv_cmp.py --input /mnt/d/bc$i/run-pass1-bin0.featureXML --ref /mnt/d/sample.csv --output /mnt/d/bc$i/features; done

for ((i=4; i<16; i+=1)); do python csv_cmp.py --input /mnt/d/bc$i/run.featureXML --ref /mnt/d/sample.csv --output /mnt/d/bc$i/features; done

#

for ((i=1; i<16; i+=1)); do python csv_cmp.py --input /mnt/d/pp$i/run.featureXML --ref /mnt/d/sample_less.csv --output /mnt/d/pp$i/features_less; done; for i in 2b 5a 5b 5c 8a 8b 8c; do python csv_cmp.py --input /mnt/d/pp$i/run.featureXML --ref /mnt/d/sample_less.csv --output /mnt/d/pp$i/features_less; done; for ((i=1; i<4; i+=1)); do python csv_cmp.py --input /mnt/d/bc$i/run-pass1-bin0.featureXML --ref /mnt/d/sample_less.csv --output /mnt/d/bc$i/features_less; done; for ((i=4; i<16; i+=1)); do python csv_cmp.py --input /mnt/d/bc$i/run.featureXML --ref /mnt/d/sample_less.csv --output /mnt/d/bc$i/features_less; done

# Translation scripts

for ((i=1; i<16; i+=1)); do python translate_features.py --input /mnt/d/pp$i/run.featureXML --output /mnt/d/pp$i/run.csv; done; for i in 2b 5a 5b 5c 8a 8b 8c; do python translate_features.py --input /mnt/d/pp$i/run.featureXML --output /mnt/d/pp$i/run.csv; done; for ((i=1; i<4; i+=1)); do python translate_features.py --input /mnt/d/bc$i/run-pass1-bin0.featureXML --output /mnt/d/bc$i/run.csv; done; for ((i=4; i<16; i+=1)); do python translate_features.py --input /mnt/d/bc$i/run.featureXML --output /mnt/d/bc$i/run.csv; done

# More comparison scripts

for ((i=1; i<16; i+=1)); do python csv_cmp.py --input /mnt/d/pp$i/run.csv --ref /mnt/d/sample.csv --output /mnt/d/pp$i/features.csv; done; for i in 2b 5a 5b 5c 8a 8b 8c; do python csv_cmp.py --input /mnt/d/pp$i/run.csv --ref /mnt/d/sample.csv --output /mnt/d/pp$i/features.csv; done; for ((i=1; i<4; i+=1)); do python csv_cmp.py --input /mnt/d/bc$i/run.csv --ref /mnt/d/sample.csv --output /mnt/d/bc$i/features.csv; done; for ((i=4; i<16; i+=1)); do python csv_cmp.py --input /mnt/d/bc$i/run.csv --ref /mnt/d/sample.csv --output /mnt/d/bc$i/features.csv; done

for ((i=1; i<16; i+=1)); do python csv_cmp.py --input /mnt/d/pp$i/run.csv --ref /mnt/d/sample_less.csv --output /mnt/d/pp$i/features_less.csv; done; for i in 2b 5a 5b 5c 8a 8b 8c; do python csv_cmp.py --input /mnt/d/pp$i/run.csv --ref /mnt/d/sample_less.csv --output /mnt/d/pp$i/features_less.csv; done; for ((i=1; i<4; i+=1)); do python csv_cmp.py --input /mnt/d/bc$i/run.csv --ref /mnt/d/sample_less.csv --output /mnt/d/bc$i/features_less.csv; done; for ((i=4; i<16; i+=1)); do python csv_cmp.py --input /mnt/d/bc$i/run.csv --ref /mnt/d/sample_less.csv --output /mnt/d/bc$i/features_less.csv; done
