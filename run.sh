#!/usr/bin/env bash

RUNNUM=14
NUMBINS=10
PEAKPICK=

INFILE="../sample"
OUTFILE="rt"
OUTDIR="../run$RUNNUM"
BLDIR="../runb"

FOUNDF="../run$RUNNUM/$OUTFILE-pass1"
OPENMSF="openms_features"
BASELNF="$BLDIR/baseline"

if [[ $PEAKPICK ]]; then
	python3 im_binning.py --infile $INFILE --outfile $OUTFILE --outdir $OUTDIR --num_bins $NUMBINS --peak_pick
else
	python3 im_binning.py --infile $INFILE --outfile $OUTFILE --outdir $OUTDIR --num_bins $NUMBINS;
fi

python3 baseline.py --infile $INFILE --outfile $OUTFILE --outdir $BLDIR --mode 1 --mz_epsilon 2.5 --im_epsilon 0.025 --num_frames 833 --window_size 17 --rt_length 2

if [ -e $OPENMSF.featureXML ]; then
	python3 cmp_features.py --found $FOUNDF --openms $OPENMSF --baseline $BASELNF --out $OUTFILE
else
	python3 cmp_features.py --found $FOUNDF --source $INFILE --baseline $BASELNF --out $OUTFILE
fi
