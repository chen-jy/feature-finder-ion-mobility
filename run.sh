#!/usr/bin/env bash

RUNNUM=1
NUMBINS=50
PEAKPICK=1

INFILE="../sample"
OUTFILE="rt"
OUTDIR="../run$RUNNUM"

FOUNDF="../run$RUNNUM/$OUTFILE"
OPENMSF="openms_features"

if [[ $PEAKPICK ]]; then
python3 im_binning.py --infile $INFILE --outfile $OUTFILE --outdir $OUTDIR --num_bins $NUMBINS --peak_pick; else
python3 im_binning.py --infile $INFILE --outfile $OUTFILE --outdir $OUTDIR --num_bins $NUMBINS;
fi

python3 cmp_features.py --found $FOUNDF --openms $OPENMSF --out $OUTFILE
