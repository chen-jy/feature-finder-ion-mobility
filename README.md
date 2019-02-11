# FragToPre
Approaches to linking fragments to their precursors in mass spec data.

Sample Run Command: 

export OMP_NUM_THREADS=1

python baseline.py --infile sample_3001_to_3098 --outfile rt --outdir run --mz_epsilon 2.5 --im_epsilon 0.025 --num_frames 833 --window_size 17 --rt_length 2
