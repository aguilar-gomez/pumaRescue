

#i input_file
#s sample_file
#a overall ancestry proportion (0.8 for FL, 0.2 for TX; value coming from OHANA)
#-p 0 100000 0.8    pulse 0 for FL, enter long time ago (100000 >> 10), 0.8 proportion
#-p 1 -5 0.2       pulse 1 for TX, estimated to enter 5 gen ago (negative value indicating estimation), 0.2 proportion
# -e 0.01 error rate per site
# --tmax maximum time of a pulse
# --pmax maximum portion of ancestry from one ancestry type
# -b 10 1000 10 bootstrap replicates each using a block size of 1000 SNPs (same as default)


ancestry_hmm -i ahmm_input \
             -s sample_ploidy4ahmm  \
             -a 2 0.80 0.20 \
             -p 0 100000 0.8 \
             -p 1 -10 0.2 \
             -e 0.01 \
             --tmax 5 \
             --pmax 1 \
             -b 10 1000
