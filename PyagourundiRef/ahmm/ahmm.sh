

#i input format 

ancestry_hmm -i ahmm_input \
             -s sample_ploidy4ahmm  \
             -a 2 0.80 0.20 \
             -p 0 100000 0.8 \
             -p 1 -10 0.2 \
             -e 0.01 \
             --tmax 5 \
             --pmax 1 \
             -b 10 1000
