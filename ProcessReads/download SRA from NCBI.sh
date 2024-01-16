for accession in SRR7537344 SRR7537345
do
fasterq-dump  --split-files $accession &
done
