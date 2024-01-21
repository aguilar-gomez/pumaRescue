nohup fasterq-dump  --split-files SRR6071635 > out.download2
nohup bwa index GCF_014898765.1_PumYag_genomic.fna > out.index &
nohup perl /home/diana/bin/countfasta.pl GCF_014898765.1_PumYag_genomic.fna > out.countfasta_outgroup &
