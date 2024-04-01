
#GATK
conda create -n gatk_env
conda activate gatk_env
conda install -c bioconda gatk4=4.2.2

#RepeatMasker
conda create -n repeatmask
conda activate repeatmask
conda install -c bioconda -c conda-forge repeatmasker


