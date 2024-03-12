conda create -n snpeff -c bioconda snpeff=5.1


conda activate snpeff


FIRSTDIR=$(dirname $(dirname $(which snpEff)))
TMP=$(dirname $(readlink -s `which snpEff`))

SNPEFFDIR=${FIRSTDIR}${TMP/../}


mkdir snpeff
cp $SNPEFFDIR/snpEff.config snpeff/


echo "
# Puma yagouaroundi genome, GCF_014898765.1_PumYag_genomic
PumYag.genome : PumaYag
" >> snpeff/snpEff.config


cp -s /space/s1/lin.yuan/puma/genome_outgroup/GCF_014898765.1_PumYag_genomic.fna ./PumYag.genome
mkdir data
cd data
mkdir genomes
cd genomes
cp -s /space/s1/lin.yuan/puma/genome_outgroup/GCF_014898765.1_PumYag_genomic.fna ./PumYag.fa


mkdir PumYag
cd PumYag
cp -s /space/s1/lin.yuan/puma/genome_outgroup/annotation/genomic.gff ./genes.gff


nohup snpEff build -Xmx4g  -noCheckCds -noCheckProtein -gff3 -c ./snpEff.config  -v PumYag > out.builddb &
