vcf=$1
ref=$2
region=$3
wordir=$4

## mkdir workdir
mkdir $wordir
touch $wordir/all.fa
## extract the samples from the vcf
module load BCFtools/1.15.1
bcftools query -l $vcf > $wordir/samples.txt

## extract the reference sequence
module load SAMtools/1.9
samtools faidx $ref $region > $wordir/ref.fa
## extract the reference sequence for each sample
while read sample;
do
seq=$(bcftools consensus -f $wordir/ref.fa -s $sample $vcf -p $sample -H2 -M '-')
echo $seq >> $wordir/all.fa
done < $wordir/samples.txt
