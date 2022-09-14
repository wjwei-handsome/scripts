vcf=$0
ref=$1
region=$2
wordir=$3

## mkdir workdir
mkdir $wordir

## extract the samples from the vcf
module load bcftools
bcftools query -l $vcf > $wordir/samples.txt

## extract the reference sequence
module load SAMtools/1.9
samtools faidx $ref $region > $wordir/ref.fa
## extract the reference sequence for each sample
while read sample;
do
bcftools consensus -f $wordir/ref.fa -s $sample $vcf -H1> $wordir/$sample.fa
done < $wordir/samples.txt

cat $wordir/*.fa > $wordir/all.fa