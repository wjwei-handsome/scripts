sample=$1
prefix=$2
cat <(echo -e "id\tcall_gt\tdp\tad") <(bcftools view -f PASS,lowad ${sample}.${prefix}.vcf|bcftools query -f '%ID\t[%GT\t%DP\t%AD]\n') > ${sample}.${prefix}.tsv
cat <(echo -e "id\treal_gt") <(bcftools query -f '%ID\t[%GT]\n' -s ${sample} ../wwj_chr10_1st20M.vcf.gz)|sed 's/|/\//g' > ${sample}.real.tsv

csvtk join -t -f '1,1' ${sample}.real.tsv ${sample}.geno.tsv > ${sample}.compare.tsv

csvtk filter2 -t -f '$real_gt=="1/1"' ${sample}.compare.tsv |wc -l
csvtk filter2 -t -f '$real_gt=="1/1" && $call_gt=="1/1"' ${sample}.compare.tsv |wc -l
csvtk filter2 -t -f '$real_gt=="0/0"' ${sample}.compare.tsv |wc -l
csvtk filter2 -t -f '$real_gt=="0/0"&& $call_gt=="0/0"' ${sample}.compare.tsv |wc -l
