import sys
from pysam import VariantFile
from hashlib import sha1

vcf_path = sys.argv[1]
out_path = sys.argv[2]

vcf_file = VariantFile(vcf_path)

with open(out_path, 'w'):
    for record in vcf_file.fetch():
        snp_id = record.id
        snp_pos = record.pos
        snp_chr = record.chrom
        snp_ref = record.ref
        snp_alt = record.alts[0] ## skip bi-allle
        snp_str = snp_chr+'\n'+snp_pos+'\n'+snp_ref+'\n'+snp_alt+'\n'
        snp_sha1 = sha1(snp_str.encode('utf-8')).hexdigest()
        out_path.write(snp_id+'\t'+snp_sha1+'\n')