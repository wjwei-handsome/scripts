from pysam import VariantFile
from rich.progress import track


def get_hash_fa(seq_path, total, mode, rm_flag=True):
    with open(seq_path) as f:
        _hash_fa = {}
        for line in track(f, total=total, description=f'Reading {mode} sequences'):
            try:
                line=line.strip()
                bed_string = line.split('\t')[0]
                if 'unk' in bed_string or 'scaf' in bed_string or 'scfd' in bed_string: # TODO: 这里可以改成正则表达式
                    bed_string = bed_string.replace('_','tmp',1)
                if 'Unknown' in bed_string:
                    bed_string = bed_string.replace('_','tmp',2)  ##sb
                if 'scaftig' in bed_string:
                    bed_string = bed_string.replace('_','tmp',4)  ##sb
                chro = '_'.join(bed_string.replace('tmp', '_').split('_')[:-1])
                start = int(bed_string.split('_')[1].split('-')[0])-1
                end = int(bed_string.split('_')[1].split('-')[1].split(':')[0])
                strand = bed_string.split('_')[1].split('-',1)[1].split(':')[1]
                tmp_coord = "@".join([chro,str(start),str(end),'J','N',strand])
                _hash_fa[tmp_coord] = line.split('\t')[1]
            except IndexError:
                print(line.split('\t')[0])
        if rm_flag:
            os.system(f"rm -f {seq_path}")
    return _hash_fa


def write_seq_vcf(outprefix,ref_hash_fa, query_hash_fa, ref_bed_total, query_bed_total, rm_flag=True, leftshift=False):

    vcf_in2 = VariantFile(outprefix+'_bedstring.vcf')
    final_vcf_path = outprefix+'_seq.vcf'
    vcf_out2 = VariantFile(final_vcf_path, 'w', header=vcf_in2.header)

    for rec in track(vcf_in2.fetch(), total=(ref_bed_total+query_bed_total)/2):
        hased_ref_coord = rec.ref
        hased_query_coord = rec.alts[0]
        ref_seq = ref_hash_fa[hased_ref_coord]
        # query_seq = query_hash_fa.get(hased_query_coord, 'N')
        query_seq = query_hash_fa[hased_query_coord]
        rec.alleles = (ref_seq, query_seq)
        if leftshift:
            rec.start = rec.start - 1
        vcf_out2.write(rec)

    vcf_in2.close()
    if rm_flag:
        os.system(f"rm -f {outprefix+'_bedstring.vcf'}")
    vcf_out2.close()

ref_hash_fa = get_hash_fa('Zl-TEO02_ref_hash_fa.tsv', 300000, 'ref', rm_flag=False)
query_hash_fa = get_hash_fa('Zl-TEO02_query_hash_fa.tsv', 300000, 'query', rm_flag=False)

write_seq_vcf('Zl-TEO02',ref_hash_fa, query_hash_fa, 300000, 300000, rm_flag=False, leftshift=True)