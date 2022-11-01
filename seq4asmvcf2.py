import argparse
from dataclasses import dataclass
import logging
import sys
from pysam import VariantFile
from rich.progress import track


from rich.console import Console
from rich.logging import RichHandler

console = Console()
logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s - %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S',
                    handlers=[RichHandler(rich_tracebacks=True,
                                          console=Console(stderr=True))])

def get_args():
    # define arguments
    parser = argparse.ArgumentParser(
        description=None, prog='seq4asmvcf2.py',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # required arguments
    required_parser = parser.add_argument_group('required arguments')
    required_parser.add_argument('-r', '--ref', action='store',
                                 type=str, dest='ref_fa_path',required=True,
                                 help='ref_fa_path')
    required_parser.add_argument('-q', '--query', action='store', dest='query_fa_path', type=str, required=True,
                                 help='query_fa_path')
    required_parser.add_argument('-v', '--vcf', action='store', dest='vcf', type=str, required=True,
                                 help='vcf file')
    required_parser.add_argument('-b', '--bed', action='store', dest='bed', type=str, required=True,
                                 help='bed file')
    required_parser.add_argument('-o', '--out', action='store', dest='out', type=str, required=True,
                                 help='output prefix')


    # options arguments

    return parser.parse_args()

args = get_args()


@dataclass
class Coord:
    chro: str
    start: int
    end: int
    strand: str

    def __repr__(self) -> str:
        return f"{self.chro};{self.start};{self.end};J;N;{self.strand}"

@dataclass
class RefQueryBed:
    ref_coord: Coord
    query_coord: Coord

    def __repr__(self) -> str:
        return f"{self.ref_coord}\n{self.query_coord}"



def dict_add_rep(dic, add_key, add_value):
    if add_key not in dic:
        dic[add_key] = add_value
    elif isinstance(dic[add_key], list):
        dic[add_key].append(add_value)
    else:
        dic[add_key] = [dic[add_key], add_value]

# REF_FA_PATH = args.ref_fa_path
# REF_FA = pybedtools.BedTool(REF_FA_PATH)

# QUERY_FA_PATH = args.query_fa_path
# QUERY_FA = pybedtools.BedTool(QUERY_FA_PATH)


# a = pybedtools.BedTool("1\t1\t5\tX\tX\t-", from_string=True)
# a = a.sequence(fi=REF_FA,s=True,fullHeader=True)
# print(open(a.seqfn).read())

# def exact_seq_from_bed(bed:Coord, fa:pybedtools.BedTool):
#     pybdt_obj = pybedtools.BedTool("\t".join(
#         [bed.chro, str(bed.start-1), str(bed.end), 'JUNK', 'NAME', bed.strand]
#         ), from_string=True)
#     _seq = pybdt_obj.sequence(fi=fa,s=True,)
#     seq = open(_seq.seqfn).read()
#     seq = seq.strip().split('\n')[1]
#     return seq


# bed_path = "data/B97_scaftig500.AsmSV.Assemblytics_structural_variants.bed"
bed_path = args.bed
bed_file_lines = sum(1 for line in open(bed_path))
with open(bed_path) as bed_file:
    _bed_coord_dict  = {}
    for line in track(bed_file, total=bed_file_lines,description='Reading bed'):
        line = line.strip()
        if line.startswith("#"):
            continue
        ref = line.split("\t")[0]
        ref_start = int(line.split("\t")[1])
        ref_coord = Coord(ref, ref_start-1, int(line.split("\t")[2]), line.split("\t")[5])
        query_coord_str = line.split("\t")[9]
        query_coord = Coord(
            query_coord_str.split(":")[0],
            int(query_coord_str.split(":")[1].split('-')[0])-1,
            int(query_coord_str.split(":")[1].split('-')[1]),
            query_coord_str.split(":")[2],
            )
        ref_query_bed = RefQueryBed(ref_coord, query_coord)
        hased_bed_ID = f"{ref}:{ref_start}"
        dict_add_rep(_bed_coord_dict, hased_bed_ID, ref_query_bed)
    # bed_coord_dict = {}
    # for k,v in _bed_coord_dict.items():
    #     if isinstance(v, list):
    #         for idx, i in enumerate(v):
    #             bed_coord_dict[f"{k}_{idx}"] = i
    #     else:
    #         bed_coord_dict[k] = v
    # del _bed_coord_dict

# vcf_path = 'data/B97_scaftig500.vcf'
vcf_path = args.vcf
vcf_in = VariantFile(vcf_path, 'r')  # auto-detect input format

## 这一步看似很傻逼，实际上它可以修复header
_ = []
for rec in vcf_in.fetch():
    _.append(rec)

bedstring_vcf_out_path = args.out+'_bedstring.vcf'
vcf_out = VariantFile(bedstring_vcf_out_path, 'w', header=vcf_in.header)

ref_bed_file_path = args.out+'_ref.bed'
query_bed_file_path = args.out+'_query.bed'
ref_bed_file = open(ref_bed_file_path, 'w')
query_bed_file = open(query_bed_file_path, 'w')

# count for bar
ref_bed_total = 0
query_bed_total = 0

for rec in track(_,description=f'Writing bedstring vcf to {bedstring_vcf_out_path}'):
    hashed_vcf_ID = f'{rec.chrom}:{rec.pos}'
    ref_query_bed = _bed_coord_dict[hashed_vcf_ID]
    if isinstance(ref_query_bed, list):
        for idx, i in enumerate(ref_query_bed):
            rec.id = f"{hashed_vcf_ID}_{idx}"
            ref_bed = i.ref_coord
            query_bed = i.query_coord
            ref = str(ref_bed)
            ref_bed_file.write(ref.replace(';','\t')+'\n')
            ref_bed_total += 1
            alt = str(query_bed)
            query_bed_file.write(alt.replace(';','\t')+'\n')
            query_bed_total += 1
            rec.alleles = (ref, alt)   ## replace sequences
            vcf_out.write(rec)
    else:
        ref_bed = ref_query_bed.ref_coord
        query_bed = ref_query_bed.query_coord
        # rec2 = rec.copy()
        ref = str(ref_bed)
        ref_bed_file.write(ref.replace(';','\t')+'\n')
        ref_bed_total += 1
        alt = str(query_bed)
        query_bed_file.write(alt.replace(';','\t')+'\n')
        query_bed_total += 1
        rec.alleles = (ref, alt)
        # print(rec)
        vcf_out.write(rec)

vcf_in.close()
vcf_out.close()
ref_bed_file.close()
query_bed_file.close()

logging.info("Done with bedstring vcf")

import os

with console.status("[bold green]Exract seq using seqkit...") as status:
    ref_bed_seq_path = args.out+'_ref_hash_fa.tsv'
    query_bed_seq_path = args.out+'_query_hash_fa.tsv'
    sig = os.system(f"seqkit subseq --bed {ref_bed_file_path} {args.ref_fa_path} | seqkit fx2tab -i > {ref_bed_seq_path}")
    if sig != 0:
        logging.error("seqkit subseq error with ref")
        sys.exit(1)
    os.system(f"rm -f {ref_bed_file_path}")
    sig2 = os.system(f"seqkit subseq --bed {query_bed_file_path} {args.query_fa_path} | seqkit fx2tab -i > {query_bed_seq_path}")
    if sig2 != 0:
        logging.error("seqkit subseq error with query")
        sys.exit(1)
    os.system(f"rm -f {query_bed_file_path}")

logging.info('Done with hashed_fa_bed')


def get_hash_fa(seq_path, total, mode):
    with open(seq_path) as f:
        _hash_fa = {}
        for line in track(f, total=total, description=f'Reading {mode} sequences'):
            try:
                line=line.strip()
                bed_string = line.split('\t')[0]
                if 'unk' in bed_string or 'scfd' in bed_string:
                    bed_string = bed_string.replace('_','tmp',1)
                if 'Unknown' in bed_string:
                    bed_string = bed_string.replace('_','tmp',2)  ##sb
                chro = '_'.join(bed_string.replace('tmp', '_').split('_')[:-1])
                start = int(bed_string.split('_')[1].split('-')[0])-1
                end = int(bed_string.split('_')[1].split('-')[1].split(':')[0])
                strand = bed_string.split('_')[1].split('-')[1].split(':')[1]
                tmp_coord = ";".join([chro,str(start),str(end),'J','N',strand])
                _hash_fa[tmp_coord] = line.split('\t')[1]
            except IndexError:
                print(line.split('\t')[0])
        os.system(f"rm -f {seq_path}")
    return _hash_fa

# with open(ref_bed_seq_path) as f:
#     ref_hash_fa = {}
#     for line in track(f, total=ref_bed_total, description='Reading ref sequences'):
#         line=line.strip()
#         bed_string = line.split('\t')[0]
#         chro = bed_string.split('_')[0]
#         start = int(bed_string.split('_')[1].split('-')[0])-1
#         end = int(bed_string.split('_')[1].split('-')[1].split(':')[0])
#         strand = bed_string.split('_')[1].split('-')[1].split(':')[1]
#         tmp_coord = ";".join([chro,str(start),str(end),'J','N',strand])
#         ref_hash_fa[tmp_coord] = line.split('\t')[1]
#     os.system(f"rm -f {ref_bed_seq_path}")

# with open(query_bed_seq_path) as f:
#     query_hash_fa = {}
#     for line in track(f, total=query_bed_total, description='Reading query suquences'):
#         line=line.strip()
#         bed_string = line.split('\t')[0]
#         chro = bed_string.split('_')[0]
#         start = int(bed_string.split('_')[1].split('-')[0])-1
#         end = int(bed_string.split('_')[1].split('-')[1].split(':')[0])
#         strand = bed_string.split('_')[1].split('-')[1].split(':')[1]
#         tmp_coord = ";".join([chro,str(start),str(end),'J','N',strand])
#         query_hash_fa[tmp_coord] = line.split('\t')[1]

# with console.status("[bold green]format scfd...") as status:
#     os.system(f"sed -i s/scfd_/scfd/g {ref_bed_seq_path}")
#     os.system(f"sed -i s/scfd_/scfd/g {query_bed_seq_path}")
#     os.system(f"sed -i s/scfd_/scfd/g {query_bed_file_path}")
#     os.system(f"sed -i s/scfd_/scfd/g {bedstring_vcf_out_path}")


ref_hash_fa = get_hash_fa(ref_bed_seq_path, ref_bed_total, 'ref')
query_hash_fa = get_hash_fa(query_bed_seq_path, query_bed_total, 'query')


vcf_in2 = VariantFile(bedstring_vcf_out_path)
final_vcf_path = args.out+'_seq.vcf'
vcf_out2 = VariantFile(final_vcf_path, 'w', header=vcf_in2.header)

for rec in track(vcf_in2.fetch(), total=(ref_bed_total+query_bed_total)/2):
    hased_ref_coord = rec.ref
    hased_query_coord = rec.alts[0]
    ref_seq = ref_hash_fa[hased_ref_coord]
    query_seq = query_hash_fa.get(hased_query_coord, 'N')
    # query_seq = query_hash_fa[hased_query_coord]
    rec.alleles = (ref_seq, query_seq)
    vcf_out2.write(rec)

vcf_in2.close()
os.system(f"rm -f {bedstring_vcf_out_path}")
vcf_out2.close()