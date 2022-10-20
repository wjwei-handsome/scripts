import argparse
from dataclasses import dataclass
from pysam import VariantFile
import pybedtools
from rich.progress import track


def get_args():
    # define arguments
    parser = argparse.ArgumentParser(
        description=None, prog='seq4asmvcf.py',
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
                                 help='output vcf file')


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
        return f"{self.chro}:{self.start}-{self.end}({self.strand})"

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

REF_FA_PATH = args.ref_fa_path
REF_FA = pybedtools.BedTool(REF_FA_PATH)

QUERY_FA_PATH = args.query_fa_path
QUERY_FA = pybedtools.BedTool(QUERY_FA_PATH)


# a = pybedtools.BedTool("1\t1\t5\tX\tX\t-", from_string=True)
# a = a.sequence(fi=REF_FA,s=True,fullHeader=True)
# print(open(a.seqfn).read())

def exact_seq_from_bed(bed:Coord, fa:pybedtools.BedTool):
    pybdt_obj = pybedtools.BedTool("\t".join(
        [bed.chro, str(bed.start-1), str(bed.end), 'JUNK', 'NAME', bed.strand]
        ), from_string=True)
    _seq = pybdt_obj.sequence(fi=fa,s=True,)
    seq = open(_seq.seqfn).read()
    seq = seq.strip().split('\n')[1]
    return seq


# bed_path = "data/B97_scaftig500.AsmSV.Assemblytics_structural_variants.bed"
bed_path = args.bed
bed_file_lines = sum(1 for line in open(bed_path))
with open(bed_path) as bed_file:
    _bed_coord_dict  = {}
    for line in track(bed_file, total=bed_file_lines):
        line = line.strip()
        if line.startswith("#"):
            continue
        ref = line.split("\t")[0]
        ref_start = int(line.split("\t")[1])
        ref_coord = Coord(ref, ref_start, int(line.split("\t")[2]), line.split("\t")[5])
        query_coord_str = line.split("\t")[9]
        query_coord = Coord(
            query_coord_str.split(":")[0],
            int(query_coord_str.split(":")[1].split('-')[0]),
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

vcf_out = VariantFile(args.out, 'w', header=vcf_in.header)

for rec in track(_):
    hashed_vcf_ID = f'{rec.chrom}:{rec.pos}'
    ref_query_bed = _bed_coord_dict[hashed_vcf_ID]
    if isinstance(ref_query_bed, list):
        for idx, i in enumerate(ref_query_bed):
            rec.id = f"{hashed_vcf_ID}_{idx}"
            ref_bed = i.ref_coord
            query_bed = i.query_coord
            ref = exact_seq_from_bed(ref_bed, REF_FA)
            alt = exact_seq_from_bed(query_bed, QUERY_FA)
            rec.alleles = (ref, alt)   ## replace sequences
            vcf_out.write(rec)
    else:
        ref_bed = ref_query_bed.ref_coord
        query_bed = ref_query_bed.query_coord
        # rec2 = rec.copy()
        ref = exact_seq_from_bed(ref_bed, REF_FA)
        alt = exact_seq_from_bed(query_bed, QUERY_FA)
        rec.alleles = (ref, alt)
        # print(rec)
        vcf_out.write(rec)

vcf_in.close()
vcf_out.close()
