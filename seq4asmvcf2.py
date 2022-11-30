#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
'''
@文件        :seq4asmvcf2.py
@说明        :
@时间        :2022/11/10 01:10:03
@作者        :wjwei
@版本        :0.01
@邮箱        :wjwei9908@gmail.com
'''



import argparse
from dataclasses import dataclass
import logging
import sys
import os
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
                                 help='formatted vcf by ..')
    required_parser.add_argument('-b', '--bed', action='store', dest='bed', type=str, required=True,
                                 help='bedpe file by ..')
    required_parser.add_argument('-o', '--out', action='store', dest='out', type=str, required=True,
                                 help='output prefix')


    # options arguments
    optional_parser = parser.add_argument_group('optional arguments')
    optional_parser.add_argument('-k', '--keep', action='store_true', dest='keep', help='if keep middle file',
                                 default=False)
    optional_parser.add_argument('-l', '--leftshift', action='store_true', dest='leftshift', help='if alignt by minimap',
                                 default=False)



    return parser.parse_args()


@dataclass
class Coord:
    chro: str
    start: str
    end: str
    strand: str

    def __repr__(self) -> str:
        return "@".join([self.chro, self.start, self.end, 'J', 'N',self.strand])
        # return f"{self.chro}@{self.start}@{self.end}@J@N@{self.strand}"


    def left_shift(self):
        self.start = str(int(self.start) - 1)
        self.end = str(int(self.end) - 1)

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



def get_bed_coord_dict(bed_path: str):
    bed_file_lines = sum(1 for line in open(bed_path))
    with open(bed_path) as bed_file:
        _bed_coord_dict  = {}
        for line in track(bed_file, total=bed_file_lines,description='Reading bedpe file'):
            line = line.strip()
            if line.startswith("#"):
                continue
            # ref = line.split("\t")[0]
            # ref_start = int(line.split("\t")[1])
            # ref_coord = Coord(ref, ref_start, int(line.split("\t")[2]), line.split("\t")[8])

            # # query_coord_str = line.split("\t")[9]
            # # ref_gap_size = int(line.split("\t")[7])

            # # query_start = int(query_coord_str.split(":")[1].split('-')[0])
            # # if ref_gap_size < 0:
            # #     query_start = int(query_coord_str.split(":")[1].split('-')[0]) + ref_gap_size
            # query_coord_str = line.split("\t")[6]
            # query_coord = Coord(
            #     query_coord_str.split("@")[1],
            #     query_coord_str.split("@")[2],
            #     int(query_coord_str.split(":")[1].split('-')[1]),
            #     query_coord_str.split(":")[],
            #     )
            # ref_query_bed = RefQueryBed(ref_coord, query_coord)
            # # hased_bed_ID = f"{ref}:{ref_start}"
            hased_bed_ID = line.split("\t")[6]
            # hased_bed_ID: INS@1@218831@218839@+@chr1@209816@210556@+
            ref_coord = Coord(
                *hased_bed_ID.split("@")[1:5]
            )
            query_coord = Coord(
                *hased_bed_ID.split("@")[5:9]
            )
            ref_query_bed = RefQueryBed(ref_coord, query_coord)
            _bed_coord_dict[hased_bed_ID] = ref_query_bed
            # dict_add_rep(_bed_coord_dict, hased_bed_ID, ref_query_bed)
    return _bed_coord_dict





def fix_input_vcf_header(VariantFile):
    ## 这一步看似很傻逼，实际上它可以修复header
    _ = []
    for rec in VariantFile.fetch():
        _.append(rec)
    return _

def write_bed_vcf(outprefix, input_vcf_records,bed_coord_dict, leftshift: bool=False):

    bedstring_vcf_out_path = outprefix+'_bedstring.vcf'
    vcf_out = VariantFile(bedstring_vcf_out_path, 'w', header=vcf_in.header)

    ref_bed_file_path = outprefix+'_ref.bed'
    query_bed_file_path = outprefix+'_query.bed'
    ref_bed_file = open(ref_bed_file_path, 'w')
    query_bed_file = open(query_bed_file_path, 'w')

    # count for bar
    ref_bed_total = 0
    query_bed_total = 0

    for rec in track(input_vcf_records,description=f'Writing bedstring vcf to {bedstring_vcf_out_path}'):
        hashed_vcf_ID = rec.id
        ref_query_bed = bed_coord_dict[hashed_vcf_ID]
        # if isinstance(ref_query_bed, list):
        #     for idx, i in enumerate(ref_query_bed):
        #         rec.id = f"{hashed_vcf_ID}_{idx}"
        #         ref_bed = i.ref_coord
        #         query_bed = i.query_coord
        #         ref = str(ref_bed)
        #         ref_bed_file.write(ref.replace(';','\t')+'\n')
        #         ref_bed_total += 1
        #         alt = str(query_bed)
        #         query_bed_file.write(alt.replace(';','\t')+'\n')
        #         query_bed_total += 1
        #         rec.alleles = (ref, alt)   ## replace sequences
        #         vcf_out.write(rec)
        # else:
        ref_bed = ref_query_bed.ref_coord
        query_bed = ref_query_bed.query_coord
        if leftshift:
            ref_bed.left_shift()
            query_bed.left_shift()
        ref = str(ref_bed)
        ref_bed_file.write(ref.replace('@','\t')+'\n')

        ref_bed_total += 1
        alt = str(query_bed)
        query_bed_file.write(alt.replace('@','\t')+'\n')
        query_bed_total += 1
        rec.alleles = (ref, alt)
        vcf_out.write(rec)

    vcf_in.close()
    vcf_out.close()
    ref_bed_file.close()
    query_bed_file.close()

    logging.info("Done with bedstring vcf")

    return ref_bed_file_path, query_bed_file_path, ref_bed_total, query_bed_total


def get_seq_tab(outprefix, ref_bed_file_path, query_bed_file_path,ref_fa_path,query_fa_path, rm_flag=True):
    with console.status("[bold green]Exract seq using seqkit...") as status:
        ref_bed_seq_path = outprefix+'_ref_hash_fa.tsv'
        query_bed_seq_path = outprefix+'_query_hash_fa.tsv'
        sig = os.system(f"seqkit subseq --bed {ref_bed_file_path} {ref_fa_path} | seqkit fx2tab -i > {ref_bed_seq_path}")
        if sig != 0:
            logging.error("seqkit subseq error with ref")
            sys.exit(1)
        sig2 = os.system(f"seqkit subseq --bed {query_bed_file_path} {query_fa_path} | seqkit fx2tab -i > {query_bed_seq_path}")
        if sig2 != 0:
            logging.error("seqkit subseq error with query")
            sys.exit(1)
        if rm_flag:
            os.system(f"rm -f {query_bed_file_path}")
            os.system(f"rm -f {ref_bed_file_path}")

    logging.info('Done with hashed_fa_bed')
    return ref_bed_seq_path, query_bed_seq_path


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


if __name__ == '__main__':
    args = get_args()
    if_keep = args.keep
    print(if_keep)
    leftshift = args.leftshift
    print(leftshift)
    bed_coord_dict = get_bed_coord_dict(args.bed)
    vcf_in = VariantFile(args.vcf, 'r')  # auto-detect input format
    vcf_records = fix_input_vcf_header(vcf_in)
    outprefix = args.out
    ref_bed_file_path, query_bed_file_path, ref_bed_total, query_bed_total = write_bed_vcf(outprefix, vcf_records, bed_coord_dict, leftshift)
    ref_bed_seq_path, query_bed_seq_path = get_seq_tab(outprefix, ref_bed_file_path, query_bed_file_path, args.ref_fa_path, args.query_fa_path,rm_flag=if_keep)
    ref_hash_fa = get_hash_fa(ref_bed_seq_path, ref_bed_total, 'ref', rm_flag=if_keep)
    query_hash_fa = get_hash_fa(query_bed_seq_path, query_bed_total, 'query', rm_flag=if_keep)
    write_seq_vcf(outprefix, ref_hash_fa, query_hash_fa, ref_bed_total, query_bed_total, rm_flag=if_keep,leftshift=leftshift)