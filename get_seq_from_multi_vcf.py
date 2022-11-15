import argparse
import logging
import sys
import os
from rich.progress import track
from pysam import VariantFile

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
        description=None, prog='get_seq_from_multi_vcf.py',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # required arguments
    required_parser = parser.add_argument_group('required arguments')
    required_parser.add_argument('-i', '--input', action='store', dest='inputvcf', type=str, required=True,)
    required_parser.add_argument('-a', '--all', action='store', dest='allvcf', type=str, required=True,)
    required_parser.add_argument('-o', '--output', action='store', dest='output', type=str, required=True,)
    required_parser.add_argument('-m', '--min', action='store', dest='min', type=int, required=True,)

    # options arguments
    optional_parser = parser.add_argument_group('optional arguments')
    optional_parser.add_argument('-f', '--force', action='store_true', dest='overwrite', help='if overwrite output file',default=False)

    return parser.parse_args()

def check_bcftools() -> bool:

    bcftools_path = os.popen("which bcftools").read().strip()
    if bcftools_path:
        return True
    else:
        return False



def store_all_id_seq(all_vcf_path: str) -> dict:

    vcf = VariantFile(all_vcf_path)
    id_seq_dict = {}
    for rec in track(vcf.fetch(),total=130533, description="Reading ALL VCF..."):
        id_seq_dict[rec.id] = f"{rec.ref}\t{rec.alts[0]}"
    return id_seq_dict


# def get_ref_alt_seq(rcd_id: str)-> str:

    # sample_name = rcd_id.split("@")[5].split('.')[0]
    # bcftools_result = os.popen(f"bcftools view -i '%ID==\"{rcd_id}\"' {sample_name}_50M.vcf |bcftools query -f '%REF\t%ALT\n' ").read().strip()
    # ref_seq, alt_seq = bcftools_result.split("\t")
    # return f"{sample_name}\t{ref_seq}\t{alt_seq}\n"

def check_output(output: str, overwrite: bool):
    if os.path.exists(output):
        if overwrite:
            os.remove(output)
        else:
            console.print(f"{output} exists, please check")
            sys.exit(1)

def get_ref_alt_seq(inputvcf: str, id_seq_dict: dict, output: str, min_sample: int):
    vcf = VariantFile(inputvcf)
    for rec in track(vcf.fetch(),total=68082, description="Processing..."):
        SVTYPE = rec.info["SVTYPE"]
        if SVTYPE == "INS":
            IDLIST = rec.info['IDLIST']
            if len(IDLIST) >= min_sample:
                with open(output, 'a') as f:
                    id_seq = "\n".join([f"{id.split('@')[5]}\t{id_seq_dict[id]}" for id in IDLIST])
                    f.write(f"{rec.id}\t{rec.ref}\t{rec.alts[0]}\t{len(IDLIST)}\n{id_seq}\n")
            else:
                continue
        else:
            continue
def main():
    args = get_args()
    # # check bcftools
    # if not check_bcftools():
    #     console.print("bcftools not found, please install bcftools first")
    #     sys.exit(1)

    ## check all_vcf if exists
    all_vcf_path = args.allvcf
    if not os.path.exists(all_vcf_path):
        console.print(f"{all_vcf_path} not found, please check")
        sys.exit(1)

    ## check input_vcf if exists
    input_vcf = args.inputvcf
    if not os.path.exists(input_vcf):
        console.print(f"{input_vcf} not found, please check")
        sys.exit(1)

    ## check output
    check_output(args.output, args.overwrite)

    all_id_seq_dict = store_all_id_seq(all_vcf_path)
    get_ref_alt_seq(input_vcf, all_id_seq_dict, args.output, args.min)


if __name__ == '__main__':
    main()