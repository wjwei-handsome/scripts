#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
'''
@文件        :combine_bed.py
@说明        :
@时间        :2022/11/30 16:21:17
@作者        :wjwei
@版本        :0.01
@邮箱        :wjwei9908@gmail.com
'''



import argparse
from dataclasses import dataclass
import logging
import json
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
        description=None, prog='combine_bed.py',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # required arguments
    required_parser = parser.add_argument_group('required arguments')
    required_parser.add_argument('-l', '--list', action='store',
                                 type=str, dest='file_list',required=True,
                                 help='bed file path list')
    required_parser.add_argument('-o', '--out', action='store', dest='out', type=str, required=True,
                                 help='output file')


    # options arguments
    # optional_parser = parser.add_argument_group('optional arguments')
    # optional_parser.add_argument('-k', '--keep', action='store_true', dest='keep', help='if keep middle file',
    #                              default=False)
    # optional_parser.add_argument('-l', '--leftshift', action='store_true', dest='leftshift', help='if alignt by minimap',
    #                              default=False)



    return parser.parse_args()


def main(file_list_path: str, out_path: str):
    with open(file_list_path) as f:
        file_list = f.read().splitlines()

    fout = open(out_path, 'a')
    for file in track(file_list):
        phone = file.split('.')[0]
        with open(file) as f:
            tmp_list = []
            for line in track(f):
                if line.startswith('chr'):
                    region_str = line.split("\t")[0]+":"+line.split("\t")[1]+"-"+line.split("\t")[2]
                    try:
                        lgp_float = float(line.split('\t')[3])
                    except:
                        lgp_float = 0
                    tmp_list
                else:
                    continue
            fout.write(phone+"\t"+region_str+"\n")

main(get_args().file_list, get_args().out)

"""
cat tmp|while read l; do q=$(echo $l|cut -f1 -d ':' ); phe=$(echo $l|cut -f2 -d ':'); aa=$(grep '^chr' Flower-2011GX-TSI.max.bed|cut -f4 |sed ":a;N;s/\n/;/g;ta");echo -e "${phe}\t${aa}" >> fff; done
"""