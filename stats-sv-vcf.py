import logging
import sys

from rich.progress import track
from rich.console import Console
from rich.logging import RichHandler

console = Console()
logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s - %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S',
                    handlers=[RichHandler(rich_tracebacks=True,
                                          console=Console(stderr=True))])


file_in = sys.argv[1]
file_out = sys.argv[2]

def get_gt_code(sample_geno):
    gt_code = []
    for i in sample_geno:
        gt_code.append(i.split('=')[1].split('/')[0])
    return "".join(gt_code)

def get_svlen(pos,svtype,end,seq):
    if svtype == "DEL":
        svlen = int(pos) - int(end)
    else:
        svlen = len(seq)
    return svlen

with open(file_in) as f:
    with open(file_out, 'a') as fout:
        for line in track(f, total=2176149):
            line=line.strip()
            chro,pos,svtype,end,seq = line.split("\t")[0:5]
            sample_geno = line.split("\t")[6:]
            gt_code = get_gt_code(sample_geno)
            svlen = get_svlen(pos,svtype,end,seq)
            fout.write(f"{chro}\t{pos}\t{svtype}\t{svlen}\t{gt_code}\n")
