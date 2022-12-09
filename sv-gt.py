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

def judge(sample_geno):
    ok_list = []
    for item in sample_geno:
        cato = item.split('=')[0].split('-')[0]
        if "=1" in item:
            ok_list.append(cato)
    okk_list = list(set(ok_list))
    if len(okk_list) == 2 and 'Zm' in okk_list:
        return True

with open(file_in) as f:
    with open(file_out, 'a') as fout:
        for line in track(f, total=2176149):
            line=line.strip()
            sv_id = line.split("\t")[0]
            sample_geno = line.split("\t")[1:]
            if judge(sample_geno):
                fout.write(line+"\n")