import argparse
import sys
import pysam
from tqdm import tqdm
from rich.progress import track

def get_args():
    # define arguments
    parser = argparse.ArgumentParser(
        description=None, prog='bam_chr.py',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # required arguments
    required_parser = parser.add_argument_group('required arguments')
    required_parser.add_argument('-m', '--mode', action='store',
                                 type=str, dest='mode',required=True,
                                 help='bam or paf')
    required_parser.add_argument('-i', '--input', action='store', dest='input', type=str, required=True,
                                 help='input file')
    required_parser.add_argument('-o', '--output', action='store', dest='output', type=str, required=True,
                                 help='output file')


    # options arguments
    option_parser = parser.add_argument_group('optional arguments')
    option_parser.add_argument('-s','--size', action='store', dest='indelsize', type=int, default=5,
                               help='allowd indel size to fix the hole')
    option_parser.add_argument('-l','--len', action='store', dest='file_lines', type=int,
                               help='allowd indel size to fix the hole')
    return parser.parse_args()




def filter_bam(bamfile_path: str, outbamfile_path: str, size: int=10000000):
    bamfile = pysam.AlignmentFile(bamfile_path, 'rb')

    chro_exact_reads = pysam.AlignmentFile(outbamfile_path, "wb", template=bamfile)

    idx_stats = bamfile.get_index_statistics()
    contig_reads = sum([x.total for x in idx_stats])

    for read in track(bamfile.fetch(), total=contig_reads):
        ref_chro = f"chr{read.reference_name}_"
        query_name = read.query_name
        if query_name.startswith(ref_chro):
            query_start = int(query_name.split("_")[1])
            query_end = int(query_name.split("_")[2])
            ref_start = read.reference_start
            query_range  = range(query_start-size, query_end+size)
            if ref_start in query_range:
                chro_exact_reads.write(read)
            # chro_exact_reads.write(read)
        else:
            continue

    chro_exact_reads.close()
    bamfile.close()

def filter_paf(paffile_path: str, outpaffile_path: str, file_lines: int, size: int=10000000):
    paffile = open(paffile_path, 'r')
    outpaffile = open(outpaffile_path, 'w')
    file_lines = file_lines if file_lines else sum(1 for line in open(paffile_path))
    for line in track(paffile, total=file_lines):
        line = line.strip()
        ref_chro = line.split("\t")[5]
        query_name = line.split()[0]
        query_chro = query_name.split('_')[0].strip('chr')
        if query_chro == ref_chro:
            query_start = int(query_name.split("_")[1])
            query_end = int(query_name.split("_")[2])
            ref_start = int(line.split("\t")[7])
            query_range  = range(query_start-size, query_end+size)
            if ref_start in query_range:
                outpaffile.write(line+'\n')
            # chro_exact_reads.write(read)
        else:
            continue

    outpaffile.close()
    paffile.close()


if __name__ == '__main__':
    args = get_args()
    if args.mode == 'bam':
        filter_bam(args.input, args.output, args.indelsize)
    if args.mode == 'paf':
        filter_paf(args.input, args.output, args.file_lines, args.indelsize)
    else:
        print('FUCK! choose mode == bam or paf')
