import logging
from pysam import VariantFile
import pandas as pd
import argparse
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
        description=None, prog='sv_gt_evo.py',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # required arguments
    required_parser = parser.add_argument_group('required arguments')
    required_parser.add_argument('-v1', '--vcf1', action='store',
                                 type=str, dest='vcf1_path',required=True,
                                 help='real_vcf_path')
    required_parser.add_argument('-v2', '--vcf2', action='store',
                                 type=str, dest='vcf2_path',required=True,
                                 help='caller_vcf_path')
    required_parser.add_argument('-s', '--sample', action='store',
                                 type=str, dest='sample_name',required=True,
                                 help='sample name')
    required_parser.add_argument('-o', '--out', action='store', dest='out', type=str, required=True,
                                 help='output prefix')

    return parser.parse_args()


args = get_args()

REAL_VCF_PATH = args.vcf1_path
CALLER_VCF_PTAH = args.vcf2_path
SAMPLE_NAME = args.sample_name
OUTPUT_TABLE_FILE = f'{args.out}.tsv'
OUTPUT_COUNT_TABLE_FILE = f'{args.out}_count.tsv'
OUTPUT_SANKEY_JSON_FILE = f'{args.out}_sankey.json'

COMAPRE_TABLE = {
    "1/1@1/1": "TP", # True Positive
    "0/0@0/0": "TN", # True Negative
    "1/1@0/1": "HP", # Heterozygous Positive
    "0/0@0/1": "HN", # Heterozygous Negative
    "1/1@1/0": "HP", # Heterozygous Positive
    "0/0@1/0": "HN", # Heterozygous Negative
    "0/0@1/1": "FP", # False Positive
    "1/1@0/0": "FN", # False Negative
    "1/1@None/None": "NP", # None Positive
    "0/0@None/None": "NN", # None Negative
}

GT_SYNBLOIC = {
    "1/1": "T",
    "0/0": "F",
    "0/1": "H",
    "1/0": "H",
    "None/None": "N"
}

# 1. get SVLEN/TYPE/ID/GT info from REAL

REAL_VCF = VariantFile(REAL_VCF_PATH)

output_infos = {}
for record in track(REAL_VCF.fetch(), description="get SV info from REAL"):
    output_info = {}
    sv_type = record.info['SVTYPE']
    if sv_type == 'DEL':
        sv_len = len(record.ref)
    else:
        sv_len = len(record.alts[0])
    sv_id = record.id
    sv_gt = record.samples[SAMPLE_NAME]['GT']
    sv_gt_fmt = "/".join([str(i) for i in sv_gt])
    sv_real_gt_syn = GT_SYNBLOIC[sv_gt_fmt]

    output_infos[sv_id] = {
                'sv_type': sv_type,
                'sv_len': sv_len,
                'sv_id': sv_id,
                'sv_real_gt': sv_gt_fmt,
                'sv_real_gt_syn': sv_real_gt_syn,
            }


# 2. get GT info from CALLER and compare
CALLER_VCF = VariantFile(CALLER_VCF_PTAH)

for record in track(CALLER_VCF.fetch(), description="get SV info from CALLER"):
    sv_id = record.id
    sv_gt = record.samples[SAMPLE_NAME]['GT']
    sv_gt_fmt = "/".join([str(i) for i in sv_gt])
    sv_call_gt_syn = GT_SYNBLOIC[sv_gt_fmt]
    output_infos[sv_id]['sv_caller_gt'] = sv_gt_fmt
    output_infos[sv_id]['sv_caller_gt_syn'] = sv_call_gt_syn
    gt_comapre = output_infos[sv_id]['sv_real_gt'] + '@' + sv_gt_fmt
    gt_judge = COMAPRE_TABLE[gt_comapre]
    output_infos[sv_id]["gt_judge"] = gt_judge

# 3. statistics and output

# stat_df = pd.DataFrame(output_infos.values())




with open(OUTPUT_TABLE_FILE, 'w') as f:
    f.write('sv_type\tsv_len\tsv_id\tsv_real_gt_syn\tsv_caller_gt_syn\tsv_real_gt\tsv_caller_gt\tgt_judge\n')
    for item in track(output_infos.values(), description="1:output table"):
        f.write(f'{item["sv_type"]}\t{item["sv_len"]}\t{item["sv_id"]}\t{item["sv_real_gt_syn"]}\t{item["sv_caller_gt_syn"]}\t{item["sv_real_gt"]}\t{item["sv_caller_gt"]}\t{item["gt_judge"]}\n')


tmp_dict = {
    "lv1": [],
    "lv2": [],
    "lv3": [],
    "lv4": [],
}
for item in track(output_infos.values(), description="2:output count table"):

    lv1 = 'REAL'
    sv_real_syn = item['sv_real_gt_syn']
    sv_type = item['sv_type']
    sv_call_syn = item['sv_caller_gt_syn']

    lv2 = f'{lv1}-{sv_real_syn}'
    lv3 = f'{lv2}-{sv_type}'
    lv4 = f'{lv3}-{sv_call_syn}-CALL'

    tmp_dict['lv1'].append(lv1)
    tmp_dict['lv2'].append(lv2)
    tmp_dict['lv3'].append(lv3)
    tmp_dict['lv4'].append(lv4)


stat_df=pd.DataFrame(tmp_dict)
qq=stat_df.groupby(['lv1','lv2','lv3','lv4']).size().reset_index(name='counts')
qq.to_csv(OUTPUT_COUNT_TABLE_FILE, index=False, sep='\t')

def genSankey(df,cat_cols=[],value_cols='',title='Sankey Diagram'):

    nodes = []
    for catCol in cat_cols:
        labelListTemp =  list(set(df[catCol].values))
        nodes = nodes + labelListTemp
    # remove duplicates from labelList
    nodes = list(dict.fromkeys(nodes))

    nodes_with_color = []
    for item in nodes:
        if item == 'REAL':
            nodes_with_color.append({'name':item,'itemStyle':{'color':'#c7522a','borderColor':'#c7522a'}})
        elif item == 'REAL-F':
            nodes_with_color.append({'name':item,'itemStyle':{'color':'#faaac7','borderColor':'#faaac7'}})
        elif item == 'REAL-T':
            nodes_with_color.append({'name':item,'itemStyle':{'color':'#dde5b4','borderColor':'#dde5b4'}})
        elif 'INS' in item:
            nodes_with_color.append({'name':item,'itemStyle':{'color':'#43b0f1','borderColor':'#43b0f1'}})
        elif 'DEL' in item:
            nodes_with_color.append({'name':item,'itemStyle':{'color':'#008585','borderColor':'#008585'}})
        else:
            pass


    # transform df into a source-target pair
    for i in range(len(cat_cols)-1):
        if i==0:
            sourceTargetDf = df[[cat_cols[i],cat_cols[i+1],value_cols]]
            sourceTargetDf.columns = ['source','target','count']
        else:
            tempDf = df[[cat_cols[i],cat_cols[i+1],value_cols]]
            tempDf.columns = ['source','target','count']
            sourceTargetDf = pd.concat([sourceTargetDf,tempDf])
        sourceTargetDf = sourceTargetDf.groupby(['source','target']).agg({'count':'sum'}).reset_index()
    # creating the sankey diagram
    raw_links = dict(
          source = sourceTargetDf['source'].to_list(),
          target = sourceTargetDf['target'].to_list(),
          value = sourceTargetDf['count'].to_list()
        )
    links = []
    for i in range(len(raw_links['source'])):
        links.append(dict(source=raw_links['source'][i],
                          target=raw_links['target'][i],
                          value=raw_links['value'][i]))
    return nodes_with_color, links


with console.status("[bold green]3. Generating Sankey Diagram Config file...") as s:
    nodes, links = genSankey(qq,cat_cols=['lv1','lv2','lv3','lv4'],value_cols='counts',title='Sankey Diagram')
    console.log(f"3. output sankey config file: {OUTPUT_SANKEY_JSON_FILE}")
sankey_json = {
    "option":
        {
            "tooltip": {"trigger": "item", "triggerOn": "mousemove"},
            "series": {
                "type": "sankey",
                "emphasis": {"focus": "adjacency"},
                "data": nodes,
                "links": links,
                "levels": [
                    {
                        "depth": 0,
                        "itemStyle": {
                            "color": '#fbb4ae'
                        },
                        "lineStyle": {
                            "color": 'source',
                            "opacity": 0.6
                        }
                    },
                    {
                        "depth": 1,
                        "itemStyle": {
                            "color": '#b3cde3'
                        },
                        "lineStyle": {
                            "color": 'source',
                            "opacity": 0.6
                        }
                    },
                    {
                        "depth": 2,
                        "itemStyle": {
                            "color": '#ccebc5'
                        },
                        "lineStyle": {
                            "color": 'source',
                            "opacity": 0.6
                        }
                    },
                    {
                        "depth": 3,
                        "itemStyle": {
                            "color": '#decbe4'
                        },
                        "lineStyle": {
                            "color": 'source',
                            "opacity": 0.6
                        }
                    }
                    ],
                "lineStyle": {"curveness": 0.5}
            }

        }
}

import json

with open(OUTPUT_SANKEY_JSON_FILE, 'w') as f:
    json.dump(sankey_json, f, indent=4)