ref_path_list = [1,3,4,5,6]
sv_node = [
    {"order":1, "node_nums":[2], "node_type": "INS"},
    {"order":2, "node_nums":[4], "node_type": "DEL"},
]
sample_gt = [
    {"name": "Sample1", "gt": [1,0]},
    {"name": "Sample2", "gt": [0,1]},
    {"name": "Sample3", "gt": [1,1]},
]
def output(ref_path_list, sv_node, sample_gt):
    for sample,gt_list in sample_gt.items():
        _ref_path_list = ref_path_list.copy()
        for idx, gt in enumerate(gt_list, start=0):
            operate_lst = sv_node[idx]["node_nums"]
            operate_sig = sv_node[idx]["node_type"]
            if gt == 1:
                if operate_sig == "INS":
                    _ref_path_list.extend(operate_lst)
                    # print(_ref_path_list)
                elif operate_sig == "DEL":
                    for del_ in operate_lst:
                        _ref_path_list.remove(del_)
                        # print(_ref_path_list)
                else:
                    raise ValueError("operate_sig error")
            else:
                pass
        _ref_path_list.sort()
        print(f"{sample}:{_ref_path_list}")

ref_path_list = [1,2,3,6,7,8,9,10,13,14]
sv_node = [
    {"order":1, "node_nums":[4,5], "node_type": "INS"},
    {"order":2, "node_nums":[8], "node_type": "DEL"},
    {"order":3, "node_nums":[11,12], "node_type": "INS"},
]
sample_gt = {
    "Sample1": [1,0,1],
    "Sample2": [1,1,0],
    "Sample3": [0,1,1],
}

# output(ref_path_list, sv_node, sample_gt)


def get_ref_path_list(raw_P_file):
    with open(raw_P_file) as f:
        first_line = f.readline()
        _ref_path_list = first_line.strip().split("\t")[2].split(',')
        ref_path_list = [int(i.strip('+')) for i in _ref_path_list]
        return ref_path_list


def get_sample_sv(raw_P_file, sv_gt_file):
    all_node_list = []
    with open(raw_P_file) as f:
        for line in f.readlines():
            if '_alt' in line:
                node_list = line.strip().split('\t')[2].split(',')
                node_list = [int(i.strip('+')) for i in node_list]
                all_node_list.append(node_list)
    ## sort all_node_list by firt_unit_node
    sorted_all_node_list = sorted(all_node_list, key=lambda x: x[0])

    sv_list = []
    sample_gt = {}
    with open(sv_gt_file) as f:
        for line in f.readlines():
            sv_type = line.strip().split('\t')[0].strip('>').strip('<')
            sv_list.append(sv_type)
            for sample in line.strip().split('\t')[1:]:
                sample_name = sample.split('=')[0]
                gt = sample.split('=')[1].split('/')[0]
                if sample_name not in sample_gt:
                    sample_gt[sample_name] = []
                    sample_gt[sample_name].append(int(gt))
                else:
                    sample_gt[sample_name].append(int(gt))

    sv_node = []
    for idx,sv in enumerate(sv_list, start=1):
        sv_node.append(
            {
                "order":idx,
                "node_nums":sorted_all_node_list[idx-1],
                "node_type": sv
                }
            )

    return sv_node, sample_gt


raw_P_file = 'raw.P'
ref_path_list = get_ref_path_list(raw_P_file)
print(ref_path_list)
sv_gt_file = 'sample.gt'
sv_node, sample_gt = get_sample_sv(raw_P_file, sv_gt_file)
print(sv_node)
print(sample_gt)
output(ref_path_list, sv_node, sample_gt)
