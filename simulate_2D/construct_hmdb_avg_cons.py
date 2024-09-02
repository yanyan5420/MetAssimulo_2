import pandas as pd
import numpy as np
import re


def estimate_mean_std(min_v, max_v):
    matrix_1 = np.mat("1, 3; 1, -3")
    matrix_2 = np.mat(max_v + "," + min_v).T
    mean = round(float(np.linalg.solve(matrix_1, matrix_2)[0]), 2)
    std = round(float(np.linalg.solve(matrix_1, matrix_2)[1]), 2)
    return mean, std


def get_hmdb_normal_avg_cons(hmdb_cons_dict, hmdb_id, bio_type):
    mean_list = []
    std_list = []
    unit_str = None

    for cons_dict in hmdb_cons_dict[hmdb_id][bio_type]:
        value = cons_dict['cons_value']
        unit = cons_dict['cons_unit']
        unit_str = unit

        pattern_1 = re.match(r"(?P<mean>\d+\.?\d+) ?\+\/\- ?(?P<std>\d+\.?\d+)", value)  # 88.0 +/- 33.0
        pattern_2 = re.match(r"(?P<min_value>\d+\.?\d+) ?\- ?(?P<max_value>\d+\.?\d+)", value)  # "30.00-400.0"
        # "190.0 (30.0-400.0)"
        pattern_3 = re.match(r"(?P<real_mean>\d+\.?\d+) ?\((?P<min_value>\d+\.?\d+) ?\- ?(?P<max_value>\d+\.?\d+)\)", value)
        pattern_4 = re.match(r"(?P<mean>\d+\.?\d+) ?\((?P<std>\d+\.?\d+)\)", value)  # "122.3(27.85)"

        if pattern_1:
            mean = float(pattern_1.group("mean"))
            std = float(pattern_1.group("std"))
            mean_list.append(mean)
            std_list.append(std)

        elif pattern_2:
            min_v = pattern_2.group("min_value")
            max_v = pattern_2.group("max_value")
            mean, std = estimate_mean_std(min_v, max_v)
            mean_list.append(mean)
            std_list.append(std)

        elif pattern_3:
            real_mean = float(pattern_3.group("real_mean"))
            min_v = pattern_3.group("min_value")
            max_v = pattern_3.group("max_value")
            mean, std = estimate_mean_std(min_v, max_v)
            mean_list.append(real_mean)
            std_list.append(std)

        elif pattern_4:
            mean = float(pattern_4.group("mean"))
            std = float(pattern_4.group("std"))
            mean_list.append(mean)
            std_list.append(std)

    avg_mean = round(float(np.mean(mean_list)), 4)
    avg_std = round(float(np.mean(std_list)), 4)

    if unit_str == "umol/mmol creatinine":

        creatinine_bridge = 11.99525
        avg_mean_um, avg_std_um = avg_mean * creatinine_bridge, avg_std * creatinine_bridge
        return round(avg_mean_um, 2), round(avg_std_um, 2), unit_str
        # hmdb_norm_csf_cons_dict[hmdb_id] = [avg_mean_um, avg_std_um, "uM"]

    else:
        avg_mean = round(float(np.mean(mean_list)), 2)
        avg_std = round(float(np.mean(std_list)), 2)
        return avg_mean, avg_std, unit_str


def get_hmdb_abnormal_avg_cons(hmdb_cons_dict, hmdb_id, bio_type):

    if bio_type == "Blood":

        mean_list = []
        std_list = []
        unit_str = None

        for cons_dict in hmdb_cons_dict[hmdb_id][bio_type]:
            if cons_dict["condition"] == "Heart Transplant":
                value = cons_dict['cons_value']
                unit = cons_dict['cons_unit']
                unit_str = unit

                pattern_1 = re.match(r"(?P<mean>\d+\.?\d+) ?\+\/\- ?(?P<std>\d+\.?\d+)", value)  # 88.0 +/- 33.0
                pattern_2 = re.match(r"(?P<min_value>\d+\.?\d+) ?\- ?(?P<max_value>\d+\.?\d+)", value)  # "30.00-400.0"
                # "190.0 (30.0-400.0)"
                pattern_3 = re.match(r"(?P<real_mean>\d+\.?\d+) ?\((?P<min_value>\d+\.?\d+) ?\- ?(?P<max_value>\d+\.?\d+)\)", value)
                pattern_4 = re.match(r"(?P<mean>\d+\.?\d+) ?\((?P<std>\d+\.?\d+)\)", value)  # "122.3(27.85)"

                if pattern_1:
                    mean = float(pattern_1.group("mean"))
                    std = float(pattern_1.group("std"))
                    mean_list.append(mean)
                    std_list.append(std)

                elif pattern_2:
                    min_v = pattern_2.group("min_value")
                    max_v = pattern_2.group("max_value")
                    mean, std = estimate_mean_std(min_v, max_v)
                    mean_list.append(mean)
                    std_list.append(std)

                elif pattern_3:
                    real_mean = float(pattern_3.group("real_mean"))
                    min_v = pattern_3.group("min_value")
                    max_v = pattern_3.group("max_value")
                    mean, std = estimate_mean_std(min_v, max_v)
                    mean_list.append(real_mean)
                    std_list.append(std)

                elif pattern_4:
                    mean = float(pattern_4.group("mean"))
                    std = float(pattern_4.group("std"))
                    mean_list.append(mean)
                    std_list.append(std)

        avg_mean = round(float(np.mean(mean_list)), 2)
        avg_std = round(float(np.mean(std_list)), 2)
        return avg_mean, avg_std, unit_str

    else:
        return 0, 0, None
