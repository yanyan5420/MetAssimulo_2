import pandas as pd
import numpy as np
from itertools import groupby
from operator import itemgetter

from simulate_2D.peak_detection_2d import peak_cluster_detection
# from simulate_2D.preprocess_1d_spectra import smooth_spectra


def calculate_peak_shift(data, x_scale, temp_pka, temp_ph, filtered_signal_subset, delta_acid_base_list):
    ph_standard = 7.4
    temp_data = np.array(data)
    shift_2d_array = np.zeros(temp_data.shape)
    step_size = x_scale[0] - x_scale[1]

    for idx, subset in enumerate(filtered_signal_subset):
        delta_acid_base = delta_acid_base_list[idx]
        delta_shift = (delta_acid_base * (10 ** (ph_standard - temp_pka) - 10 ** (temp_ph - temp_pka))) / \
                      ((1 + 10 ** (ph_standard - temp_pka)) * (1 + 10 ** (temp_ph - temp_pka)))
        if abs(delta_shift) > 0.5:
            if delta_shift > 0:
                delta_shift = 0.5
            else:
                delta_shift = -0.5
        shift_size = round(delta_shift / step_size)
        shift_subset = [s - shift_size for s in subset]
        print("shift_size:", shift_size, "delta_shift:", delta_shift, "step_size:", step_size,
              "temp_ph:", temp_ph, "temp_pka:", temp_pka, "delta_acid_base:", delta_acid_base)
        shift_2d_array[:, shift_subset] = temp_data[:, subset]
    return shift_2d_array


def construct_shift_data_for_all_repli(mixture_list, x_scale, norm_data_dict, meta_subset_dict,
                                       mixture_pka_dict, cons_ph_table_data, group_flag):

    ph_data_dict = dict()
    all_ph_dict = list(filter(lambda d: d["meta_name"] == 'pH', cons_ph_table_data))[0]
    print("all_ph_dict: ", all_ph_dict)
    print()
    group_ph_dict = dict(filter(lambda i: i[0].startswith(group_flag+"_replicate"), all_ph_dict.items()))
    print("group_ph_dict: ", group_ph_dict)
    print()
    for repli_name, repli_ph in group_ph_dict.items():
        shift_data_dict = dict()
        for meta_name in mixture_list:
            temp_pka = mixture_pka_dict[meta_name]
            temp_data = np.array(norm_data_dict[meta_name])
            filtered_signal_subset, delta_acid_base_list = meta_subset_dict[meta_name]
            temp_data_shift = calculate_peak_shift(temp_data, x_scale, temp_pka, float(repli_ph), filtered_signal_subset,
                                                delta_acid_base_list)
            shift_data_dict[meta_name] = temp_data_shift
        # smooth_data_dict = smooth_spectra(shift_data_dict, 0.05)
        ph_data_dict[repli_name] = shift_data_dict
    return group_ph_dict, ph_data_dict


def construct_shift_data_continuous_for_all_repli(mixture_list, x_scale, norm_data_dict, meta_subset_dict,
                                                  mixture_pka_dict, cons_ph_table_data):

    ph_data_dict = dict()
    all_ph_dict = list(filter(lambda d: d["meta_name"] == 'pH', cons_ph_table_data))[0]
    temp_ph_dict = all_ph_dict.copy()
    del temp_ph_dict['meta_name']
    del temp_ph_dict['hmdb_id']
    for repli_name, repli_ph in temp_ph_dict.items():
        shift_data_dict = dict()
        for meta_name in mixture_list:
            temp_pka = mixture_pka_dict[meta_name]
            temp_y = norm_data_dict[meta_name]
            filtered_signal_subset, delta_acid_base_list = meta_subset_dict[meta_name]
            temp_y_shift = calculate_peak_shift(temp_y, x_scale, temp_pka, float(repli_ph), filtered_signal_subset,
                                                delta_acid_base_list)
            shift_data_dict[meta_name] = temp_y_shift
        # smooth_data_dict = smooth_spectra(shift_data_dict, 0.05)
        ph_data_dict[repli_name] = shift_data_dict

    return temp_ph_dict, ph_data_dict



