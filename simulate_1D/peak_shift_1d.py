import pandas as pd
import numpy as np
from itertools import groupby
from operator import itemgetter
import copy

# from simulate_1D.peak_detection_1d import peak_cluster_detection
from simulate_1D.preprocess_1d_spectra import smooth_spectra


def calculate_peak_shift(y, ppm_scale, temp_pka, temp_ph, filtered_signal_subset, delta_acid_base_list):
    ph_standard = 7.4
    y = np.array(y)
    y_shift_list = np.zeros(len(y))
    step_size = ppm_scale[0] - ppm_scale[1]

    for idx, subset in enumerate(filtered_signal_subset):
        delta_acid_base = delta_acid_base_list[idx]
        delta_shift = (delta_acid_base * (10 ** (ph_standard - temp_pka) - 10 ** (temp_ph - temp_pka))) / \
                      ((1 + 10 ** (ph_standard - temp_pka)) * (1 + 10 ** (temp_ph - temp_pka)))
            # (delta_acid_base * (10 ** (ph_standard - temp_ph) - 10 ** (temp_ph - temp_pka))) / \
            #           ((1 + 10 ** (ph_standard - temp_pka)) * (1 + 10 ** (temp_ph - temp_pka)))
        if abs(delta_shift) > 0.5:
            if delta_shift > 0:
                delta_shift = 0.5
            else:
                delta_shift = -0.5
        shift_size = round(delta_shift / step_size)
        shift_subset = [s - shift_size for s in subset]
        print("shift_size:", shift_size, "delta_shift:", delta_shift, "step_size:", step_size,
              "temp_ph:", temp_ph, "temp_pka:", temp_pka, "delta_acid_base:", delta_acid_base)
        # print(len(shift_subset), shift_subset)
        # print(len(subset), subset)
        # print("shift_subset: ", shift_subset, "subset: ", subset)
        y_shift_list[shift_subset] = y[subset]

    return y_shift_list


def construct_shift_data_for_all_repli(mixture_list, ppm_scale, norm_data_dict, meta_subset_dict,
                                       mixture_pka_dict, cons_ph_table_data, group_flag):

    ph_data_dict = dict()
    all_ph_dict = list(filter(lambda d: d["meta_name"] == 'pH', cons_ph_table_data))[0]
    group_ph_dict = dict(filter(lambda i: i[0].startswith(group_flag+"_replicate"), all_ph_dict.items()))
    for repli_name, repli_ph in group_ph_dict.items():
        shift_data_dict = dict()
        for meta_name in mixture_list:
            temp_pka = mixture_pka_dict[meta_name]
            temp_y = norm_data_dict[meta_name]
            # print("temp_y: ", temp_y)
            filtered_signal_subset, delta_acid_base_list = meta_subset_dict[meta_name]
            print(meta_name)
            temp_y_shift = calculate_peak_shift(temp_y, ppm_scale, temp_pka, float(repli_ph), filtered_signal_subset,
                                                delta_acid_base_list)
            shift_data_dict[meta_name] = temp_y_shift
        smooth_data_dict = smooth_spectra(shift_data_dict, 0.05)
        ph_data_dict[repli_name] = smooth_data_dict

    return group_ph_dict, ph_data_dict


def construct_shift_data_continuous_for_all_repli(mixture_list, ppm_scale, norm_data_dict, meta_subset_dict,
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
            temp_y_shift = calculate_peak_shift(temp_y, ppm_scale, temp_pka, float(repli_ph), filtered_signal_subset,
                                                delta_acid_base_list)
            shift_data_dict[meta_name] = temp_y_shift
        smooth_data_dict = smooth_spectra(shift_data_dict, 0.05)
        ph_data_dict[repli_name] = smooth_data_dict

    return temp_ph_dict, ph_data_dict
