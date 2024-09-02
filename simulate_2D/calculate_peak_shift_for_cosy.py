import nmrglue as ng
import pandas as pd
import numpy as np
from pathlib import Path
import os
import glob
import pathlib
import plotly.express as px
import plotly.graph_objects as go
from scipy.signal import convolve2d
from scipy.stats import truncnorm
from scipy.signal import find_peaks
from itertools import groupby
from operator import itemgetter
import matplotlib.pyplot as plt
import copy


def get_projection_f1(matrix):
    p_f1 = []
    for i in range(matrix.shape[1]):
        p = np.max(matrix[i, :])
        p_f1.append(p)
    return p_f1


def peak_cluster_detection(y, find_peak_thres):
    # find multiplets width
    signals = np.zeros(len(y))
    mean = np.mean(y)
    for idx, value in enumerate(y):
        if value > mean:
            signals[idx] = 1
        else:
            signals[idx] = 0

    signal_index_list = np.where(signals == 1)[0]
    signal_width_subset = []
    for k, g in groupby(enumerate(signal_index_list), lambda ix: ix[0] - ix[1]):
        signal_width_subset.append(list(map(itemgetter(1), g)))

    # check signal width subsets
    # find peaks
    peaks_index_list, _ = find_peaks(y, height=max(y) * find_peak_thres)
    filtered_signal_subset = []
    for subset in signal_width_subset:
        if bool(set(peaks_index_list) & set(subset)) and len(subset) > 3:
            filtered_signal_subset.append(subset)

    # sample delta_acid_base for each peak cluster
    delta_acid_base_list = np.random.normal(-0.118, 0.204, len(filtered_signal_subset))

    return peaks_index_list, filtered_signal_subset, delta_acid_base_list


def calculate_peak_shift(x_scale, temp_pka, temp_ph, filtered_signal_subset, delta_acid_base_list):
    ph_standard = 7.4
    step_size = x_scale[0] - x_scale[1]
    shift_subset_list = []

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
        # shift_size = -100
        shift_subset = [s - shift_size for s in subset]
        print(subset, shift_subset)
        print("shift_size:", shift_size, "delta_shift:", delta_shift, "step_size:", step_size,
              "temp_ph:", temp_ph, "temp_pka:", temp_pka, "delta_acid_base:", delta_acid_base)
        print()

        shift_subset_list.append([shift_size, subset, shift_subset])

    return shift_subset_list


def modified_shift_on_f2(shifted_temp_data, shift_subset_list):
    # shifted_temp_data = data.copy()

    empty_temp_data = np.zeros(shifted_temp_data.shape)

    for shift_size, subset, shift_subset in shift_subset_list:
        # shift on the F2 axis
        empty_temp_data[:, shift_subset] = shifted_temp_data[:, subset]

    return empty_temp_data


def modified_shift_on_f1(shifted_temp_data, shift_subset_list):
    empty_data = np.zeros(shifted_temp_data.shape)

    shifted_temp_data = np.transpose(shifted_temp_data)[::-1, :]

    for shift_size, subset, shift_subset in shift_subset_list:
        empty_data[:, shift_subset] = shifted_temp_data[:, subset]

    empty_data = np.transpose(empty_data)[:, ::-1]
    return empty_data


def get_shifted_data_for_pure_compounds(meta, norm_data_dict, x_scale, temp_pka, temp_ph):
    temp_data = np.array(norm_data_dict[meta])
    temp_p_f1 = get_projection_f1(temp_data)

    peaks_index_list, filtered_signal_subset, delta_acid_base_list = peak_cluster_detection(temp_p_f1, 0.05)
    shifted_subset = calculate_peak_shift(x_scale, temp_pka, temp_ph, filtered_signal_subset, delta_acid_base_list)

    temp_data_copy = copy.deepcopy(temp_data)
    shifted_data_on_f2 = modified_shift_on_f2(temp_data_copy, shifted_subset)
    final_shifted_data = modified_shift_on_f1(shifted_data_on_f2, shifted_subset)

    return final_shifted_data


def get_shifted_data_for_each_replicate(mixture_list, norm_data_dict, x_scale, mixture_pka_dict, temp_ph):
    shifted_data_dict = dict()
    for meta_name in mixture_list:
        temp_pka = mixture_pka_dict[meta_name]
        shift_data = get_shifted_data_for_pure_compounds(meta_name, norm_data_dict, x_scale, temp_pka, temp_ph)
        shifted_data_dict[meta_name] = np.array(shift_data)
    return shifted_data_dict


def sum_mixture_data_for_each_replicate(repli_name, mixture_dict, shift_data_dict, cons_table_rows, protons_df, snr):
    sum_data = 0
    for meta_name, hmdb_id in mixture_dict.items():
        temp_dict = list(filter(lambda t: t['meta_name'] == meta_name, cons_table_rows))[0]
        temp_cons = float(temp_dict[repli_name])
        temp_protons = int(protons_df.loc[hmdb_id, "number_of_protons"])
        temp_data = np.array(shift_data_dict[meta_name])
        temp_intensity = temp_data * temp_cons * temp_protons
        sum_data = sum_data + temp_intensity

    # max_level = np.max(np.array(sum_data))
    # noise = np.random.normal(0, max_level / snr, sum_data.shape)
    #
    # sum_data = sum_data + noise
    return sum_data


# ------------------------------ shift data for group -------------------------------
def get_shifted_data_for_all_replicates(group_flag, group_repli_ph_dict, mixture_list, norm_data_dict, x_scale,
                                        mixture_pka_dict):

    group_ph_dict = dict(filter(lambda i: i[0].startswith(group_flag + "_replicate"), group_repli_ph_dict.items()))
    group_shifted_data_dict = dict()
    for repli_name, repli_ph in group_ph_dict.items():
        temp_ph = float(repli_ph)
        shift_data_dict = get_shifted_data_for_each_replicate(mixture_list, norm_data_dict, x_scale, mixture_pka_dict, temp_ph)
        group_shifted_data_dict[repli_name] = shift_data_dict
    return group_ph_dict, group_shifted_data_dict


def get_mixture_data_for_all_replicates(group_flag, group_repli_ph_dict, mixture_dict, norm_data_dict, mixture_pka_dict,
                                        x_scale, cons_ph_table_data, protons_df, snr):
    mixture_list = list(mixture_dict.keys())
    group_ph_dict, group_shifted_data_dict = get_shifted_data_for_all_replicates(group_flag, group_repli_ph_dict,
                                                                                 mixture_list, norm_data_dict, x_scale,
                                                                                 mixture_pka_dict)

    repli_mix_data_dict = dict()
    for repli_name, shift_data_dict in group_shifted_data_dict.items():
        temp_sum_data = sum_mixture_data_for_each_replicate(repli_name, mixture_dict, shift_data_dict, cons_ph_table_data,
                                                            protons_df, snr)
        repli_ph = group_ph_dict[repli_name]
        repli_mix_data_dict[repli_name.lstrip(group_flag+"_")] = [repli_ph, temp_sum_data]

    return repli_mix_data_dict


# ------------------------------ shift data for continuous -------------------------------
def conti_get_shifted_data_for_all_replicates(conti_repli_ph_dict, mixture_list, norm_data_dict, x_scale,
                                              mixture_pka_dict):

    conti_ph_dict = copy.deepcopy(conti_repli_ph_dict)
    del conti_ph_dict['meta_name']
    del conti_ph_dict['hmdb_id']

    conti_shifted_data_dict = dict()
    for repli_name, repli_ph in conti_ph_dict.items():
        temp_ph = float(repli_ph)
        shift_data_dict = get_shifted_data_for_each_replicate(mixture_list, norm_data_dict, x_scale, mixture_pka_dict,
                                                              temp_ph)
        conti_shifted_data_dict[repli_name] = shift_data_dict
    return conti_ph_dict, conti_shifted_data_dict


def conti_get_mixture_data_for_all_replicates(conti_repli_ph_dict, mixture_dict, norm_data_dict, mixture_pka_dict,
                                              x_scale, cons_ph_table_data, protons_df, snr):
    mixture_list = list(mixture_dict.keys())
    conti_ph_dict, conti_shifted_data_dict = conti_get_shifted_data_for_all_replicates(conti_repli_ph_dict,
                                                                                       mixture_list, norm_data_dict,
                                                                                       x_scale, mixture_pka_dict)

    repli_mix_data_dict = dict()
    for repli_name, shift_data_dict in conti_shifted_data_dict.items():
        temp_sum_data = sum_mixture_data_for_each_replicate(repli_name, mixture_dict, shift_data_dict,
                                                            cons_ph_table_data,
                                                            protons_df, snr)
        repli_ph = conti_ph_dict[repli_name]
        repli_mix_data_dict[repli_name] = [repli_ph, temp_sum_data]

    return repli_mix_data_dict

