import pandas as pd
import numpy as np
from scipy import signal

from simulate_1D.peak_shift_1d import construct_shift_data_for_all_repli, construct_shift_data_continuous_for_all_repli


def sum_mixture_with_peak_shift_with_albumin_for_each_repli(repli_name, mixture_dict, shift_data_dict,
                                                            cons_table_rows, protons_df, albumin_norm_data_dict_1,
                                                            albumin_level, snr, wins):
    sum_data = 0
    for meta_name, hmdb_id in mixture_dict.items():
        temp_dict = list(filter(lambda t: t['meta_name'] == meta_name, cons_table_rows))[0]
        temp_cons = float(temp_dict[repli_name])
        # temp_protons = int(protons_df.loc[hmdb_id, "number_of_protons"])
        try:
            temp_protons = int(protons_df.loc[hmdb_id, "number_of_protons"])
        except:
            temp_protons = 1
        temp_data = np.array(shift_data_dict[meta_name])
        temp_intensity = temp_data * temp_cons * temp_protons
        sum_data = sum_data + temp_intensity

    sum_data = sum_data + albumin_norm_data_dict_1["Albumin"] * 678 * albumin_level
    return sum_data


def simulate_mixture_with_peak_shift_with_albumin_for_all_repli(mixture_dict, mixture_list, ppm_scale, norm_data_dict,
                                                                meta_subset_dict, mixture_pka_dict, cons_ph_table_data,
                                                                group_flag, protons_df, albumin_norm_data_dict_1,
                                                                albumin_level, snr, wins):
    replicate_dict = dict()
    group_ph_dict, ph_data_dict = construct_shift_data_for_all_repli(mixture_list, ppm_scale, norm_data_dict,
                                                                     meta_subset_dict, mixture_pka_dict,
                                                                     cons_ph_table_data, group_flag)
    for repli_name, shift_data_dict in ph_data_dict.items():
        temp_sum_data = sum_mixture_with_peak_shift_with_albumin_for_each_repli(repli_name, mixture_dict, shift_data_dict,
                                                            cons_ph_table_data, protons_df, albumin_norm_data_dict_1,
                                                            albumin_level, snr, wins)
        repli_ph = group_ph_dict[repli_name]
        replicate_dict[repli_name.lstrip(group_flag+"_")] = [repli_ph, temp_sum_data]

    replicate_dict["replicate_mean"] = np.mean(list(map(lambda x: x[1], replicate_dict.values())), axis=0)
    return replicate_dict



def sum_mixture_for_each_repli(repli_name, mixture_list, shift_data_dict, cons_table_rows, protons_df, snr, wins):
    sum_data = 0
    for meta_name in mixture_list:
        temp_dict = list(filter(lambda t: t['meta_name'] == meta_name, cons_table_rows))[0]
        temp_cons = float(temp_dict[repli_name])
        temp_protons = 1
            # int(protons_df.loc[meta_name, "number"])
        temp_data = np.array(shift_data_dict[meta_name])

        temp_intensity = temp_data * temp_cons * temp_protons
        sum_data = sum_data + temp_intensity

    # add noise, SNR = 1000
    noise_std = max(sum_data) / snr
    print(noise_std, len(sum_data))
    noise_y = np.random.normal(0, noise_std, len(sum_data))

    win = signal.windows.hann(wins)
    smooth_noise_y = signal.convolve(noise_y, win, mode='same') / sum(win)

    final_sum_data = sum_data + smooth_noise_y

    return final_sum_data


def simulate_mixture_with_peak_shift_for_all_repli(mixture_list, ppm_scale, norm_data_dict, meta_subset_dict,
                                                   mixture_pka_dict, cons_ph_table_data,
                                                   group_flag, protons_df, snr, wins):
    replicate_dict = dict()
    group_ph_dict, ph_data_dict = construct_shift_data_for_all_repli(mixture_list, ppm_scale, norm_data_dict,
                                                                     meta_subset_dict, mixture_pka_dict,
                                                                     cons_ph_table_data, group_flag)
    for repli_name, shift_data_dict in ph_data_dict.items():
        temp_sum_data = sum_mixture_for_each_repli(repli_name, mixture_list, shift_data_dict, cons_ph_table_data,
                                                   protons_df, snr, wins)
        repli_ph = group_ph_dict[repli_name]
        replicate_dict[repli_name.lstrip(group_flag+"_")] = [repli_ph, temp_sum_data]

    # replicate_dict["replicate_mean"] = np.mean(list(map(lambda x: x[1], replicate_dict.values())), axis=0)
    return replicate_dict


def sum_mixture_continuous_with_peak_shift_for_each_repli(repli_name, mixture_dict, shift_data_dict, cons_table_rows,
                                                          protons_df, albumin_norm_data_dict_1, albumin_level, snr, wins):
    sum_data = 0
    for meta_name, hmdb_id in mixture_dict.items():
        temp_dict = list(filter(lambda t: t['meta_name'] == meta_name, cons_table_rows))[0]
        temp_cons = float(temp_dict[repli_name])
        try:
            temp_protons = int(protons_df.loc[hmdb_id, "number_of_protons"])
        except:
            temp_protons = 1
        temp_data = np.array(shift_data_dict[meta_name])
        temp_intensity = temp_data * temp_cons * temp_protons
        sum_data = sum_data + temp_intensity

    sum_data = sum_data + albumin_norm_data_dict_1["Albumin"] * 678 * albumin_level

    # add noise, SNR = 1000
    noise_std = max(sum_data) / snr
    print(noise_std, len(sum_data))
    noise_y = np.random.normal(0, noise_std, len(sum_data))

    win = signal.windows.hann(wins)
    smooth_noise_y = signal.convolve(noise_y, win, mode='same') / sum(win)

    final_sum_data = sum_data + smooth_noise_y
    return final_sum_data


def simulate_mixture_continuous_with_peak_shift_for_all_repli_with_albumin(mixture_dict, mixture_list, ppm_scale, norm_data_dict,
                                                              meta_subset_dict, mixture_pka_dict, cons_ph_table_data,
                                                              protons_df, albumin_norm_data_dict_1, albumin_level,
                                                              snr, wins):
    replicate_dict = dict()
    conti_ph_dict, ph_data_dict = construct_shift_data_continuous_for_all_repli(mixture_list, ppm_scale, norm_data_dict,
                                                            meta_subset_dict, mixture_pka_dict, cons_ph_table_data)
    for repli_name, shift_data_dict in ph_data_dict.items():
        temp_sum_data = sum_mixture_continuous_with_peak_shift_for_each_repli(repli_name, mixture_dict, shift_data_dict,
                                                                              cons_ph_table_data, protons_df,
                                                                              albumin_norm_data_dict_1, albumin_level,
                                                                              snr, wins)
        repli_ph = conti_ph_dict[repli_name]
        replicate_dict[repli_name] = [repli_ph, temp_sum_data]

    # replicate_dict["replicate_mean"] = np.mean(list(map(lambda x: x[1], replicate_dict.values())), axis=0)
    return replicate_dict



def simulate_mixture_continuous_with_peak_shift_for_all_repli(mixture_list, ppm_scale, norm_data_dict, meta_subset_dict,
                                                   mixture_pka_dict, cons_ph_table_data, protons_df, snr, wins):
    replicate_dict = dict()
    conti_ph_dict, ph_data_dict = construct_shift_data_continuous_for_all_repli(mixture_list, ppm_scale, norm_data_dict,
                                                            meta_subset_dict, mixture_pka_dict, cons_ph_table_data)
    for repli_name, shift_data_dict in ph_data_dict.items():
        temp_sum_data = sum_mixture_for_each_repli(repli_name, mixture_list, shift_data_dict, cons_ph_table_data,
                                                   protons_df, snr, wins)
        repli_ph = conti_ph_dict[repli_name]
        replicate_dict[repli_name] = [repli_ph, temp_sum_data]

    # replicate_dict["replicate_mean"] = np.mean(list(map(lambda x: x[1], replicate_dict.values())), axis=0)
    return replicate_dict
