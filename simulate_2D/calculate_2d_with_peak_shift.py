import pandas as pd
import numpy as np

from simulate_2D.peak_shift_2d import construct_shift_data_for_all_repli, construct_shift_data_continuous_for_all_repli


def sum_mixture_for_each_repli(repli_name, mixture_list, shift_data_dict, cons_table_rows, protons_df, snr):
    sum_data = 0
    for meta_name in mixture_list:
        temp_dict = list(filter(lambda t: t['meta_name'] == meta_name, cons_table_rows))[0]
        temp_cons = float(temp_dict[repli_name])
        temp_protons = 1
            # int(protons_df.loc[meta_name, "number"])
        temp_data = np.array(shift_data_dict[meta_name])

        temp_intensity = temp_data * temp_cons * temp_protons
        sum_data = sum_data + temp_intensity

    # max_level = np.max(np.array(sum_data))
    # noise = np.random.normal(0, max_level / snr, sum_data.shape)
    #
    # sum_data = sum_data + noise

    return sum_data


def simulate_mixture_with_peak_shift_for_all_repli(mixture_list, x_scale, norm_data_dict, meta_subset_dict,
                                                   mixture_pka_dict, cons_ph_table_data,
                                                   group_flag, protons_df, snr):
    replicate_dict = dict()
    group_ph_dict, ph_data_dict = construct_shift_data_for_all_repli(mixture_list, x_scale, norm_data_dict, meta_subset_dict,
                                       mixture_pka_dict, cons_ph_table_data, group_flag)
    for repli_name, shift_data_dict in ph_data_dict.items():
        temp_sum_data = sum_mixture_for_each_repli(repli_name, mixture_list, shift_data_dict,
                                                   cons_ph_table_data, protons_df, snr)
        repli_ph = group_ph_dict[repli_name]
        replicate_dict[repli_name.lstrip(group_flag+"_")] = [repli_ph, temp_sum_data]
    # replicate_dict["replicate_mean"] = np.mean(list(map(lambda x: x[1], replicate_dict.values())), axis=0)
    return replicate_dict


def get_shift_p_jres_for_all_repli(replicate_dict):
    num_replicates = len(replicate_dict)

    shift_p_jres_dict = dict()
    for n in range(num_replicates):
        temp_repli = 'replicate_' + str(n + 1)
        str_ph, data = replicate_dict[temp_repli]
        temp_data = np.array(data)
        shift_p_jres_scale = np.zeros(temp_data.shape[1])
        for i in range(temp_data.shape[1]):
            temp_y = max(temp_data[:, i])
            shift_p_jres_scale[i] = temp_y
        shift_p_jres_dict[temp_repli] = [str_ph, shift_p_jres_scale]

    return shift_p_jres_dict


def simulate_mixture_continuous_with_peak_shift_for_all_repli(mixture_list, ppm_scale, norm_data_dict, meta_subset_dict,
                                                   mixture_pka_dict, cons_ph_table_data, protons_df, snr):
    replicate_dict = dict()
    conti_ph_dict, ph_data_dict = construct_shift_data_continuous_for_all_repli(mixture_list, ppm_scale, norm_data_dict,
                                                            meta_subset_dict, mixture_pka_dict, cons_ph_table_data)
    for repli_name, shift_data_dict in ph_data_dict.items():
        temp_sum_data = sum_mixture_for_each_repli(repli_name, mixture_list, shift_data_dict, cons_ph_table_data,
                                                   protons_df, snr)
        repli_ph = conti_ph_dict[repli_name]
        replicate_dict[repli_name] = [repli_ph, temp_sum_data]

    # replicate_dict["replicate_mean"] = np.mean(list(map(lambda x: x[1], replicate_dict.values())), axis=0)
    return replicate_dict
