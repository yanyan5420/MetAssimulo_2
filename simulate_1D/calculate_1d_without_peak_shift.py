import numpy as np
from scipy import signal


def sum_mixture_spectra_with_albumin_for_each_repli(n, mixture_dict, format_norm_data_dict, cons_table_rows, protons_df,
                                     albumin_norm_data_dict_1, albumin_level, snr, wins, group_flag):

    sum_data = 0
    for meta_name, hmdb_id in mixture_dict.items():
        temp_dict = list(filter(lambda t: t['meta_name'] == meta_name, cons_table_rows))[0]
        temp_cons = float(temp_dict[group_flag + "_replicate_" + str(n + 1)])
        try:
            temp_protons = int(protons_df.loc[hmdb_id, "number_of_protons"])
        except:
            temp_protons = 1
        temp_data = np.array(format_norm_data_dict[meta_name])
        temp_intensity = temp_data * temp_cons * temp_protons
        sum_data = sum_data + temp_intensity

    sum_data = sum_data + albumin_norm_data_dict_1["Albumin"] * 678 * albumin_level

    return sum_data


def sum_mixture_spectra_with_albumin_for_all_repli(num_replicates, mixture_dict, format_norm_data_dict,
                                   cons_table_rows, protons_df, albumin_norm_data_dict_1, albumin_level, snr, wins, group_flag):
    replicate_dict = dict()
    for n in range(num_replicates):
        temp_sum_data = sum_mixture_spectra_with_albumin_for_each_repli(n, mixture_dict, format_norm_data_dict,
                                                                        cons_table_rows, protons_df,
                                                                        albumin_norm_data_dict_1, albumin_level,
                                                                        snr, wins, group_flag)
        replicate_dict['replicate_' + str(n + 1)] = temp_sum_data

    mean_repli_data = sum(replicate_dict.values()) / num_replicates
    replicate_dict["replicate_mean"] = mean_repli_data

    return replicate_dict


def sum_mixture_for_each_repli(n, mixture_list, format_norm_data_dict, cons_table_rows, protons_df, snr, wins, group_flag):
    sum_data = 0
    for meta_name in mixture_list:
        temp_dict = list(filter(lambda t: t['meta_name'] == meta_name, cons_table_rows))[0]
        temp_cons = float(temp_dict[group_flag+"_replicate_"+str(n+1)])
        temp_protons = 1
            # int(protons_df.loc[meta_name, "number"])
        temp_data = np.array(format_norm_data_dict[meta_name])

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


def simulate_mixture_for_all_repli(num_replicates, mixture_list, format_norm_data_dict,
                                   cons_table_rows, protons_df, snr, wins, group_flag):
    replicate_dict = dict()
    for n in range(num_replicates):
        temp_sum_data = sum_mixture_for_each_repli(n, mixture_list, format_norm_data_dict,
                                                   cons_table_rows, protons_df, snr, wins, group_flag)
        replicate_dict['replicate_' + str(n + 1)] = temp_sum_data

    mean_repli_data = sum(replicate_dict.values()) / num_replicates
    replicate_dict["replicate_mean"] = mean_repli_data

    return replicate_dict


def conti_sum_mixture_for_each_repli(n, mixture_dict, format_norm_data_dict, cons_table_rows, protons_df,
                                     albumin_norm_data_dict_1, albumin_level, snr, wins):
    sum_data = 0
    for meta_name, hmdb_id in mixture_dict.items():
        temp_dict = list(filter(lambda t: t['meta_name'] == meta_name, cons_table_rows))[0]
        temp_cons = float(temp_dict["replicate_"+str(n+1)])
        try:
            temp_protons = int(protons_df.loc[hmdb_id, "number_of_protons"])
        except:
            temp_protons = 1
        temp_data = np.array(format_norm_data_dict[meta_name])
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


def simulate_continuous_mixture_for_all_repli(num_replicates, mixture_dict, format_norm_data_dict,
                                              cons_table_rows, protons_df, albumin_norm_data_dict_1, albumin_level,
                                              snr, wins):
    replicate_dict = dict()
    for n in range(num_replicates):
        temp_sum_data = conti_sum_mixture_for_each_repli(n, mixture_dict, format_norm_data_dict,
                                                         cons_table_rows, protons_df,
                                                         albumin_norm_data_dict_1, albumin_level,
                                                         snr, wins)
        replicate_dict['replicate_' + str(n + 1)] = temp_sum_data

    mean_repli_data = sum(replicate_dict.values()) / num_replicates
    replicate_dict["replicate_mean"] = mean_repli_data

    return replicate_dict


