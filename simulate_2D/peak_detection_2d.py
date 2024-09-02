import pandas as pd
import numpy as np
from scipy.signal import find_peaks
from itertools import groupby
from operator import itemgetter


def get_p_jres_dict(mixture_list, final_data_dict):
    p_jres_dict = dict()
    for meta_name in mixture_list:
        temp_data = np.array(final_data_dict[meta_name])
        p_jres_scale = np.zeros(temp_data.shape[1])
        for i in range(temp_data.shape[1]):
            temp_y = max(temp_data[:, i])
            p_jres_scale[i] = temp_y
        p_jres_dict[meta_name] = p_jres_scale
    return p_jres_dict


def peak_cluster_detection(y):
    """
    y: the projection of JRes (pJRes)
    """
    # find peaks first
    peaks_index_list, _ = find_peaks(y, height=max(y)*0.1)

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
    filtered_signal_subset = []
    for subset in signal_width_subset:
        if bool(set(peaks_index_list) & set(subset)):
            filtered_signal_subset.append(subset)

    # sample delta_acid_base for each peak cluster
    delta_acid_base_list = np.random.normal(-0.118, 0.204, len(filtered_signal_subset))

    return peaks_index_list, filtered_signal_subset, delta_acid_base_list


def get_peak_cluster_acid_base_list(mixture_list, processed_data_dict):
    meta_subset_dict = dict()
    for meta_name in mixture_list:
        temp_y = processed_data_dict[meta_name]
        peaks_index_list, filtered_signal_subset, delta_acid_base_list = peak_cluster_detection(temp_y)
        meta_subset_dict[meta_name] = [filtered_signal_subset, delta_acid_base_list]
    return meta_subset_dict

