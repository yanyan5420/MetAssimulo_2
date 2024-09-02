import pandas as pd
import numpy as np
import math
from scipy import interpolate
from itertools import groupby
from operator import itemgetter


def remove_water_calibration(data_dict, ppm_scale, range_water, range_calibration):
    """
    range_water: a sorted tuple/list, e.g. [4.5, 5.0]
    range_calibration: a sorted tuple/list, e.g. [-0.3, 0.3]
    """
    removed_data_dict = dict()

    for meta_name, data in data_dict.items():
        temp_df = pd.DataFrame(data)
        temp_df.index = ppm_scale

        for ppm in list(temp_df.index):
            # range for water suppression is 4.5-5.0
            if (range_water[0] < ppm < range_water[1]) | (range_calibration[0] < ppm < range_calibration[1]):
                temp_df.loc[ppm, :] = 0
        removed_data_dict[meta_name] = np.array(temp_df).ravel()

    return removed_data_dict


# baseline correction
def baseline_correction(removed_data_dict, bins, thres_perc):
    corrected_data_dict = dict()

    for meta_name, data in removed_data_dict.items():
        temp_y = data
        binsize = math.floor(len(temp_y) / bins)

        bin_0 = 0
        bin_1 = binsize
        bin_median_list = [0, ]
        x_list = [0, ]

        while bin_1 <= len(temp_y):
            temp_bin = temp_y[bin_0:bin_1]
            temp_median = np.median(temp_bin)
            temp_x = (bin_0 + bin_1) / 2
            bin_median_list.append(temp_median)
            x_list.append(temp_x)

            bin_0 = bin_0 + binsize
            bin_1 = bin_1 + binsize

        threshold = max(temp_y) * thres_perc
        in_threshold_index_list = list(np.where(temp_y < threshold)[0])

        new_y = interpolate.CubicSpline(x_list, bin_median_list)(range(len(in_threshold_index_list)))
        corrected_y = temp_y[in_threshold_index_list] - new_y
        temp_y[in_threshold_index_list] = corrected_y
        corrected_data_dict[meta_name] = temp_y

    return corrected_data_dict


# smooth_noises around baseline
def smooth_each_spectrum(y, thres_perc):
    temp_y = y.copy()
    threshold = max(temp_y) * thres_perc
    in_threshold_index_list = list(np.where(temp_y < threshold)[0])

    temp_smooth_df = pd.DataFrame(temp_y)
    temp_smooth_df.loc[:, 'smooth_y'] = temp_y

    in_threshold_subset = []
    for k, g in groupby(enumerate(in_threshold_index_list), lambda ix: ix[0] - ix[1]):
        in_threshold_subset.append(list(map(itemgetter(1), g)))

    for subset in in_threshold_subset:
        temp_y_subset = temp_y[subset]
        temp_smooth_subset = np.zeros(len(temp_y_subset))
        forward_window_size = round(len(temp_y_subset) / 3)

        for idx, value in enumerate(temp_y_subset):
            if idx < forward_window_size:
                start = idx - idx
                end = idx + idx + 1
                temp_l = temp_y_subset[start:end]
                temp_mean = np.mean(temp_l)
            elif (len(temp_y_subset) - 1 - idx) < forward_window_size:
                start = idx - (len(temp_y_subset) - 1 - idx)
                end = idx + (len(temp_y_subset) - 1 - idx) + 1
                temp_l = temp_y_subset[start:end]
                temp_mean = np.mean(temp_l)
            else:
                start = idx - forward_window_size
                end = idx + forward_window_size + 1
                temp_l = temp_y_subset[start:end]
                temp_mean = np.mean(temp_l)

            temp_smooth_subset[idx] = temp_mean

        temp_smooth_df.loc[subset, 'smooth_y'] = temp_smooth_subset

    smooth_y = np.array(temp_smooth_df['smooth_y']).ravel()
    return smooth_y


def smooth_spectra(corrected_data_dict, thres_perc):
    smooth_data_dict = dict()

    for meta_name, data in corrected_data_dict.items():
        temp_y = data
        smooth_y = smooth_each_spectrum(temp_y, thres_perc)
        smooth_data_dict[meta_name] = smooth_y

    return smooth_data_dict


# normalise spectrum for each metabolite in the mixture list
def norm_spectra(smooth_data_dict):
    norm_data_dict = dict()
    for meta_name, data in smooth_data_dict.items():
        temp_sum = np.sum(data)
        norm_data = data / temp_sum
        # temp_max = np.max(data)
        # norm_data = data / temp_max
        norm_data_dict[meta_name] = norm_data
    return norm_data_dict

