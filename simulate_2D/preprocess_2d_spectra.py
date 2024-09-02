import pandas as pd
import numpy as np
from scipy.signal import convolve2d


def remove_water_calibration(data_dict, x_scale, y_scale, range_water, range_calibration):
    """
    range_water: a sorted tuple/list, e.g. [4.5, 5.0]
    range_calibration: a sorted tuple/list, e.g. [-0.3, 0.3]
    """
    removed_data_dict = dict()

    for meta_name, data in data_dict.items():
        temp_df = pd.DataFrame(data)
        temp_df.columns = x_scale
        temp_df.index = y_scale

        for ppm in list(temp_df.columns):
            # range for water suppression is 4.5-5.0
            if (float(range_water[0]) < ppm < float(range_water[1])) | (float(range_calibration[0]) < ppm < float(range_calibration[1])):
                temp_df.loc[:, ppm] = 0
        removed_data_dict[meta_name] = np.array(temp_df)

    return removed_data_dict


def filter_noise(removed_data_dict, threshold):
    filtered_data_dict = dict()

    for meta_name, data in removed_data_dict.items():
        # print(meta_name, " filter noise: max is ", np.max(data))
        temp_max = np.max(data) * float(threshold)
        data[np.where(data < temp_max)] = 0
        filtered_data_dict[meta_name] = data

    return filtered_data_dict


def smooth_data(filtered_data_dict, m, n):
    """
    m, n: the size of window = m*n, e.g. 3x3
    """
    win = np.ones((int(m), int(n)))

    smooth_data_dict = dict()
    for meta_name, data in filtered_data_dict.items():
        temp_data = convolve2d(data, win, mode='same', boundary='symm')
        smooth_data_dict[meta_name] = temp_data

    return smooth_data_dict


def normalize_data(smooth_data_dict):
    norm_data_dict = dict()
    for meta_name, data in smooth_data_dict.items():
        temp_data = np.array(data)
        temp_sum = np.sum(temp_data)
        norm_data = temp_data / temp_sum
        norm_data_dict[meta_name] = norm_data

        # print(meta_name, " normalize: max is ", np.max(norm_data))
        # print(norm_data.shape)

    return norm_data_dict
