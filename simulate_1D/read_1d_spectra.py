import nPYc
import pandas as pd
import numpy as np
import pathlib


def read_1d_data(file_path, sop, pulseProgram):
    base_path = pathlib.Path(__file__).resolve().parents[1]
    data_path = base_path.joinpath(file_path)
    nmrData = nPYc.NMRDataset(str(data_path), pulseProgram=pulseProgram, sop=sop, variableSize=64000)

    ppm_scale = np.array(nmrData.featureMetadata).ravel()

    meta_list = list(nmrData.sampleMetadata['Sample File Name'].apply(lambda x: x.split('/')[0]))
    data_dict = dict()
    for idx, meta_name in enumerate(meta_list):
        data = nmrData.intensityData[idx]
        data_dict[meta_name] = data

    return data_dict, ppm_scale


