import pandas as pd
import numpy as np


def format_data_names(norm_data_dict):
    format_norm_data_dict = dict()
    for meta_name, data in norm_data_dict.items():
        lower_meta_name = meta_name.lower()
        if lower_meta_name[1] == '_' or lower_meta_name[2] == '_':
            lower_meta_name = lower_meta_name.replace("_", "-", 1)
        lower_meta_name = ' '.join(lower_meta_name.split("_"))

        if "acid" in lower_meta_name:
            lower_meta_name = ' '.join(lower_meta_name.split("acid")).strip(' ') + ' acid'

        format_norm_data_dict[lower_meta_name] = data

    return format_norm_data_dict


def db_names_match_hmdb(format_norm_data_dict, hmdb_dict):
    db_names_ids_dict = dict()
    for name in format_norm_data_dict.keys():
        id_list = []
        for idx, name_list in hmdb_dict.items():
            if name in name_list:
                id_list.append(idx)
        db_names_ids_dict[name] = id_list

    return db_names_ids_dict


def cons_names_match_hmdb(cons_df, hmdb_dict):
    cons_names_ids_dict = dict()
    for name in cons_df.index:
        id_list = []
        for idx, name_list in hmdb_dict.items():
            if name in name_list:
                id_list.append(idx)
        cons_names_ids_dict[name] = id_list
    return cons_names_ids_dict


def db_match_cons(norm_data_dict, cons_df, hmdb_dict1):
    format_norm_data_dict = format_data_names(norm_data_dict)
    db_names_ids_dict = db_names_match_hmdb(format_norm_data_dict, hmdb_dict1)
    cons_names_ids_dict = cons_names_match_hmdb(cons_df, hmdb_dict1)

    match_data_dict = dict()
    for name, data in format_norm_data_dict.items():
        db_id_list = db_names_ids_dict[name]
        for cons_name, cons_id_list in cons_names_ids_dict.items():
            if (name == cons_name) or (set(db_id_list).intersection(set(cons_id_list))):
                match_data_dict[cons_name] = data

    return match_data_dict


def db_names_match_hmdb_names(norm_data_dict, hmdb_dict):
    format_norm_data_dict = format_data_names(norm_data_dict)
    db_names_ids_dict = db_names_match_hmdb(format_norm_data_dict, hmdb_dict)

    match_data_dict = dict()
    for name, data in format_norm_data_dict.items():
        db_id_list = db_names_ids_dict[name]
        if db_id_list != []:
            match_data_dict[name] = data

    return match_data_dict


# match input mixture list with db names
def input_match_db(mixture_list, match_data_dict, hmdb_dict):
    # format input mixtures
    format_mixture_list = []
    for meta_name in mixture_list:
        lower_meta_name = meta_name.lower()
        if lower_meta_name[1] == '_' or lower_meta_name[2] == '_':
            lower_meta_name = lower_meta_name.replace("_", "-", 1)
        lower_meta_name = ' '.join(lower_meta_name.split("_"))

        if "acid" in lower_meta_name:
            lower_meta_name = ' '.join(lower_meta_name.split("acid")).strip(' ') + ' acid'
        format_mixture_list.append(lower_meta_name)

    # get the db match hmdb ids dict
    db_names_ids_dict = db_names_match_hmdb(match_data_dict, hmdb_dict)

    # match input mixtures with db names
    mixture_match_db_list = []
    for meta_name in format_mixture_list:
        if meta_name in match_data_dict.keys():
            mixture_match_db_list.append(meta_name)
        else:
            idx_list = []
            for idx, name_list in hmdb_dict.items():
                if meta_name in name_list:
                    idx_list.append(idx)
            for db_name, db_idx_list in db_names_ids_dict.items():
                if set(idx_list).intersection(set(db_idx_list)):
                    mixture_match_db_list.append(db_name)

    return mixture_match_db_list


def input_corr_match_db(corr_df, match_data_dict, hmdb_dict):
    # format input correlation names
    format_corr_names_dict = dict()
    for meta_name in list(corr_df.columns):
        lower_meta_name = meta_name.lower()
        if lower_meta_name[1] == '_' or lower_meta_name[2] == '_':
            lower_meta_name = lower_meta_name.replace("_", "-", 1)
        lower_meta_name = ' '.join(lower_meta_name.split("_"))

        if "acid" in lower_meta_name:
            lower_meta_name = ' '.join(lower_meta_name.split("acid")).strip(' ') + ' acid'
        format_corr_names_dict[meta_name] = lower_meta_name

    # get the db match hmdb ids dict
    db_names_ids_dict = db_names_match_hmdb(match_data_dict, hmdb_dict)

    # match input corr_names with db names
    corr_match_db_dict = dict()
    for orig_name, format_name in format_corr_names_dict.items():
        if format_name in match_data_dict.keys():
            corr_match_db_dict[orig_name] = format_name
        else:
            idx_list = []
            for idx, name_list in hmdb_dict.items():
                if format_name in name_list:
                    idx_list.append(idx)
            for db_name, db_idx_list in db_names_ids_dict.items():
                if set(idx_list).intersection(set(db_idx_list)):
                    corr_match_db_dict[orig_name] = db_name

    corr_df.rename(columns=corr_match_db_dict, index=corr_match_db_dict, inplace=True)
    return corr_df


def format_input_mixture(mixture_list):
    # format input mixtures
    format_mixture_list = []
    for meta_name in mixture_list:
        lower_meta_name = meta_name.lower()
        if lower_meta_name[1] == '_' or lower_meta_name[2] == '_':
            lower_meta_name = lower_meta_name.replace("_", "-", 1)
        lower_meta_name = ' '.join(lower_meta_name.split("_"))

        if "acid" in lower_meta_name:
            lower_meta_name = ' '.join(lower_meta_name.split("acid")).strip(' ') + ' acid'
        format_mixture_list.append(lower_meta_name)

    return format_mixture_list


def input_cons_match_db(cons_df, match_data_dict, hmdb_dict):
    # format input concentration df names
    format_cons_names_dict = dict()
    for meta_name in list(cons_df.iloc[:,0]):
        lower_meta_name = meta_name.lower()
        if lower_meta_name[1] == '_' or lower_meta_name[2] == '_':
            lower_meta_name = lower_meta_name.replace("_", "-", 1)
        lower_meta_name = ' '.join(lower_meta_name.split("_"))

        if "acid" in lower_meta_name:
            lower_meta_name = ' '.join(lower_meta_name.split("acid")).strip(' ') + ' acid'
        format_cons_names_dict[meta_name] = lower_meta_name

    # get the db match hmdb ids dict
    db_names_ids_dict = db_names_match_hmdb(match_data_dict, hmdb_dict)

    # match input corr_names with db names
    cons_match_db_dict = dict()
    for orig_name, format_name in format_cons_names_dict.items():
        if format_name in match_data_dict.keys():
            cons_match_db_dict[orig_name] = format_name
        else:
            idx_list = []
            for idx, name_list in hmdb_dict.items():
                if format_name in name_list:
                    idx_list.append(idx)
            for db_name, db_idx_list in db_names_ids_dict.items():
                if set(idx_list).intersection(set(db_idx_list)):
                    cons_match_db_dict[orig_name] = db_name
    cons_df.iloc[:, 0] = list(cons_match_db_dict.values())
    return cons_df