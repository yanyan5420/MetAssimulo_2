import pandas as pd
import numpy as np
import pathlib

# from simulate_1D.construct_hmdb_avg_cons import get_hmdb_normal_avg_cons


def pos_normal_no_correlated(mean, std, num_replicates):
    x = np.random.normal(mean,std, num_replicates)
    if np.all(x>0):
        return x
    else:
        return pos_normal_no_correlated(mean, std, num_replicates)


def pos_normal_correlated(mean_list, cov_df, num_replicates):
    x = np.random.multivariate_normal(mean_list, cov_df, num_replicates)
    if np.all(x>0):
        return x
    else:
        return pos_normal_correlated(mean_list, cov_df, num_replicates)


def _getAplus(A):
    eigval, eigvec = np.linalg.eig(A)
    Q = np.matrix(eigvec)
    xdiag = np.matrix(np.diag(np.maximum(eigval, 0)))
    return Q*xdiag*Q.T


def _getPs(A, W=None):
    W05 = np.matrix(W**.5)
    return  W05.I * _getAplus(W05 * A * W05) * W05.I


def _getPu(A, W=None):
    Aret = np.array(A.copy())
    Aret[W > 0] = np.array(W)[W > 0]
    return np.matrix(Aret)


def nearPD(A, nit=10):
    n = A.shape[0]
    W = np.identity(n)
    deltaS = 0
    Yk = A.copy()
    for k in range(nit):
        Rk = Yk - deltaS
        Xk = _getPs(Rk, W=W)
        deltaS = Xk - Rk
        Yk = _getPu(Xk, W=W)
    return np.array(Yk)


def simulate_concentrations(table_rows, num_replicates, correlated_flag, corr_df, group_flag):
    if corr_df is None:
        correlated_flag = False

    if correlated_flag:
        mixture_list = list(map(lambda d: d["meta_name"], table_rows))
        corr_metabolites = list(corr_df.columns)
        not_corr_metabolites = list(set(mixture_list) - set(corr_metabolites))

        mean_list = []
        std_list = []
        for meta_name in corr_metabolites:
            temp_dict = list(filter(lambda d: d['meta_name'] == meta_name, table_rows))[0]
            avg_mean = float(temp_dict["mean_" + group_flag])
            avg_std = float(temp_dict["std_" + group_flag])
            mean_list.append(avg_mean)
            std_list.append(avg_std)

        # find the nearest correlation matrix
        nearest_corr_array = nearPD(np.array(corr_df), nit=10)
        cov_df = np.diag(std_list) * nearest_corr_array * np.diag(std_list)
        meta_cons_array = np.round(pos_normal_correlated(mean_list, cov_df, num_replicates), decimals=2).T
        cons_dict = dict(zip(corr_metabolites, meta_cons_array))

        for meta_name in not_corr_metabolites:
            temp_dict = list(filter(lambda d: d['meta_name'] == meta_name, table_rows))[0]
            avg_mean = float(temp_dict["mean_" + group_flag])
            avg_std = float(temp_dict["std_" + group_flag])
            sampled_cons = np.round(pos_normal_no_correlated(avg_mean, avg_std, num_replicates), decimals=2)
            cons_dict[meta_name] = np.array(sampled_cons)

    else:
        cons_dict = dict()
        mixture_list = list(map(lambda d: d["meta_name"], table_rows))
        for meta_name in mixture_list:
            temp_dict = list(filter(lambda d: d['meta_name'] == meta_name, table_rows))[0]
            avg_mean = float(temp_dict["mean_"+group_flag])
            avg_std = float(temp_dict["std_"+group_flag])
            sampled_cons = np.round(pos_normal_no_correlated(avg_mean, avg_std, num_replicates), decimals=2)
            cons_dict[meta_name] = np.array(sampled_cons)

    cons_df = pd.DataFrame.from_dict(cons_dict, orient='index')
    cons_df.reset_index(inplace=True)
    cons_df.columns = ["meta_name"] + [group_flag+"_replicate_"+str(i+1) for i in range(num_replicates)]
    return cons_df


def simulate_continuous_concentrations(table_rows, y_num_repli, sample_y):
    cons_dict_list = []
    y_sample_mean = np.mean(sample_y)
    y_sample_std = np.mean(sample_y)
    # meta_name that has correlation with Y:
    meta_with_y_list = list(filter(lambda dic: float(dic["a"]) != 0, table_rows))
    meta_not_y_list = list(filter(lambda dic: float(dic["a"]) == 0, table_rows))
    for d in meta_with_y_list:
        random_error = np.random.normal(0, float(d["std_error"]), y_num_repli)
        x_mean = float(d["mean"])
        x_std = float(d["std"])
        scale_a = (float(d["std"]) / float(y_sample_std)) * float(d["a"])
        origin_x = np.round(((sample_y - y_sample_mean) / y_sample_std - random_error) / scale_a * x_std + x_mean, 2)
        origin_x = origin_x.clip(min=0)
        print(d["meta_name"], origin_x, random_error)
        temp_dict = dict({"meta_name": d["meta_name"], "hmdb_id": d["hmdb_id"]},
                         **dict(zip(["replicate_{}".format(i + 1) for i in range(y_num_repli)], origin_x)))
        cons_dict_list.append(temp_dict)

    for di in meta_not_y_list:
        sampled_cons = np.round(pos_normal_no_correlated(float(di["mean"]), float(di["std"]), y_num_repli), 2)
        temp_dict = dict({"meta_name": di["meta_name"], "hmdb_id": di["hmdb_id"]},
                         **dict(zip(["replicate_{}".format(i + 1) for i in range(y_num_repli)], sampled_cons)))
        cons_dict_list.append(temp_dict)

    print(cons_dict_list)

    return cons_dict_list, meta_with_y_list, meta_not_y_list

