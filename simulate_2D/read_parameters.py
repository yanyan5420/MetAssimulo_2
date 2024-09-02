def read_param(parameter_file_name):
    param_dict = dict()
    with open(parameter_file_name, 'r') as file:
        for line in file:
            if line.startswith('"'):
                k = line.strip('\n').split('\t')[0].strip('""')
                v = line.strip('\n').split('\t')[1].strip('""')
                param_dict[k] = v
    return param_dict
