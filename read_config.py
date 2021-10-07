from config import param
import pandas as pd


def read_config(conf):
    params = param.copy()
    l = conf.readlines()
    for i in l:
        if i == '\n':
            continue
        i = i.split(' ')
        name_param = i[0].lower()
        value_param = i[1]
        if len(i) == 5:
            value_param = [float(i[1]), float(i[2]), float(i[3]), float(i[4])]
        else:
            value_param = float(value_param)
        if name_param in params.keys():
            params[name_param] = value_param
        if name_param == 'start_year':
            start_year = int(value_param)
        if name_param == 'end_year':
            end_year = int(value_param)
    years = list(range(start_year, end_year + 1, 1))
    shedule = pd.DataFrame({'year': [], 'nsk': [], 'qmax': []})
    for i in l:
        if i == '\n':
            continue
        i = i.split(' ')
        name_param = i[0].lower()
        if name_param == 'year':
            shedule = shedule.append({'year': int(i[1]), 'nsk': int(i[3]), 'qmax': int(i[5])}, ignore_index=True)
    shedule = pd.merge(pd.DataFrame(years, columns=['year']), shedule, how='left', on='year')
    shedule.fillna(method='ffill', inplace=True)
    return params, shedule

#con = open('config.txt', 'r')
#p, sh = read_config(con)
#print(sh)