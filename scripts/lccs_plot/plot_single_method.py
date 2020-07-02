import os
import re
import numpy as np
import matplotlib.pylab as plt
from scipy.spatial import ConvexHull
from itertools import chain, product
from scipy.interpolate import interp1d

from collections import defaultdict
from plot_sigmod import parse_res
                    
def get_c(res):
    return int(res[0][0])

def get_time(res):
    return float(res[1][2])

def get_recall(res):
    return float(res[1][3])

def get_params(res):
    return (res[0][0], ) + tuple(res[0][1].items())

def plot_records(filename):
    data_dict = defaultdict(list)
    
    for record in parse_res(filename):
        params = get_params(record)
        print('record=', record)
        param_key = params[1]
        print('param_key=', params)
        
        c = get_c(record)
        t = get_time(record)
        recall = get_recall(record)
        data_dict[param_key] += [[c, t, recall]]
        
    #use marker to encode p
    #use color to encode l
    
    marker_p = {
        0 : 'o', 
        0.5: 'x', 
        1 : 's', 
        2 : '^', 
        4 : 'd', 
        8 : '*', 
#         16: '*', 
    }
    ls = [8, 16, 32, 64, 128, 256, 512]
    markers = ['o','x','s','^','d','*', 'p']
#     ls = [128]
#     colors = [np.random.rand(3, ) for l in ls]
    colors = ['b', 'g', 'r', 'c', 'm', 'y', 'tab:blue', 'tab:orange', 'tab:purple', 'tab:pink', 'tab:brown', 'tab:gray']
    #let color be random colors
    
    for (param_key, data_arr), (marker, color) in zip(data_dict.items(), product(markers, colors)):
        data_lp = np.array(data_arr)
        # print(marker, color, param_key, data_lp)
        plt.semilogy(data_lp[:, -1], data_lp[:, -2], marker=marker, label=str(param_key), color=color, markerfacecolor='none')
    plt.xlim(0, 100)
    plt.legend(ncol=6)
    plt.show()
#     color_ls = [random.randin]
    #scheme :    L, nprobe, ncheck, time, recall


# filename = 'results/Sift_srs_[02-17_13_15].out'
# filename = 'results/Sift_c2lsh_[10-16_11_19].out'
# filename = 'results/deep_qalsh_[02-18_09_52].out'
filename = 'results/Gist_mp_lccs_[10-15_04_55].out'
plot_records(filename)