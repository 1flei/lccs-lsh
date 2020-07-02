import os
import re
import numpy as np
import matplotlib.pylab as plt
from scipy.spatial import ConvexHull
from itertools import chain, product
from scipy.interpolate import interp1d

from collections import defaultdict
from plot_sigmod import parse_res, get_file_prefix, get_latest_res
                    
def get_time(res):
    return float(res[1][2])

def get_recall(res):
    return float(res[1][3])

def get_k(res):
    return int(res[1][0])

def get_params(res):
    return (res[0][0], ) + tuple(res[0][1].items())

def get_params_raw(res):
    return res[0]

def get_ratio(res):
    return float(res[1][1])

def plot_records(filename):
    ks = [1, 2, 5, 10, 20, 50, 100]
    
    data = []
    
    param_dict = {}
    param_dict_inv = {}
    def paramToIdx(param):
        if param not in param_dict:
            param_dict[param] = len(param_dict)
            param_dict_inv[param_dict[param]] = param
        return param_dict[param]
        
    
    for record in parse_res(filename, ks):
        params = get_params(record)
        print(record)
        
        t = get_time(record)
        recall = get_recall(record)
        ratio = get_ratio(record)
        k = get_k(record)
        print(params, t, recall)
        data += [[paramToIdx(params), k, t, ratio, recall]]
    data = np.array(data)
    print(data)
    
    markers = ['o', 's', '^', 'd', '*', 'v', '<', '>', 'x']
    colors = ['b', 'g', 'r', 'c', 'm', 'y', 'tab:blue', 'tab:orange', 'tab:purple', 'tab:pink', 'tab:brown', 'tab:gray']    #12
    linestyles = ['-', '--', '-.']
    #let color be random colors
    
    for pid, (marker, color, linestyle) in enumerate(product(markers, colors, linestyles) ):
        data_pid = data[data[:, 0]==pid]
        print(data_pid)
#         print(param_dict_inv[pid])
        if len(data_pid)==0:
            break
        
        plt.semilogx(data_pid[:, 1], data_pid[:, -1], marker=marker, label='%s'%(param_dict_inv[pid], ), c=color, ls=linestyle, markerfacecolor='none')
    plt.xlim(1, 100)
    plt.legend(ncol=4, fontsize=6, framealpha=0.5)
    plt.show()
#     color_ls = [random.randin]
    #scheme :    L, nprobe, ncheck, time, recall

chosen_datasets = ['Sift']
methods = ['c2lsh']
# methods = ['lccs', 'e2lsh', 'mp_lccs', 'mplsh', 'c2lsh']

# from mpl_toolkits.mplot3d import Axes3D

if __name__ == '__main__':
    distances = ['l2', 'angle']
    
    for method in methods:
        for dataset in chosen_datasets:
            for distance in distances:
                data = []
                filename_prefix = get_file_prefix(dataset, method, distance)
                filename = get_latest_res(filename_prefix)
                print(filename, method, dataset, distance)

                plot_records(filename)

