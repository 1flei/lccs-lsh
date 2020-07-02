import os
import re
import numpy as np
import matplotlib.pylab as plt
from scipy.spatial import ConvexHull
from itertools import chain
from scipy.interpolate import interp1d

from collections import defaultdict
from plot_sigmod import parse_res
                    
def get_l(res):
    return int(res[0][1]['L'])

def get_p(res):
    return float(res[0][1]['p'])

def get_c(res):
    return int(res[0][0])

def get_time(res):
    return float(res[1][2])

def get_recall(res):
    return float(res[1][3])

def plot_records(filename):
    data = []
    
    for record in parse_res(filename):
        print(record)
        l = get_l(record)
        p = get_p(record)
        c = get_c(record)
        t = get_time(record)
        recall = get_recall(record)
        print(l, p, c, t, recall)
        data += [[l, p, c, t, recall]]
    data = np.array(data)
        
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
#     ls = [128]
#     colors = [np.random.rand(3, ) for l in ls]
    colors = ['b', 'g', 'r', 'c', 'm', 'y', 'tab:blue', 'tab:orange', 'tab:purple', 'tab:pink', 'tab:brown', 'tab:gray']
    #let color be random colors
    
    for p, marker in marker_p.items():
        for color, l in zip(colors, ls):
            data_lp = data[np.logical_and(data[:, 0]==l, data[:, 1]==p)]
            # print(l, p, data_lp)
            
            plt.semilogy(data_lp[:, -1], data_lp[:, -2], marker=marker, label='l=%d, p=%.1f'%(l, p), c=color, markerfacecolor='none')
    plt.xlim(0, 100)
    plt.legend(ncol=6)
    plt.show()
#     color_ls = [random.randin]
    #scheme :    L, nprobe, ncheck, time, recall


filename = 'results/mp_lccs/glove100_mp_lccs_[10-20_17_47].out'
# filename = 'results/Gist_mp_lccs_[10-15_04_55].out'
plot_records(filename)