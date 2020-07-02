import os
import re
import numpy as np
import matplotlib.pylab as plt
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

from scipy.spatial import ConvexHull
from itertools import chain, count
from scipy.interpolate import interp1d

from collections import defaultdict
from plot_utility import *

import glob, os

possible_datasets = ['Sift', 'Gist', 'glove', 'Trevi', 'Tiny', 'Sift10M', 'NUSW', 'deep', 'glove100', 'Msong']

dataset_labels_map = {
    'Sift': 'Sift', 
    'Gist': 'Gist', 
    'glove': 'Glove', 
    'Trevi': 'Trevi', 
    'Tiny': 'Tiny', 
    'Sift10M': 'Sift10M', 
    'Mnist784': 'Mnsit', 
    'NUSW': 'NUSW', 
    'deep': 'Deep', 
    'glove100': 'GloVe', 
    'Msong': 'Msong', 
}

method_labels_map = {'l2': {
        'lccs' : 'LCCS-LSH', 
        'e2lsh' : 'E2LSH', 
        'c2lsh' : 'C2LSH', 
        'mplsh' : 'Multi-Probe LSH', 
        'mp_lccs' : "MP-LCCS-LSH", 
        'srs' : 'SRS', 
        'qalsh' : 'QALSH'
    }, 
    'angle': {
        'lccs' : 'LCCS-LSH', 
        'e2lsh' : 'E2LSH', 
        'c2lsh' : 'C2LSH', 
        'mplsh' : 'FALCONN', 
        'mp_lccs' : "MP-LCCS-LSH"
    }
}

method_labels_map_wo_distance = {
    'lccs' : 'LCCS-LSH', 
    'e2lsh' : 'E2LSH', 
    'c2lsh' : 'C2LSH', 
    'mplsh' : 'Multi-Probe LSH (FALCONN)', 
    'mp_lccs' : "MP-LCCS-LSH"
}

distance_labels_map = {
    'l2': 'Euclidean', 
    'angle': 'Angular'
}

# distances = ['angle']
# l2
distance = 'l2'
methods = ['lccs', 'mp_lccs', 'e2lsh', 'mplsh', 'c2lsh', 'srs', 'qalsh']

# angle
# distance = 'angle'
# methods = ['lccs', 'mp_lccs', 'e2lsh', 'mplsh', 'c2lsh']

datasets = ['Msong', 'Sift', 'Gist', 'glove100', 'deep']
dataset_labels = [dataset_labels_map[dataset] for dataset in datasets]
# methods = ['lccs', 'e2lsh', 'mplsh', 'c2lsh']
method_labels = [method_labels_map[distance][method] for method in methods]

method_colors = ['red', 'purple', 'blue', 'green', 'darkgoldenrod', 'c', 'y']
method_markers = ['o', '^', 's', 'd', '*', 'p', 'x']

result_folder_path = 'results'

def get_latest_res(file_prefix):
    pat = re.compile(r'%s_\[.*\]'%(file_prefix))
    
    latestFile = None
    for file in os.listdir('results'):
        if pat.match(file):
            if latestFile is None:
                latestFile = file
            elif latestFile < file:
                latestFile = file
    if latestFile is None:
        return None
    return os.path.join('results/', latestFile)

def get_file_prefix(datasetName, methodName, distanceName):
    mmap = {
        ('e2lsh', 'l2'): 'e2lsh', 
        ('c2lsh', 'l2'): 'c2lsh', 
        ('lccs', 'l2'): 'lccs', 
        ('mplsh', 'l2'): 'mplsh_lshkit', 
        ('mp_lccs', 'l2'): 'mp_lccs', 
        
        ('e2lsh', 'angle'): 'polytope_e2', 
        ('c2lsh', 'angle'): 'polytope_c2', 
        ('lccs', 'angle'): 'polytope_lccs', 
        ('mplsh', 'angle'): 'falconn', 
        ('mp_lccs', 'angle'): 'polytope_mplccs', 
        
        ('srs', 'l2'): 'srs', 
        ('qalsh', 'l2'): 'qalsh'
    }
    return '%s_%s'%(datasetName, mmap[(methodName, distanceName)])

# e2lsh  r=700.000000, K=10 L=100
# Indexing Time: 5.903143 Seconds
# memory-usage: 17223505864 (bytes)
# 
# check_k=64
# 1    1.000000    0.068175    51.099998
# 2    1.044417    0.068178    46.349998
# 5    1.080480    0.069063    39.419998
# 10    1.108910    0.070709    34.900002
# 
# check_k=128
# 1    1.000000    0.082642    66.699997
# 2    1.005960    0.082610    63.000000
# 5    1.037351    0.083090    54.580002
# 10    1.060527    0.084648    49.049999

def parse_res(filename, chosen_top_ks=[10]):
    setting_pattern = re.compile(r'\S+\s+.*=.*')
    
    setting_l = re.compile(r'.*(L)=(\d+).*')
    setting_p = re.compile(r'.*(p)=(\d+).*')
    setting_r = re.compile(r'.*(r)=(\d+\.\d+).*')
    setting_k = re.compile(r'.*(K)=(\d+).*')
    setting_threshold = re.compile(r'.*(threshold)=(\d+).*')
    param_settings = [setting_l, setting_p, setting_r, setting_k, setting_threshold]
    
    
    indexing_time_pattern = re.compile(r'Indexing Time: (\d+\.\d+)')
    memory_usage_pattern = re.compile(r'memory-usage: (\d+).*')
    check_k_pattern = re.compile(r'check_k=(\d+)')
    records_pattern = re.compile(r'(\d+)\s*(\d+\.\d+)\s*(\d+\.\d+)\s*(\d+\.\d+)')
    est_memory_usage_pattern = re.compile(r'Estimated Memory Usage: (\d+) Bytes')
    
    with open(filename, 'r') as f:
        for line in f:
            res = setting_pattern.match(line)
            if res:
                params = {}
                for param_setting in param_settings:
                    tmp_res = param_setting.match(line)
                    if tmp_res is not None:
#                         print(tmp_res.groups())
                        params[tmp_res.group(1)] = tmp_res.group(2)
            res = indexing_time_pattern.match(line)
            if res:
                indexing_time = float(res.group(1))
#                 print('idx_time=', indexing_time)
            res = memory_usage_pattern.match(line)
            if res:
                memory_usage = int(res.group(1))
#                 print('idx_time=', indexing_time)
            res = check_k_pattern.match(line)
            if res:
                check_k = int(res.group(1))
#                 print('check_k=', check_k)
            res = est_memory_usage_pattern.match(line)
            if res:
                memory_usage = int(res.group(1))
            
            res = records_pattern.match(line)
            if res:
                top_k = int(res.group(1))
                avg_ratio = float(res.group(2))
                qtime_in_ms = float(res.group(3))
                recall = float(res.group(4))
#                 print(top_k, avg_ratio, qtime_in_ms, recall)
                
                if top_k in chosen_top_ks:
                    yield ((check_k, params), (top_k, avg_ratio, qtime_in_ms, recall, indexing_time, memory_usage))
                    
def getratio(res):
    return res[1]   
def getrecall(res):
    return res[3]   
def gettime(res):
    return res[2]
def getindexingtime(res):
    return res[4]
def getindexsize(res):
    return res[5]
def get_l(res):
    return int(res[0][1]['L'])
def get_c(res):
    return int(res[0][0])
def get_time(res):
    return float(res[1][2])
def get_recall(res):
    return float(res[1][3])
def get_p(res):
    return int(res[0][1]['p'])

#fine the best parameter setting at given recall level
def best_config_at_recall_level(xys, settings, ress, recallThreshold):
    besttime = 1e9
    bestsetting = settings[0]
    bestres = ress[0]
    
    for xy, setting, res in zip(xys, settings, ress):
        recall = xy[1]
        time = xy[0]
        
        if recall > recallThreshold and time < besttime:
            besttime = time
            bestsetting = setting
            bestres = res
    return bestsetting, bestres

def config_more_than_recall_level(xys, settings, ress, recallThreshold):
    ret_settings = []
    ret_ress = []
    ret_times = []
    
    for xy, setting, res in zip(xys, settings, ress):
        recall = xy[1]
        time = xy[0]
        
        if recall > recallThreshold:
#             print(recall, time, setting, res)
            ret_times += [time]
            ret_settings += [setting]
            ret_ress += [res]
    return ret_times, ret_settings, ret_ress

def lower_bound_curve(xys, extra_pnts = [[0, 0]]):
#     print(xys, extra_pnts)
#     xys = np.append(xys, extra_pnts, axis=0)
    
    eps = np.random.normal(size=xys.shape) * 1e-4
    xys += eps
#     print(xys)
#     xys = np.array(sorted(xys, key=lambda x:x[1]) )
#     print(xys)
    hull = ConvexHull(xys)
#     print(hull.vertices)
    
    hull_vs = xys[hull.vertices]
    
    v1s = []
    maxv0 = [-1, -1]
    for v0, v1 in zip(hull_vs, chain(hull_vs[1:], hull_vs[:1])):
#         print(v0, v1)
        if v0[1] > v1[1] and v0[0] > v1[0]:
#             plt.semilogy([v0[1], v1[1]], [v0[0], v1[0]],  'k-')
            v1s = np.append(v1s, v1, axis=-1)
            if v0[1] > maxv0[1]:
                maxv0 = v0
                
    
    vs = np.array(np.append(maxv0, v1s)).reshape(-1, 2)
#     print(vs)
#     plt.semilogy(vs[:, 1], vs[:, 0], 'k-')
        
    f = interp1d(vs[:, 1], vs[:, 0])
    
    minx = np.min(vs[:, 1])+1e-6
    maxx = np.max(vs[:, 1])-1e-6
    x = np.arange(minx, maxx, 1)
    y = list(map(f, x))
#     print(x, y)
#     plt.semilogy(x, y, 'k-')
    return x, y

def lower_bound_curve2(xs, ys):
#     print('xs, ys', xs, ys)
    xys = np.zeros(shape=(len(xs), 2))
    xys[:, 0] = xs
    xys[:, 1] = ys
    
    if len(xs)>2 and xs[-1]>0:
        hull = ConvexHull(xys)
        
        hull_vs = xys[hull.vertices]
        ret_vs = []
        
#         print("hull_vs: ", hull_vs)
        
        pflg = False
        for v0, v1, v2 in zip(chain(hull_vs[-1:], hull_vs[:-1]), hull_vs, chain(hull_vs[1:], hull_vs[:1])):
    #         print(v0, v1)
            if v0[0] < v1[0]:
    #             plt.semilogy([v0[1], v1[1]], [v0[0], v1[0]],  'k-')
                ret_vs = np.append(ret_vs, v1, axis=-1)
            elif v1[0] < v2[0]:
                ret_vs = np.append(ret_vs, v1, axis=-1)
        ret_vs = ret_vs.reshape((-1, 2))
        ret_vs = np.array(sorted(ret_vs, key=lambda x:x[0]) )
        return ret_vs
    return xys

    
def plot_indexing_phase(datasets, methods, method_labels, fig_width=0.55+3.333*len(datasets), fig_height=6.2):
    plt_helper = PlotHelper(plt, fig_width, fig_height)
    plt_helper.plot_subplots_adjust()
    n_datasets = len(datasets)
    
    recall_level = 50
#     recall_level = 85
    
#     y_lim_dataset = {
#         ('Mnist784', ):(0.1, 10),
#     }
#     y_ticks_dataset = {
#         'Mnist784':(0.1, 10),
#     }
    
    for di, (dataset, dataset_label) in enumerate(zip(datasets, dataset_labels)):
        ax_size = plt.subplot(2, n_datasets, di+1)
        plt.xlabel('Index size (GB)')
        
        plt.title(dataset_label)
        ax_time = plt.subplot(2, n_datasets, n_datasets+di+1)
#         plt.xlim(0, 100)
        plt.xlabel('Indexing time (s)')
            
#             if dataset in y_lim_dataset:
#                 ax_size.set_ylim(y_lim_dataset[dataset])
#                 ax_time.set_ylim(y_lim_dataset[dataset])
#             if dataset in x_lim_size_dataset:
#                 ax_size.set_xlim(x_lim_size_dataset[dataset])
#             if dataset in x_lim_time_dataset:
#                 ax_time.set_xlim(x_lim_time_dataset[dataset])
                
            
#         plt.ylim(ymin=0)
        
        miny = 1e9
        maxy = -1e9
        
        for method_idx, method, method_label, method_color, method_marker in zip(count(), methods, method_labels, method_colors, method_markers):
            filename_prefix = get_file_prefix(dataset, method, distance)
            filename = get_latest_res(filename_prefix)
            if filename is None:
                continue
            print('-------------', filename, '------------')
            xys = []
            settings = []
            ress = []
            
            index_timesize_dict = defaultdict(list)
            for setting, res in parse_res(filename):
                qtime = gettime(res)
                qrecall = getrecall(res)
                index_time = getindexingtime(res)
                index_size = getindexsize(res)
                
                index_timesize_dict[(index_time, index_size)] += [[qrecall, qtime]]

#             print(xys)

            index_times, index_sizes, qtimes_at_50recall = [], [], []
            for (index_time, index_size), qrecall_times in index_timesize_dict.items():
#                 print(index_size, qrecall_times)
                qrecall_times = np.array(qrecall_times+[[0, 0]])
                qrecalls = qrecall_times[:, 0]
                qtimes = qrecall_times[:, 1]
                
                if np.max(qrecalls) > recall_level:                    
                    
#                     print(qrecalls, qtimes)
                    f = interp1d(qrecalls, qtimes)
                    time_at_50recall = f(recall_level)
                
#                     print('iit', index_time, index_size, time_at_50recall)
                    index_times += [index_time]
                    index_sizes += [index_size]
                    qtimes_at_50recall += [time_at_50recall]
            
            index_times = np.array(index_times)
            index_sizes = np.array(index_sizes)
            qtimes_at_50recall = np.array(qtimes_at_50recall)
            
            
            index_size_qtimes = lower_bound_curve2(index_sizes/1e9, qtimes_at_50recall)
            if len(index_size_qtimes)>0:
                print('min_qtime=', np.min(index_size_qtimes[:, 1]))
                
                ax_size.semilogy(index_size_qtimes[:, 0], index_size_qtimes[:, 1], '-', color=method_color, marker=method_marker, label=method_label if di==0 else "", 
                            markerfacecolor='none', markersize=10)
                
                
                index_time_qtimes = lower_bound_curve2(index_times, qtimes_at_50recall)
                print(method, index_time_qtimes)
                ax_time.semilogy(index_time_qtimes[:, 0], index_time_qtimes[:, 1], '-', color=method_color, marker=method_marker, label="", 
                            markerfacecolor='none', markersize=10, zorder=len(methods)-method_idx)
                
                miny = min(miny, np.min(index_time_qtimes[:, 1]) )
                maxy = max(maxy, np.max(index_time_qtimes[:, 1]) ) 
                miny = min(miny, np.min(index_size_qtimes[:, 1]) )
                maxy = max(maxy, np.max(index_size_qtimes[:, 1]) ) 
                
        plt_helper.set_y_axis_close(ax_time, miny, maxy)
        plt_helper.set_y_axis_close(ax_size, miny, maxy)
        if di==0:
            ax_size.set_ylabel('Query time (ms)')
            ax_time.set_ylabel('Query time (ms)')
                
#             print('xys_=', xys_)
#             if len(xys_) > 1:
#                 query_time, indexing_time = lower_bound_curve(xys_)
#             else:
#                 query_time, indexing_time = xys_[:, 0], xys_[:, 1]
#             
# #             ax.plot(xys[:, 1], xys[:, 0], '.', color=method_color, marker=method_marker, label=method_label)
# #             ax.semilogy(xys[:, 1], xys[:, 0], '.', color=method_color, marker=method_marker, label=method_label)
#             ax.semilogy(indexing_time, query_time, '-', color=method_color, marker=method_marker, label=method_label, markevery=10, 
#                         markerfacecolor='none', markersize=10)
    
#     plt.figlegend(fontsize=16, bbox_to_anchor=(0.07,0.85,0.75,0.2), loc="center",
#                 mode="expand", borderaxespad=0, ncol=len(methods))
    plt_helper.plot_fig_legend(ncol=len(methods))    
    plt_helper.plot_and_save('saved_figure/indexing_time_recall_%s'%distance)
    
    
                    
def plot_methods(datasets, methods, method_labels, fig_width = 3.*len(datasets), fig_height = 2.7 + 0.8):
    plt_helper = PlotHelper(plt, fig_width, fig_height)
    plt_helper.plot_subplots_adjust()
    
    for di, (dataset, dataset_label) in enumerate(zip(datasets, dataset_labels)):
        ax = plt.subplot(1, len(datasets), di+1)
        plt.title(dataset_label)
        plt.xlim(0, 100)
        plt.xlabel('Recall (%)')
        if di==0:
            plt.ylabel('Query time (ms)')
#         plt.ylim(ymin=0)

        miny = 1e9
        maxy = -1e9
        
        for method_idx, method, method_label, method_color, method_marker in zip(count(), methods, method_labels, method_colors, method_markers):
            filename_prefix = get_file_prefix(dataset, method, distance)
            filename = get_latest_res(filename_prefix)
            if filename is None:
                continue
            print(filename)
            xys = []
            settings = []
            ress = []
            for setting, res in parse_res(filename):
#                 print(setting, res)
                xys += [[gettime(res), getrecall(res)]]
                ress += [res]
                settings += [setting]
            xys = np.array(xys)
#             print(xys)

            lower_x, lower_y = lower_bound_curve(xys)
#             print(lower_x, lower_y)
            miny = min(miny, np.min(lower_y) )
            maxy = max(maxy, np.max(lower_y) ) 
            ax.semilogy(lower_x, lower_y, '-', color=method_color, marker=method_marker, label=method_label if di==0 else "", markevery=10, 
                        markerfacecolor='none', markersize=7, zorder=len(methods)-method_idx)
            
#         print(miny, maxy)
        plt_helper.set_y_axis_log10(ax, miny, maxy)
    
    plt_helper.plot_fig_legend(ncol=len(methods))
    plt_helper.plot_and_save('saved_figure/time_recall_%s'%distance)
    plt.show()
                    
def plot_methods_recall_ratio(datasets, methods, method_labels):
    n_datasets = len(datasets)
    fig = plt.figure(figsize=(0.55+3.333*len(datasets), 6.))
    plt.rcParams.update({'font.size': 13})
    plt.subplots_adjust(bottom=0.15, top=0.85, hspace=0.35)
    
    for di, (dataset, dataset_label) in enumerate(zip(datasets, dataset_labels)):
        ax_recall = plt.subplot(2, n_datasets, di+1)
        plt.xlabel('Recall (%)')
        plt.title(dataset_label)
        if di==0:
            plt.ylabel('Query time (ms)')
        plt.xlim(0, 100)
        
        plt.title(dataset_label)
        ax_ratio = plt.subplot(2, n_datasets, n_datasets+di+1)
        plt.xlabel('Ratio')
        ax_ratio.set_xlim(1., 1.5)
        ax_ratio.set_xticks([1., 1.1, 1.2, 1.3, 1.4, 1.5])
#         plt.xlim(0, 100)
        if di==0:
            plt.ylabel('Query time (ms)')
#         plt.ylim(ymin=0)
        
        for method, method_label, method_color, method_marker in zip(methods, method_labels, method_colors, method_markers):
            filename_prefix = get_file_prefix(dataset, method, distance)
            filename = get_latest_res(filename_prefix)
            if filename is None:
                continue
            print(filename)
            time_recalls = []
            time_ratios = []
            settings = []
            ress = []
            for setting, res in parse_res(filename):
#                 print(setting, res)
                time_recalls += [[gettime(res), getrecall(res)]]
                time_ratios += [[gettime(res), getratio(res)]]
                ress += [res]
                settings += [setting]
            time_recalls = np.array(time_recalls)
            time_ratios = np.array(time_ratios)
#             print(xys)
            
            lower_x, lower_y = lower_bound_curve(time_recalls)
            ax_recall.semilogy(lower_x, lower_y, '-', color=method_color, marker=method_marker, label=method_label if di==0 else "", markevery=10, 
                        markerfacecolor='none', markersize=10)
            
#             lower_x, lower_y = lower_bound_curve2(time_ratios[:, 0], time_ratios[:, 1])
            qtime_ratio_lower = lower_bound_curve2(time_ratios[:, 0], time_ratios[:, 1])
            ax_ratio.semilogy(qtime_ratio_lower[:, 1], qtime_ratio_lower[:, 0], '-', color=method_color, marker=method_marker, label="", markevery=5, 
                        markerfacecolor='none', markersize=10)
    
    plt.figlegend(fontsize=16, bbox_to_anchor=(0.07,0.85,0.75,0.2), loc="center",
                mode="expand", borderaxespad=0, ncol=len(methods))
    plt.savefig('saved_figure/time_recall_ratio_%s.png'%distance, bbox_inches='tight', pad_inches=0.27, dpi=fig.dpi)
    plt.savefig('saved_figure/time_recall_ratio_%s.eps'%distance, bbox_inches='tight', pad_inches=0.27, dpi=fig.dpi)
    plt.show()
    
if __name__ == '__main__':
    plot_methods(datasets, methods, method_labels)
    # plot_methods_recall_ratio(datasets, methods, method_labels)
    plot_indexing_phase(datasets, methods, method_labels)
#     plot_indexing_time(datasets, methods, method_labels)