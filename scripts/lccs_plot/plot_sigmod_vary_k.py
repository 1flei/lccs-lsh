from plot_sigmod import *
import matplotlib.pyplot as plt
from plot_one_method_multiple_parameters import get_k, get_params, get_ratio, get_recall, get_time, get_params_raw

import numpy as np
from itertools import product

from plot_utility import *

# from mpl_toolkits.mplot3d import Axes3D

# chsen_datasets = ['Sift', 'deep']
# methods = ['lccs', 'mp_lccs', 'e2lsh', 'mplsh', 'c2lsh']
# distances = ['l2', 'angle']

chosen_params = {
    #Sift l2
    ('Sift', 'lccs', 'l2'): (16, ('L', 256)), 
#     ('Sift', 'mp_lccs', 'l2'): (1, ('L', 256), ('p', 1)), 
    ('Sift', 'mp_lccs', 'l2'): (16, ('L', 256), ('p', 0)), 
    ('Sift', 'e2lsh', 'l2'): (16384, ('L', 64), ('K', 6)), 
#     ('Sift', 'c2lsh', 'l2'): (3, ('L', 8)), 
    ('Sift', 'c2lsh', 'l2'): (16384, ('L', 64), ('threshold', 8) ), 
    ('Sift', 'mplsh', 'l2'): (8, ('L', 16), ('K', 8)), 
    
    #Sift angle
#     ('Sift', 'lccs', 'angle'): (2, ('L', 256)), 
    ('Sift', 'lccs', 'angle'): (2, ('L', 256)), 
#     ('Sift', 'mp_lccs', 'angle'): (1, ('L', 128), ('p', 1)), 
    ('Sift', 'mp_lccs', 'angle'): (2, ('L', 256), ('p', 0)), 
    ('Sift', 'e2lsh', 'angle'): (4096, ('L', 64), ('K', 6)), 
#     ('Sift', 'c2lsh', 'angle'): (4, ('L', 32)), 
    ('Sift', 'c2lsh', 'angle'): (4096, ('L', 16), ('threshold', 6)), 
    ('Sift', 'mplsh', 'angle'): (16, ('L', 16), ('K', 6)), 
    
    
    #Deep l2
    ('deep', 'lccs', 'l2'): (16, ('L', 512)), 
    ('deep', 'mp_lccs', 'l2'): (4, ('L', 256), ('p', 1)), 
    ('deep', 'e2lsh', 'l2'): (65536, ('L', 64), ('K', 5)), 
    ('deep', 'c2lsh', 'l2'): (5, ('L', 16)), 
    ('deep', 'mplsh', 'l2'): (32, ('L', 32), ('K', 9)), 
    
    #Deep angle
    ('deep', 'lccs', 'angle'): (1, ('L', 512)), 
    ('deep', 'mp_lccs', 'angle'): (1, ('L', 128), ('p', 1)), 
    ('deep', 'e2lsh', 'angle'): (16384, ('L', 128), ('K', 4)), 
    ('deep', 'c2lsh', 'angle'): (6, ('L', 64)), 
    ('deep', 'mplsh', 'angle'): (8, ('L', 64), ('K', 4)), 
    
    
    #Gist l2
    ('Gist', 'lccs', 'l2'): (64, ('L', 512)), 
    ('Gist', 'mp_lccs', 'l2'): (16, ('L', 256), ('p', 1)), 
    ('Gist', 'e2lsh', 'l2'): (65536, ('L', 64), ('K', 6)), 
    ('Gist', 'c2lsh', 'l2'): (65536, ('L', 128), ('threshold', 10)), 
    ('Gist', 'mplsh', 'l2'): (32, ('L', 32), ('K', 9)), 
    
    #Gist angle
    ('Gist', 'lccs', 'angle'): (16, ('L', 512)), 
    ('Gist', 'mp_lccs', 'angle'): (4, ('L', 256), ('p', 1)), 
    ('Gist', 'e2lsh', 'angle'): (65546, ('L', 64), ('K', 3)), 
    ('Gist', 'c2lsh', 'angle'): (16384, ('L', 32), ('threshold', 10)), 
    ('Gist', 'mplsh', 'angle'): (4, ('L', 64), ('K', 5)), 
}


def is_param_chosen(param, dataset, method, distance):
    if (dataset, method, distance) not in chosen_params:
        return False
    
    chosen_param = chosen_params[(dataset, method, distance)]
    
    if chosen_param[0] != int(param[0]):
        return False
    
    for k, v in chosen_param[1:]:
#         print(int(param[1][k]), v)
        if int(param[1][k]) != v:
            return False
    return True

def plot_method_varying_k(chosen_datasets, distances, methods, chosen_top_ks=[1, 2, 5, 10, 20, 50, 100], fig_width=6, fig_height=8, legend_col=5):
#     chosen_top_ks = [1, 2, 5, 10]
    plt_n_cols = len(chosen_datasets) * len(distances)
    plt_n_rows = 3
    
#     fig = plt.figure(figsize=(0.55+3.33*plt_n_cols, 8.))
#     plt.rcParams.update({'font.size': 13})
#     plt.subplots_adjust(bottom=0.15, top=0.85, hspace = 0.35)
    
    plt_helper = PlotHelper(plt, fig_width, fig_height)
    plt_helper.plot_subplots_adjust(top_space=1.6, left_space=0.9, wspace=0.3)
    
    for row_id, (distance, dataset) in enumerate(product(distances, chosen_datasets)):
        ax_recall = plt.subplot(plt_n_rows, plt_n_cols, row_id+1)
        ax_ratio  = plt.subplot(plt_n_rows, plt_n_cols, row_id+1+plt_n_cols)
        ax_time   = plt.subplot(plt_n_rows, plt_n_cols, row_id+1+plt_n_cols*2)
        
        plt.xlabel('$k$')
        if row_id%plt_n_cols==0:
            ax_recall.set_ylabel(r'Recall (%)')
            ax_ratio.set_ylabel(r'Ratio')
            ax_time.set_ylabel(r'Query Time (ms)')
            
#         ax_recall.set_title('%s-%s'%(dataset_labels_map[dataset], distance_labels_map[distance]))
        if len(distances)==1:
            ax_recall.set_title('%s'%dataset_labels_map[dataset])
        else:
            ax_recall.set_title('%s-%s'%(dataset_labels_map[dataset], distance_labels_map[distance]))
            
        mink, maxk = min(chosen_top_ks), max(chosen_top_ks)
        ax_recall.set_xlim(mink, maxk)
        ax_ratio.set_xlim(mink, maxk)
#         ax_ratio.set_ylim(0.995, 1.07)
        ax_time.set_xlim(mink, maxk)
        
            
#             if dataset in y_lim_dataset:
#                 ax_size.set_ylim(y_lim_dataset[dataset])
#                 ax_time.set_ylim(y_lim_dataset[dataset])
#             if dataset in x_lim_size_dataset:
#                 ax_size.set_xlim(x_lim_size_dataset[dataset])
#             if dataset in x_lim_time_dataset:
#                 ax_time.set_xlim(x_lim_time_dataset[dataset])
        
        mint = 1e9
        maxt = -1e9
        minratio = 1e9
        maxratio = -1e9
        for method_idx, method, method_color, method_marker in zip(count(), methods, method_colors, method_markers):  
            data = []
            
#             method_label = method_labels_map[distance][method]
            method_label = method_labels_map_wo_distance[method]
            
            filename_prefix = get_file_prefix(dataset, method, distance)
            filename = get_latest_res(filename_prefix)
            print(filename, method, dataset, distance)
            
            if filename is None:
                continue
            
            for record in parse_res(filename, chosen_top_ks=chosen_top_ks):
                k = get_k(record)
                t = get_time(record)
                recall = get_recall(record)
                ratio = get_ratio(record)
                param = get_params_raw(record)
                
                if not is_param_chosen(param, dataset, method, distance):
#                     if method=='mp_lccs' and distance=='l2':
#                         print(param, dataset, method, distance)
                    continue
#                 print(param)
                
                
                data += [[k, t, recall, ratio]]
                
            data = np.array(data)
            print(data)
            
            if len(data)>0:
                ax_recall.semilogx(data[:, 0], data[:, 2], label=method_label, color=method_color, marker=method_marker, markerfacecolor='none', markersize=7, zorder=len(methods)-method_idx)
                ax_ratio.semilogx(data[:, 0], data[:, 3], label=method_label, color=method_color, marker=method_marker, markerfacecolor='none', markersize=7, zorder=len(methods)-method_idx)
                ax_time.loglog(data[:, 0], data[:, 1], label=method_label, color=method_color, marker=method_marker, markerfacecolor='none', markersize=7, zorder=len(methods)-method_idx)
    
                mint = min(mint, np.min(data[:, 1]) )
                maxt = max(maxt, np.max(data[:, 1]) ) 
                minratio = min(minratio, np.min(data[:, 3]) )
                maxratio = max(maxratio, np.max(data[:, 3]) )
        
        ax_recall.set_xticks(chosen_top_ks, minor=False)
        ax_ratio.set_xticks(chosen_top_ks, minor=False)
        ax_time.set_xticks(chosen_top_ks, minor=False)
        labels = ax_recall.set_xticklabels(['%d'%k for k in chosen_top_ks], minor=False, rotation=0)
#         labels[-2].set_x(labels[-2].get_position()[0] - 0.02)
#         print(labels, labels[-1].get_position())
        labels = ax_ratio.set_xticklabels(['%d'%k for k in chosen_top_ks], minor=False, rotation=0)
#         labels[-1].set_x(labels[-1].get_position()[0] + 0.3)
        labels = ax_time.set_xticklabels(['%d'%k for k in chosen_top_ks], minor=False, rotation=0)
#         labels[-1].set_x(labels[-1].get_position()[0] + 0.1)
        
        ax_recall.set_ylim(20, 80)
        plt_helper.set_y_axis_ratio(ax_ratio, minratio, maxratio)
        plt_helper.set_y_axis_close(ax_time, mint, maxt, 0.1)
        
    plt_helper.plot_fig_legend(ncol=legend_col, legend_width=0.8, legend_height=0.2, columnspacing=-0.5)
    if len(distances)==1:
        plt_helper.plot_and_save('saved_figure/varying_k_%s'%(distance))
    else:
        plt_helper.plot_and_save('saved_figure/varying_k')
#     plt.figlegend(fontsize=16, bbox_to_anchor=(0.07,0.82,0.75,0.2), loc="center",
#                 mode="expand", borderaxespad=0, ncol=len(methods))
#     plt.savefig('saved_figure/varying_k.png', bbox_inches='tight', pad_inches=0.27, dpi=fig.dpi)
#     plt.savefig('saved_figure/varying_k.eps', bbox_inches='tight', pad_inches=0.27, dpi=fig.dpi)
#     plt.show()
    

if __name__ == '__main__':
    chosen_datasets = ['Sift']
    distances = ['l2', 'angle']
    methods = ['lccs', 'e2lsh', 'mplsh', 'mp_lccs', 'c2lsh']
#     methods = ['lccs', 'e2lsh', 'mplsh', 'c2lsh']
#     method_labels = [method_labels_map[method] for method in methods]
    method_colors = ['red', 'blue', 'green', 'purple', 'darkgoldenrod', 'yellow']
    method_markers = ['o', 's', 'd', '^', '*']
#     plot_method_varying_k(chosen_datasets, distances, methods)
    plot_method_varying_k(chosen_datasets, ['l2', 'angle'], methods, fig_width=6, fig_height=7.2, legend_col=2)
#     plot_method_varying_k(chosen_datasets, ['l2'], methods, fig_width=6, fig_height=8, legend_col=3)
#     plot_method_varying_k(chosen_datasets, ['angle'], methods, fig_width=6, fig_height=8, legend_col=3)