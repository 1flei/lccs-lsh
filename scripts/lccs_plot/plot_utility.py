import numpy as np
import matplotlib.ticker as mticker
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

# def plot_fig_legend(plt, legend_width=0.75, legend_height=0.2, ncol=5):
#     bbox = (0.5-legend_width/2., 0.98-legend_height, legend_width, legend_height)
#     plt.figlegend(fontsize=16, loc='upper center', bbox_to_anchor=bbox, 
#                 mode="expand", borderaxespad=0.05, ncol=ncol)
# 
#     
# def plot_subplots_adjust(plt, fig_width, fig_height, left_space=0.8, bottom_space=0.55, top_space=0.8, right_space=0.25):
#     fig = plt.figure(figsize=(fig_width, fig_height))
#     plt.rcParams.update({'font.size': 13})
#     
#     bottom = bottom_space / fig_height
#     top = (fig_height-top_space) / fig_height
#     left = left_space / fig_width
#     right = (fig_width - right_space) / fig_width 
#     plt.subplots_adjust(bottom=bottom, top=top, left=left, right=right, wspace=0.22)
    
    
    
class PlotHelper:
    def __init__(self, plt, fig_width, fig_height):
        self.plt = plt
        self.fig_width  = fig_width
        self.fig_height = fig_height
        self.minx = None
        self.maxx = None
        
    def plot_fig_legend(self, legend_width=0.9, legend_height=0.2, ncol=5, columnspacing=None):
        bbox = (0.5-legend_width/2., 0.98-legend_height, legend_width, legend_height)
        
#         handles, labels = self.plt.get_legend_handles_labels()
        self.plt.figlegend(fontsize=16, loc='upper center', bbox_to_anchor=bbox, 
                    mode="expand", borderaxespad=0.05, ncol=ncol, columnspacing=columnspacing)
    
        
    def plot_subplots_adjust(self, left_space=0.8, bottom_space=0.55, top_space=0.9, right_space=0.25, wspace=0.24, hspace=0.3):
        fig = self.plt.figure(figsize=(self.fig_width, self.fig_height))
        self.plt.rcParams.update({'font.size': 13})
        
        bottom = bottom_space / self.fig_height
        top = (self.fig_height-top_space) / self.fig_height
        left = left_space / self.fig_width
        right = (self.fig_width - right_space) / self.fig_width 
        self.plt.subplots_adjust(bottom=bottom, top=top, left=left, right=right, wspace=wspace, hspace=hspace)
        
#     def semilogy(self, ax, x, y, *args, **kwargs):
#         ax.semilogy(x, y, *args, **kwargs)
#         self.minx = np.min(x)
#         self.maxx = np.max(x)

    def get_ticks(self, miny, maxy, possible_ticks=[1, 10, 100], error=0):
        for ipt, pt in enumerate(possible_ticks[::-1]):
            if pt <= miny-error:
                pt_lower = len(possible_ticks)-ipt-1
                break
            
        for ipt, pt in enumerate(possible_ticks):
            if pt >= maxy+error:
                pt_upper = ipt
                break
#         print(possible_ticks[pt_lower], possible_ticks[pt_upper])
        return possible_ticks[pt_lower:pt_upper+1]

    def set_y_axis_ratio(self, ax, miny, maxy):
        possible_ticks = [1., 1.02, 1.04, 1.06, 1.08]
        ticks = self.get_ticks(miny, maxy, possible_ticks)
        
        #disable minor ticks 
        ax.set_yticks([], minor=True)
#         f = mticker.ScalarFormatter(useOffset=False, useMathText=True)
#         ticks_labels = ["${}$".format(f._formatSciNotation('%1.10e' % t)) for t in ticks]
        ax.set_yticks(ticks)
#         ax.set_yticklabels(ticks_labels)
        ax.set_ylim(ticks[0], ticks[-1])
            
#         y0 = 10**np.floor(np.log(miny)/np.log(10))
#         y1 = 10**np.ceil(np.log(maxy)/np.log(10))
#         print('y0=', y0, 'y1=', y1)
#         ax.set_ylim((y0, y1))

    def set_y_axis_log10(self, ax, miny, maxy):
        possible_ticks = [10**(i-10) for i in range(20)]
        ticks = self.get_ticks(miny, maxy, possible_ticks)
        
        #disable minor ticks 
        ax.set_yticks([], minor=True)
#         f = mticker.ScalarFormatter(useOffset=False, useMathText=True)
#         ticks_labels = ["${}$".format(f._formatSciNotation('%1.10e' % t)) for t in ticks]
        ax.set_yticks(ticks)
#         ax.set_yticklabels(ticks_labels)
        ax.set_ylim(ticks[0], ticks[-1])
            
#         y0 = 10**np.floor(np.log(miny)/np.log(10))
#         y1 = 10**np.ceil(np.log(maxy)/np.log(10))
#         print('y0=', y0, 'y1=', y1)
#         ax.set_ylim((y0, y1))

    def set_y_axis_close(self, ax, miny, maxy, error=0.):
        possible_ticks = []
        for i in range(20):
            possible_ticks += [1*10**(i-10), 2*10**(i-10), 5*10**(i-10)]
        ticks = self.get_ticks(miny, maxy, possible_ticks, error)
        
        #disable minor ticks 
        ax.set_yticks([], minor=True)
        ticks_labels = ['%d'%t for t in ticks]
        ax.set_yticks(ticks)
        ax.set_yticklabels(ticks_labels)
        ax.set_ylim(ticks[0], ticks[-1])
#         ax.set_ticks()
            
        
#     def set_y_axis_close(self, ax, miny, maxy):
#         logminy = np.log(miny)/np.log(10)
#         logmaxy = np.log(maxy)/np.log(10)
#         
#         y0 = 10**np.floor(logminy)
#         y1 = 10**np.ceil(logmaxy)
#         y0_10 = y0
#         y1_10 = y1/10
#         
#         miny_ = max(y0, np.floor(miny / y0_10) * y0_10)
# #         maxy_ = min(y1, np.ceil(maxy / y1_10) * y1_10)
#         maxy_ = y1
#         
#         print(miny, maxy, miny_, maxy_)
#         
#         y_ticks = []
#         minor_ticks = []
#         if miny_ != y0:
#             y_ticks += [miny_]
#         else:
#             y_ticks += [y0]
#         for i in range(0, 1000):
#             if y0*10**i < maxy_:
#                 if y0*10**i > miny:
#                     y_ticks += [y0*10**i]
#             else:
#                 break
#             
#             if y0*10**i * 2 < maxy and y0*10**i * 2 > miny:
#                 minor_ticks += [y0*10**i * 2]
#             if y0*10**i * 5 < maxy and y0*10**i * 5 > miny:
#                 minor_ticks += [y0*10**i * 5]
# #         if maxy_ != y1:
# #             y_ticks += [maxy_]
#         y_ticks += [maxy_]
#             
#         print(y_ticks, ['%d'%(int(y_tick), ) for y_tick in y_ticks])
#         
# #         yticks = [10**i for i in range(int(np.log(y0)/np.log(10)), int(np.log(y1)/np.log(10))+1)]
# #         print(yticks)
# #         ax.set_yticks(y_ticks, ['%d'%(int(y_tick), ) for y_tick in y_ticks
#         ax.set_ylim((miny_, maxy_))
#         print('existing y_ticks', ax.get_yticks())
#         ax.set_yticks(y_ticks)
#         ax.set_yticklabels(['%d'%(int(y_tick), ) for y_tick in y_ticks])
#         print('existing y_ticks', ax.get_yticks())
# #         minor_ticks = ax.get_yticks(minor=True)
#         
#         ax.set_yticks(minor_ticks, minor=True)
#         ax.set_yticklabels(['%d'%(int(y_tick), ) for y_tick in minor_ticks], minor=True)
#         ax.set_ylim((y0, y1))
        
    def plot_and_save(self, filename):
        self.plt.savefig('%s.png'%filename)
        self.plt.savefig('%s.eps'%filename)
        self.plt.show()