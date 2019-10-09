from dataset_config import *
import os
import sys


def get_dataset_path(dataset, isbinary=True):
    suffix = 'dsb' if isbinary else 'ds'
    return '../data/%s/%s.%s'%(dataset.name, dataset.name, suffix)
def get_query_path(dataset, isbinary=True):
    suffix = 'qb' if isbinary else 'q'
    return '../data/%s/%s.%s'%(dataset.name, dataset.name, suffix)


def run_ground_truth(datasets, dist='l2', isbinary=True):
    dist_alg_dict = {
        'l2': 'ground_truth_l2', 
        'angle': 'ground_truth_angle', 
        'cs': 'ground_truth_cosine_similarity'
    }

    def get_ground_truth_filename(ds, dist):
        return '../data/%s/%s.%s'%(ds.name, ds.name, dist)

    def getArgs(ds):
        return '-n %s --qn %s -d %s -D %s -Q %s -O %s'%(ds.n, ds.qn, ds.d, get_dataset_path(ds, isbinary), get_query_path(ds, isbinary), get_ground_truth_filename(ds, dist))
        
    for ds in datasets:
        args = getArgs(ds)

        binary_input = '--binary_input' if isbinary else ''
        cmd = './lcsb -A %s %s %s'%(dist_alg_dict[dist], args, binary_input)
        print(cmd)
        os.system(cmd)

datasets = [Sift()]
if __name__ == '__main__':
    run_ground_truth(datasets=datasets, dist='l2')
    run_ground_truth(datasets=datasets, dist='angle')
    # run_ground_truth(datasets=datasets, dist='cs', isbinary=False)