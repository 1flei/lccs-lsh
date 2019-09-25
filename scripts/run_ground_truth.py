from dataset_config import *
import os
import sys

datasets = [Glove()]

def get_dataset_path(dataset):
    return '../data/%s/%s.dsb'%(dataset.name, dataset.name)
def get_query_path(dataset):
    return '../data/%s/%s.qb'%(dataset.name, dataset.name)


def run_ground_truth(datasets=datasets, dist='l2'):
    dist_alg_dict = {
        'l2': 'ground_truth_l2', 
        'angle': 'ground_truth_angle', 
        'cs': 'ground_truth_cosine_similarity'
    }

    def get_ground_truth_filename(ds, dist):
        return '../data/%s/%s.%s'%(ds.name, ds.name, dist)

    def getArgs(ds):
        return '-n %s --qn %s -d %s -D %s -Q %s -O %s'%(ds.n, ds.qn, ds.d, get_dataset_path(ds), get_query_path(ds), get_ground_truth_filename(ds, dist))
        
    for ds in datasets:
        args = getArgs(ds)
        cmd = './lcsb -A %s %s --binary_input'%(dist_alg_dict[dist], args)
        print(cmd)
        os.system(cmd)

if __name__ == '__main__':
    run_ground_truth(datasets=datasets, dist='l2')
    run_ground_truth(datasets=datasets, dist='angle')
    # run_ground_truth(datasets=datasets, dist='cs')