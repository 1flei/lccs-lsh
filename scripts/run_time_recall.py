from dataset_config import *
from method_config import *
import os
from time import gmtime, strftime

#a decorator
def list_expansion(f):
    def g(*args, **kwargs):
        args_copy = list(args)
        for i, arg in enumerate(args):
            if type(arg) == list:
                for element in arg:
                    args_copy[i] = element
                    g(*args_copy, **kwargs)
                return 
        f(*args)
    return g


def get_dataset_path(dataset):
    return '../data/%s/%s.dsb'%(dataset.name, dataset.name)
def get_query_path(dataset):
    return '../data/%s/%s.qb'%(dataset.name, dataset.name)
def get_grount_truth_path(dataset, dist_name='l2'):
    return '../data/%s/%s.%s'%(dataset.name, dataset.name, dist_name)

def get_output_filename(dataset, method, curtime=''):
    return '../results/%s_%s_[%s].out'%(dataset.name, method, curtime)

def get_dataset_params(ds):
    return '-n %d -q %d -d %d -D %s -Q %s'%(ds.n, ds.qn, ds.d, 
        get_dataset_path(ds), 
        get_query_path(ds))

def get_dataset_distance_params(ds, distance):
    if distance=='l2':
        return '-r %f -G %s'%(ds.r, get_grount_truth_path(ds, distance))
    elif distance == 'angle':
        return '--normalized -G %s'%(get_grount_truth_path(ds, distance))
    return ''

def get_common_params(ds, method, maxqn=100, curtime=None):
    if curtime is None:
        curtime = strftime("%m-%d_%H_%M", gmtime())
    actual_qn = min(maxqn, ds.qn)
    return '-n %d -q %d -d %d -D %s -Q %s -O %s --binary_input'%\
        (ds.n, actual_qn, ds.d, 
        get_dataset_path(ds), 
        get_query_path(ds), 
        get_output_filename(ds, method, curtime))



def get_method_name(method, distance):
    method_name_dict = {
        ('lccs',    'l2'):    'lccs', 
        ('lccs_mp', 'l2'):    'mp_lccs',   
        ('e2lsh',   'l2'):    'e2lsh',
        ('mplsh',   'l2'):    'mplsh_lshkit',
        ('c2lsh',   'l2'):    'c2lsh',
        ('lccs',    'angle'): 'polytope_lccs',
        ('lccs_mp', 'angle'): 'polytope_mplccs',
        ('e2lsh',   'angle'): 'polytope_e2',
        ('mplsh',   'angle'): 'falconn',
        ('c2lsh',   'angle'): 'polytope_c2',

        ('srs',   'l2'): 'srs',
        ('qalsh',   'l2'): 'qalsh',
    }
    return method_name_dict[(method, distance)]

@list_expansion
def run_alg(method_obj, dataset, distance, curtime=None, binary_name='./lccs'):
    method = method_obj.name
    # print(method, dataset, distance)
    method_name = get_method_name(method, distance)
    method_name_arg = '-A %s'%method_name
    common_args = get_common_params(dataset, method_name, curtime=curtime)

    distance_dataset_arg = get_dataset_distance_params(dataset, distance)

    for method_arg in method_obj.for_param(distance):
        cmd = '%s %s %s %s %s'%(binary_name, method_name_arg, common_args, distance_dataset_arg, method_arg)
        print(cmd)
        os.system(cmd)


possible_algs = [LCCS(), LCCS_MP(), C2LSH(), E2LSH(), MPLSH()]
possible_distances = ['l2', 'angle']
# possible_datasets = [MNIST784(), Sift(), Sift10M(), Gist(), Trevi(), Glove(), Glove100(), Msong(), Deep()]
possible_datasets = [Sift(), Gist(), Glove100(), Msong(), Deep()]

if __name__ == '__main__':
    run_alg([LCCS(), MPLSH(), LCCS_MP(), E2LSH(), C2LSH(), SRS(), QALSH()], [Sift(), Gist(), Glove100(), Msong(), Deep()], 'l2')
    run_alg([LCCS(), MPLSH(), LCCS_MP(), E2LSH(), C2LSH()], [Sift(), Gist(), Glove100(), Msong(), Deep()], 'angle')