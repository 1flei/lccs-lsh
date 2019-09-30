from dataset_config import *
import os
from time import gmtime, strftime

# datasets = [Gist(), Sift(), MNIST(), Sift10M()]
datasets = [MNIST784(), Sift(), Gist()]
# datasets = [MNIST784()]
# datasets = [Sift10M()]

class AlgConfig:
    def __init__(self, **kwargs):
        self.params = kwargs

def get_dataset_path(dataset):
    return '../data/%s/%s.dsb'%(dataset.name, dataset.name)
def get_query_path(dataset):
    return '../data/%s/%s.qb'%(dataset.name, dataset.name)
def get_grount_truth_path(dataset, dist_name='l2'):
    return '../data/%s/%s.%s'%(dataset.name, dataset.name, dist_name)

def get_output_filename(dataset, method, curtime=''):
    return '../results/%s_%s_[%s].out'%(dataset.name, method, curtime)

def run_lcsb(datasets=datasets, Ls=[100], step=5, method='lcsb'):
    curtime = strftime("%m-%d_%H_%M", gmtime())
    def getArgs(ds, L, step):
        return '-A %s -n %d -q %d -d %d -D %s -Q %s -O %s -G %s -r %d -L %d --step %d --binary_input'%\
            (method, ds.n, ds.qn, ds.d, 
            get_dataset_path(ds), 
            get_query_path(ds), 
            get_output_filename(ds, method, curtime), 
            get_grount_truth_path(ds), 
            ds.r, 
            L, 
            step)
    
    for ds in datasets:
        for L in Ls:
            args = getArgs(ds, L, step)
            cmd = './lcsb '+args
            print(cmd)
            os.system(cmd)

def run_e2lsh(datasets=datasets, Ls=[100], Ks=[10]):
    method = 'e2lsh'
    curtime = strftime("%m-%d_%H_%M", gmtime())
    def getArgs(ds, L, K):
        return '-A %s -n %d -q %d -d %d -D %s -Q %s -O %s -G %s -r %d -L %d -K %d --binary_input'%\
            (method, ds.n, ds.qn, ds.d, 
            get_dataset_path(ds), 
            get_query_path(ds), 
            get_output_filename(ds, method, curtime), 
            get_grount_truth_path(ds), 
            ds.r, 
            L, 
            K)
    
    for ds in datasets:
        for L in Ls:
            for K in Ks:
                args = getArgs(ds, L, K)
                cmd = './lcsb '+args
                print(cmd)
                os.system(cmd)

def run_e2_srp(datasets=datasets, Ls=[100], Ks=[10]):
    method = 'srp_e2'
    curtime = strftime("%m-%d_%H_%M", gmtime())
    def getArgs(ds, L, K):
        return '-A %s -n %d -q %d -d %d -D %s -Q %s -O %s -G %s -L %d -K %d --binary_input'%\
            (method, ds.n, ds.qn, ds.d, 
            get_dataset_path(ds), 
            get_query_path(ds), 
            get_output_filename(ds, method, curtime), 
            get_grount_truth_path(ds, 'angle'), 
            L, 
            K)
    
    for ds in datasets:
        for L in Ls:
            for K in Ks:
                args = getArgs(ds, L, K)
                cmd = './lcsb '+args
                print(cmd)
                os.system(cmd)

def run_e2_polytope(datasets=datasets, Ls=[100], Ks=[10]):
    method = 'polytope_e2'
    curtime = strftime("%m-%d_%H_%M", gmtime())
    def getArgs(ds, L, K):
        return '-A %s -n %d -q %d -d %d -D %s -Q %s -O %s -G %s -L %d -K %d --binary_input'%\
            (method, ds.n, ds.qn, ds.d, 
            get_dataset_path(ds), 
            get_query_path(ds), 
            get_output_filename(ds, method, curtime), 
            get_grount_truth_path(ds, 'angle'), 
            L, 
            K)
    
    for ds in datasets:
        for L in Ls:
            for K in Ks:
                args = getArgs(ds, L, K)
                cmd = './lcsb '+args
                print(cmd)
                os.system(cmd)
                
def run_c2_polytope(datasets=datasets, Ls=[100], Ks=[10]):
    method = 'polytope_e2'
    curtime = strftime("%m-%d_%H_%M", gmtime())
    def getArgs(ds, L, K):
        return '-A %s -n %d -q %d -d %d -D %s -Q %s -O %s -G %s -L %d -K %d --binary_input'%\
            (method, ds.n, ds.qn, ds.d, 
            get_dataset_path(ds), 
            get_query_path(ds), 
            get_output_filename(ds, method, curtime), 
            get_grount_truth_path(ds, 'angle'), 
            L, 
            K)
    
    for ds in datasets:
        for L in Ls:
            for K in Ks:
                args = getArgs(ds, L, K)
                cmd = './lcsb '+args
                print(cmd)
                os.system(cmd)

def run_lccs_polytope(datasets=datasets, Ls=[100], step=5, method = 'polytope_lccs'):
    curtime = strftime("%m-%d_%H_%M", gmtime())
    def getArgs(ds, L, step):
        return '-A %s -n %d -q %d -d %d -D %s -Q %s -O %s -G %s -L %d --step %d --binary_input'%\
            (method, ds.n, ds.qn, ds.d, 
            get_dataset_path(ds), 
            get_query_path(ds), 
            get_output_filename(ds, method, curtime), 
            get_grount_truth_path(ds, 'angle'), 
            L, 
            step)
    
    for ds in datasets:
        for L in Ls:
            args = getArgs(ds, L, step)
            cmd = './lcsb '+args
            print(cmd)
            os.system(cmd)

def run_c2lsh(datasets=datasets, Ls=[100]):
    method = 'c2lsh'
    curtime = strftime("%m-%d_%H_%M", gmtime())
    def getArgs(ds, L):
        return '-A %s -n %d -q %d -d %d -D %s -Q %s -O %s -G %s -r %d -L %d --binary_input'%\
            (method, ds.n, ds.qn, ds.d, 
            get_dataset_path(ds), 
            get_query_path(ds), 
            get_output_filename(ds, method, curtime), 
            get_grount_truth_path(ds), 
            ds.r, 
            L)
    
    for ds in datasets:
        for L in Ls:
            args = getArgs(ds, L)
            cmd = './lcsb '+args
            print(cmd)
            os.system(cmd)

def run_srp_scan(datasets=datasets, Ls=[100]):
    method = 'srp_scan'
    curtime = strftime("%m-%d_%H_%M", gmtime())
    def getArgs(ds, L):
        return '-A %s -n %d -q %d -d %d -D %s -Q %s -O %s -G %s -L %d --binary_input'%\
            (method, ds.n, ds.qn, ds.d, 
            get_dataset_path(ds), 
            get_query_path(ds), 
            get_output_filename(ds, method, curtime), 
            get_grount_truth_path(ds, 'angle'), 
            L)
    
    for ds in datasets:
        for L in Ls:
            args = getArgs(ds, L)
            cmd = './lcsb '+args
            print(cmd)
            os.system(cmd)

def run_mplsh(datasets=datasets, Ls=[1, 2, 3, 4], Ks=list(range(4, 16)), rratio=1., method='mplsh_lshkit'):
    curtime = strftime("%m-%d_%H_%M", gmtime())
    def getArgs(ds, L, K):
        return '-A %s -n %d -q %d -d %d -D %s -Q %s -O %s -G %s -r %d -L %d -K %d --binary_input'%\
            (method, ds.n, ds.qn, ds.d, 
            get_dataset_path(ds), 
            get_query_path(ds), 
            get_output_filename(ds, method, curtime), 
            get_grount_truth_path(ds), 
            ds.r * rratio, 
            L, 
            K)
    
    for ds in datasets:
        for L in Ls:
            for K in Ks:
                args = getArgs(ds, L, K)
                cmd = './lcsb '+args
                print(cmd)
                os.system(cmd)

def run_falconn(datasets=datasets, Ls=[10, 20, 30, 40, 50], nBitss=[18], method='falconn'):
    curtime = strftime("%m-%d_%H_%M", gmtime())
    def getArgs(ds, L, nBits):
        return '-A %s -n %d -q %d -d %d -D %s -Q %s -O %s -G %s -L %d --nHashBits %d --binary_input'%\
            (method, ds.n, ds.qn, ds.d, 
            get_dataset_path(ds), 
            get_query_path(ds), 
            get_output_filename(ds, method, curtime), 
            get_grount_truth_path(ds, 'angle'), 
            L, 
            nBits)
    
    for ds in datasets:
        for L in Ls:
            for nBits in nBitss:
                args = getArgs(ds, L, nBits)
                cmd = './lcsb '+args
                print(cmd)
                os.system(cmd)

if __name__ == '__main__':
    # datasets = [Sift()]
    # datasets = [Glove()]
    datasets = [Trevi()]

    # run_lcsb(datasets=datasets, Ls=[8, 13, 21, 34, 55, 89, 144, 233], method='lcsb_reorder')
    run_lcsb(datasets=datasets, Ls=[8, 16, 32, 64, 128, 256, 512], method='lcsb', step=1)
    run_e2lsh(datasets=datasets, Ls=[300], Ks=list(range(5, 16)))
    run_mplsh(datasets=datasets, Ls=[1, 2, 3, 4], rratio=1, method='mplsh_lshkit')
    # run_e2_srp(datasets=datasets, Ls=[300], Ks=list(range(5, 16)))
    # run_c2lsh(datasets=datasets, Ls=[32, 64, 128, 256, 512])

    # run_lcsb(datasets=datasets, Ls=[32, 64, 128, 256, 512, 1024, 2048, 4096])
    run_lccs_polytope(datasets=datasets, Ls=[8, 16, 32, 64, 128, 256, 512, 1024], step=1)
    run_falconn(datasets=datasets, nBitss=[18])
    run_e2_polytope(datasets=datasets, Ls=[50, 100, 150, 200], Ks=list(range(2, 9)))