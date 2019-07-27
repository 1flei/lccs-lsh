from dataset_config import *
from weight_gen import *
import os
import sys

datasets = [MNIST784()]
weights= ['identical', 'negative', 'uniform', 'binary', 'normal', 'onehot']

def run_ground_truth(datasets=datasets, weightsets=weights):
    def getArgs(ds, ws):
        return '-n %s --qn %s -d %s -D %s -Q %s -O %s -W %s'%(ds.n, ds.qn, ds.d, ds.ds, ds.qs, ws.gt, ws.ws)
        
    for ds in datasets:
        for w in weightsets:
            ws = WeightSet(ds, w)
            checkWeight(ds, ws)
            args = getArgs(ds, ws)
            cmd = './alsh -A ground_truth_weighted_distance '+args
            print(cmd)
            os.system(cmd)

if __name__ == '__main__':
    run_ground_truth(datasets=datasets, weightsets=weights)