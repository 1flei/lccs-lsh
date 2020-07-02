from dataset_config import *
import os
import sys
from run_ground_truth import run_ground_truth
from run_time_recall import run_alg


if __name__ == '__main__':
    datasets = [MNIST()]
    run_ground_truth(datasets=datasets, dist='l2')
    run_alg([LCCS()], [MNIST()], 'l2')