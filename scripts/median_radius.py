import csv
import numpy as np


def get_filename(dataset_name):
    return '../data/%s/%s.l2'%(dataset_name, dataset_name)

def read(dataset_name):
    filename = get_filename(dataset_name)
    print(dataset_name)
    with open(filename, 'r') as f:
        reader = csv.reader(f, delimiter=' ')
        
        header = next(reader)
        print(header)
        xs = []
        for row in reader:
            # print(row)
            x = list(map(lambda x:float(x), row[:-1]))
            xs += [x]
        # print(xs)
        xs = np.array(xs)[:, 1::2]
        # print(xs)
        medians = np.median(xs, axis=0)
        print(medians[:10])


# read('Mnist')
# read('Sift')
# read('Gist')
# read('Mnist784')
# read('Trevi')
# read('NUSW')
# read('deep')
# read('glove')
# read('glove100')
read('Msong')