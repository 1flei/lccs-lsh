import csv
import struct
import numpy as np

def to_binary(dataset_name):
    inDataFilename = '../data/%s/%s.ds'%(dataset_name, dataset_name)
    inQueryFilename = '../data/%s/%s.q'%(dataset_name, dataset_name)
    outDataFilename = '../data/%s/%s.dsb'%(dataset_name, dataset_name)
    outQueryFilename = '../data/%s/%s.qb'%(dataset_name, dataset_name)

    print(dataset_name)
    with open(inDataFilename, 'r') as fin, open(outDataFilename, 'wb') as fout:
        reader = csv.reader(fin, delimiter=' ')
        
        for row in reader:
            # print(row)
            x = list(map(lambda x:float(x), row[1:]))

            fout.write(struct.pack('f'*len(x), *x))          
    with open(inQueryFilename, 'r') as fin, open(outQueryFilename, 'wb') as fout:
        reader = csv.reader(fin, delimiter=' ')
        
        for row in reader:
            # print(row)
            x = list(map(lambda x:float(x), row[1:]))

            fout.write(struct.pack('f'*len(x), *x))           


# to_binary('Mnist')
# to_binary('Sift')
# to_binary('Gist')
to_binary('Mnist784')