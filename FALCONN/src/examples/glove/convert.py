#!/usr/bin/python

import sys
import struct
import numpy as np
import re

pat = re.compile(r'-?\d+(>?\.\d+)?')

def isfloat(value):
  try:
    float(value)
    return True
  except ValueError:
    return False

matrix = []
with open('dataset/glove.840B.300d.txt', 'r') as inf:
    with open('dataset/glove.840B.300d.dat', 'wb') as ouf:
        counter = 0
        for line in inf:
            try:
                row = [float(x) for x in line.split()[1:]]
            except:
                print('exceptation:  at line', counter)
                print(line)
                row = [float(x) for x in line.split()[1:] if isfloat(x) ]
            
            if len(row) != 300:
                row = [float(x) for x in line.split()]
                if len(row)!=300:
                    print('mismatch!!!! at line', counter)
                    print(line)
                    print('row=', row)
            assert len(row) == 300
            ouf.write(struct.pack('i', len(row)))
            ouf.write(struct.pack('%sf' % len(row), *row))
            counter += 1
            matrix.append(np.array(row, dtype=np.float32))
            if counter % 10000 == 0:
                sys.stdout.write('%d points processed...\n' % counter)
np.save('dataset/glove.840B.300d', np.array(matrix))
