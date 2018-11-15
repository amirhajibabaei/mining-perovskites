import pandas
from ase.db import connect
from abx3 import ABX3
import sys
import time
import numpy

csv_data = pandas.read_csv('ABX3cubDBnewbigRtol.csv')
def get_csv(J):
    A, B, X = csv_data['A'][J], csv_data['B'][J], csv_data['X'][J]
    a = csv_data['RA'][J] + csv_data['RC'][J]
    return A, B, X, a

try:
    k0 = int(sys.argv[1])
except IndexError:
    k0=0

for k in range( k0, csv_data.shape[0] ): 
    A, B, X, a = get_csv(k)
    abx3 = ABX3( A, B, X, a, process='creation')
    time.sleep(0.01)
    abx3.savestate('initial')

