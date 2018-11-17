import pandas
from abx3 import ABX3
import sys
import numpy
import os

vlad1 = 'ABX3cubDBnewbigRtol.csv'
vlad2 = 'ABX3cubDBGPAWRtol.csv'
csv_data = pandas.read_csv( vlad2 )

def get_csv(J):
    A, B, X = csv_data['A'][J], csv_data['B'][J], csv_data['X'][J]
    a = csv_data['RA'][J] + csv_data['RC'][J]
    return A, B, X, a

def relaxations():
    for k in range( csv_data.shape[0] ): 
        A, B, X, a = get_csv(k)
        abx3 = ABX3( A, B, X, a, process='creation')
        if abx3.atoms is not None: 
            print( "Found new job: {} indexed {} pid {}".format( A+B+X+'3', k, os.getpid() ) )
        abx3.step1()
        abx3.savestate('relaxed')

if __name__=='__main__':
    relaxations()
