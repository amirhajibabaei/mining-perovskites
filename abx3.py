from ase import Atoms
import numpy
from ase.optimize.bfgs import BFGS
from ase.constraints import UnitCellFilter, StrainFilter
from gpaw import GPAW, PW, restart
import os, errno

class ABX3():
    """
    A-B-X3 cubic perovskites
    creates instance of "ase Atoms object" from 
    A:
    B: 
    X:
    alatt: lattice constant
    """
    def __init__(self, A, B, X, alatt, prefix='~/abx3_cached/' ):
        cell = numpy.identity(3,dtype=float) * alatt
        formula = A + B + 3*X
        positions = [ (0,0,0), 
                      (cell[0]+cell[1]+cell[2])/2,
                      (cell[0]+cell[1])/2,
                      (cell[0]+cell[2])/2,
                      (cell[1]+cell[2])/2
                    ]
        self.atoms = Atoms( formula, positions=positions, cell=cell, pbc=[1,1,1] )
        # handling files
        self.path = os.path.expanduser(prefix) + A+B+X+'3/'
        try:
            os.makedirs(self.path)
        except OSError as err:
            if err.errno != errno.EEXIST:
                raise

    def get_atoms(self):
        return self.atoms

    def done(self,stage):
        try: 
            self.atoms, _ = restart( self.path + "abx3_"+stage+".gpw" )
            return True
        except FileNotFoundError:
            return False

    def record(self,stage):
        self.atoms.calc.write( self.path + "abx3_"+stage+".gpw", 'all' )

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    def attach_accurate_calc(self):
        self.atoms.calc = GPAW(
                               xc='RPBE',
                               mode=PW(500, dedecut='estimate'), # default is ecut = 340
                               kpts=(8, 8, 8),
                               convergence={'eigenstates': 1.e-10},
                               )

    def attach_fast_calc(self):
        self.atoms.calc = GPAW(
                               xc='LDA',
                               mode=PW(),
                               kpts=(4, 4, 4),
                               )

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    def relax(self, fmax=0.005, Filter=None, stage="Relaxation"):
        if not self.done( stage ):
            self.atoms.calc.set( txt = self.path+'relaxation.txt' )
            if Filter:
                relaxation = BFGS( Filter( self.atoms ) )
            else:
                relaxation = BFGS( self.atoms )
            relaxation.run( fmax=fmax )
            self.record( stage )

    def step1(self):
        self.attach_fast_calc()
        self.relax( fmax=0.05, Filter=StrainFilter, stage='quickRelaxation' )
        self.attach_accurate_calc()
        self.relax( Filter=UnitCellFilter, stage="ucRelaxation" )


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if __name__=='__main__':
    abx3 = ABX3('Sn','Ti','O', 3.0)
    abx3.step1()
