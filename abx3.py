"""
In addition to default keywords, two keywords are defined
process -> unique for every calculation
state
"""
import numpy
from ase import Atoms
from ase.optimize.bfgs import BFGS
from ase.constraints import UnitCellFilter, StrainFilter
from gpaw import GPAW, PW, restart, KohnShamConvergenceError
import os, errno
from ase.db import connect

class ABX3():
    """
    A-B-X3 cubic perovskites
    creates an instance of "ase Atoms object" 
    accessible via abx3.atoms from elements
    names A, B, and X from a saved state or 
    a given lattice constant(s)
    """
    def __init__(self, A, B, X, state, process='any', prefix='~/abx3_cached/'):
        # handling files
        self.prefix = os.path.expanduser( prefix )
        self.path =  self.prefix + A + B + X +'3/'
        try:
            os.makedirs(self.path)
        except OSError as err:
            if err.errno != errno.EEXIST:
                raise
        # unique process
        self.A, self.B, self.X = A, B, X
        self.formula = A + B + X + '3'
        self.process = '_'.join( [self.formula, process] )
        self.db = connect( self.prefix + 'database.db' )
        self.pid = self.db.reserve( process = self.process )
        # make atoms
        if self.pid is None:
            self.atoms = None
        else:
            self.make_atoms( A, B, X, state )

    def make_atoms(self, A, B, X, state):
        if isinstance(state,str):
            self.atoms = self.db.get_atoms( A_ion=A, B_ion=B, X_ion=X, state=state, attach_calculator=True )
        else:
            if isinstance(state,float):
                cell = numpy.diag( 3*[state] )
            elif isinstance(state,numpy.ndarray):
                if state.shape==(3,):
                    cell = numpy.diag( state )
                if state.shape==(3,3):
                    cell = state
            positions = [ (0,0,0), 
                          (cell[0]+cell[1]+cell[2])/2,
                          (cell[0]+cell[1])/2,
                          (cell[0]+cell[2])/2,
                          (cell[1]+cell[2])/2
                        ]
            self.atoms = Atoms( A+B+3*X, positions=positions, cell=cell, pbc=[1,1,1] )
        # if setups not found
        spath = os.environ['GPAW_SETUP_PATH'] + '/'
        for elm in [A,B,X]:
            if not ( os.path.exists(spath+elm+'.LDA.gz') and 
                     os.path.exists(spath+elm+'.RPBE.gz') ):
                self.savestate('nosetup')
                break

    def savestate(self,state):
        if self.atoms is not None:
            sid = self.db.write( self.atoms, process=self.process, A_ion=self.A, B_ion=self.B, X_ion=self.X, state=state )
            self.end_process()

    def __del__(self):
        self.savestate('interrupted')

    def end_process(self):
        del self.db[self.pid]
        self.pid = None
        self.atoms = None


    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    def done(self,stage):
        if self.atoms is not None:
            try: 
                self.atoms, _ = restart( self.path + "abx3_"+stage+".gpw" )
                return True
            except FileNotFoundError:
                return False

    def record(self,stage):
        if self.atoms is not None:
            self.atoms.calc.write( self.path + "abx3_"+stage+".gpw", 'all' )

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    def attach_accurate_calc(self):
        if self.atoms is not None:
            self.atoms.calc = GPAW(
                               xc='RPBE',
                               mode=PW(500, dedecut='estimate'), # default is ecut = 340
                               kpts=(8, 8, 8),
                               convergence={'eigenstates': 1.e-10},
                               )

    def attach_fast_calc(self):
        if self.atoms is not None:
            self.atoms.calc = GPAW(
                               xc='LDA',
                               mode=PW(),
                               kpts=(4, 4, 4),
                               )

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    def relax(self, fmax=0.005, Filter=None, stage="Relaxation", maxsteps=100):
        if self.atoms is not None:
            if not self.done( stage ):
                self.atoms.calc.set( txt = self.path+'relaxation.txt' )
                if Filter:
                    relaxation = BFGS( Filter( self.atoms ) )
                else:
                    relaxation = BFGS( self.atoms )
                try:
                    relaxation.run( fmax=fmax, steps=maxsteps )
                except KohnShamConvergenceError:
                    self.savestate('notconverged')
                if relaxation.get_number_of_steps()==maxsteps: self.savestate('maxsteps')
                self.record( stage )

    def step1(self):
        if self.atoms is not None:
            #self.attach_fast_calc()
            #self.relax( fmax=0.05, Filter=StrainFilter, stage='quickRelaxation' )
            self.attach_accurate_calc()
            self.relax( Filter=UnitCellFilter, stage="ucRelaxation" )


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if __name__=='__main__':
    #
    abx3 = ABX3('Sn','Ti','O', 3.0, process='create')
    abx3.savestate('ini')
    #
    abx3 = ABX3('Sn','Ti','O', 'ini',process='relaxing')
    abx3.step1()
    abx3.savestate("relaxed")

