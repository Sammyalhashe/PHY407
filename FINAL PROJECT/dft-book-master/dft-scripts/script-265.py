from jasp import *
JASPRC['queue.nodes']=4
from ase import Atom, Atoms
atoms = Atoms([Atom('O',[5, 5, 5], magmom=1)],
             cell=(6, 6, 6))
with jasp('molecules/O_s-4nodes',
          encut=300,
          xc='PBE',
          ispin=2,
          ismear=0,
          sigma=0.001,
          setups={'O':'_s'},  # specifies O_s potential
          atoms=atoms) as calc:
    calc.calculate()