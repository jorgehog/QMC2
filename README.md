QMC2
==================
QMC2 is an efficient diffusion monte-carlo implementation coded in C++ supported by Python scripts.

To compile, go to the QMC2/qmakeQMC2 folder and run

```
qmake -makefile QMC.pro
make
```

Dependencies are

 - armadillo (http://arma.sourceforge.net)
 - MPI (can be turned off in the QMC2/src/defines.h)
 - boost (http://www.boost.org)




##VERY basic intro to how things work:


> ###Adding a new system:
-Implement your potential as a subclass of Potentials. See mainfile on how to add it to the solver.
-Implement your basis functions as BasisFunctions subclasses. Do use the OrbitalGenerator (SymPy) if things are extensive.
-Set up your basis in an Orbitals subclass, see the code for examples. This creates the mapping from single particle states to the reference Slater determinant (for Fermions).

Running your system:
-Select your orbitals
-Load either Fermions or Bosons in a System object.
-Push Potentials into your System object (add_potential(Potential* pot))
-Select a jastrow factor
-Select either BruteForce of Importance as your sampling method.
-REMEMBER TO ADD THE COULUMB INTERACTION IF YOU WANT IT (NOT DEFAULT)
-Mash these objects into a systemObject struct and feed it to a VMC/DMC constructor.
-Good to go!

For more details regarding closed form expressions for derivatives, parameterspace energy minimization etc, see the actual code for examples and comments.
#########


#########
TODO list:
#########

#Scripting#:
 -Change the input to be parsed by libconfig and not throught the CML. Goal: Deprecate runQMC.py
 -Get DCViz out of the repository and try/except everywhere it is used (remove dependency)
 -Cleanup?

#Coding#

Optimization:
-Compile different executables for 2 and three dimensions?
-have an option to fix the particle number (and initialize by arma::mat::fixed<N,N>?
-have an option to hardCode and compile a specific iniFile?
-Overall make the cache use more efficient.

Implementation:
-Slater orbitals for atoms (not sure if we need more SymPy here, as everything is given in terms of the hydrogenic ones which are already in)
See http://en.wikipedia.org/wiki/Slater-type_orbital
-Hartree-Fock for atoms. Preliminary code is implemented, however, not functioning. Need values for coulumb integrals for p orbitals and so on (higher l)
-OrbitalsWrapper for making a single system into a N-system. Exists already for N=2. Got a scheme for generalizing it, but never got to it..
-Generalize the Wrapper to different types of atoms (not n_p_local = n_p/N, but for heterogeneous distributions). This is just standard coding..

