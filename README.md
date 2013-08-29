QMC2
==================
QMC2 is an efficient diffusion monte-carlo implementation coded in C++ supported by Python scripts.

-----------------------

###Compiling

To compile, go to the QMC2/qmakeQMC2 folder and run

```
qmake -makefile QMC.pro
make
```

External library dependencies:

 - armadillo (http://arma.sourceforge.net)
 - boost (http://www.boost.org)
 - MPI (can be turned off in the QMC2/src/defines.h)

----------------------------

###Adding a new system:
> - Implement your potential as a subclass of Potentials. See mainfile on how to add it to the solver.
> - Implement your basis functions as BasisFunctions subclasses. Do use the OrbitalGenerator (SymPy) if things are extensive.
> - Set up your basis in an Orbitals subclass, see the code for examples. This creates the mapping from single particle states to the reference Slater determinant (for Fermions).

###Seting up mainfile:

Here follows an example of a mainfile structure.

```
//...

//Initializing the system

struct GeneralParams gP;
struct SystemObjects sO;
struct VariationalParams vP;

gP.dim = 3;
gP.n_p = 2;

sO.SP_basis = new DiTransform(gP, vP);

System* system = new Fermions(gP, sO.SP_basis);
system->add_potential(new DiAtomCore(gP));
system->add_potential(new Coulomb(gP));

sO.SYSTEM = system;

sO.sample_method = new Importance(gP);
sO.jastrow = new Pade_Jastrow(gP, vP);

//Initializing the solver

struct DMCparams dmcParams;
dmcParams.dt = 0.001;
dmcParams.n_b = 100;
dmcParams.n_c = 1000;
dmcParams.n_w = 1E5;
dmcParams.therm = 1000;

struct VMCparams vmcParams;
vmcParams.dt = 0.01;
vmcParams.n_c = 1E8;

vP.alpha = 1.0;
vP.beta = 0.5;

//initialize parallelization (or not!)

struct ParParams parParams;

parParams.parallel = false;
parParams.node = 0;
parParams.n_nodes = 1;
parParams.is_master = true;


//create a VMC solver object initialized to prepare DMC walkers.
//adding an error estimator which computes the variance of the samples.
VMC vmc(gP, vmcParams, sO, parParams, dmcParams.n_w);
vmc->set_error_estimator(new SimpleVar(parParams));
vmc->run_method();

//and the same for DMC. The last argument is the VMC-object which now have prepped walkers ready.
DMC dmc(gP, dmcParams, sO, parParams, &vmc);
dmc->set_error_estimator(new SimpleVar(parParams));
dmc->run_method();

//...

```

For more details regarding closed form expressions for derivatives, parameterspace energy minimization etc, see the actual code for examples and comments.



TODO list:
-------------

###Scripting:
> - Change the input to be parsed by libconfig and not throught the CML. Goal: Deprecate runQMC.py
> - Get DCViz out of the repository and try/except everywhere it is used (remove dependency)
> - Cleanup?


###Optimization:
> - Compile different executables for 2 and three dimensions?
> - have an option to fix the particle number (and initialize by arma::mat::fixed<N,N>?
> - have an option to hardCode and compile a specific iniFile?
> - Overall make the cache use more efficient.

###Implementation:
> - Slater orbitals for atoms (not sure if we need more SymPy here, as everything is given in terms of the hydrogenic ones which are already in)
>See http://en.wikipedia.org/wiki/Slater-type_orbital
> - Hartree-Fock for atoms. Preliminary code is implemented, however, not functioning. Need values for coulumb integrals for p orbitals and so on (higher l)
> - OrbitalsWrapper for making a single system into a N-system. Exists already for N=2. Got a scheme for generalizing it, but never got to it..
> - Generalize the Wrapper to different types of atoms (not n_p_local = n_p/N, but for heterogeneous distributions). This is just standard coding..

