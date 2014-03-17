#include <QMC2.h>

#include "sampleforce.h"

#include <armadillo>
#include <iomanip>
#include <vector>

using namespace QMC2;

void forceLoop(ASGD & asgd, VMC & vmc, DMC & dmc, int n_p, bool is_master);

int main(int argc, char** argv) {


    using namespace std;

    //Initialize structs for holding constructor arguments.
    struct ParParams parParams;
    struct VMCparams vmcParams;
    struct DMCparams dmcParams;
    struct VariationalParams variationalParams;
    struct GeneralParams generalParams;
    struct MinimizerParams minimizerParams;
    struct SystemObjects systemObjects;

    if (argc == 2) {
        generalParams.runpath = argv[1];
    }

    minimizerParams.SGDsamples = 400;

    dmcParams.therm=1000;
    dmcParams.n_c = 2000;
    vmcParams.n_c = 1E6;
    dmcParams.dt = 0.0001;

//    variationalParams.alpha = 1;//0.922925;
    variationalParams.beta  =  0.6106;

    //Setting up parallel parameters
    initMPI(parParams, argc, argv);

    scaleWithProcs(parParams, generalParams, minimizerParams, vmcParams, dmcParams);

    //Setting up the solver parameters
    generalParams.n_p = 16;
    generalParams.dim = 3;

    double R = 2.282;
    BodyDef Oxygen1;
    Oxygen1.n_p_local = 8;
    Oxygen1.origin = {R/2, 0,  0};

    BodyDef Oxygen2;
    Oxygen2.n_p_local = 8;
    Oxygen2.origin = {-R/2,  0,  0};

    std::vector<BodyDef> bodies;
    //    bodies.push_back(Silicon);
    bodies.push_back(Oxygen1);
    bodies.push_back(Oxygen2);


    arma::mat C(generalParams.n_p, generalParams.n_p);
    C.eye();

    OrbitalsFactory expBasisFactory(ATOMS);
    OrbitalsFactory factory(EXPANDED);
    factory.basisForExpanded = &expBasisFactory;
    factory.C = C.t();

    NBodyTransform* Oxygen = new NBodyTransform(generalParams, variationalParams, bodies, factory);
    systemObjects.SP_basis = Oxygen;
    Fermions system(generalParams, Oxygen);
//    system.add_potential(new MolecularCoulomb(generalParams, Oxygen));
    system.add_potential(new Coulomb(generalParams));
    systemObjects.system = &system;

    systemObjects.jastrow = new Pade_Jastrow(generalParams, variationalParams);
    //    systemObjects.jastrow = new No_Jastrow();
    systemObjects.sample_method = new Importance(generalParams);

    if (parParams.is_master) std::cout << "seed: " << generalParams.random_seed << std::endl;

    minimizerParams.alpha = {};
    minimizerParams.beta = {0.5};

    //Creating the solver objects
    bool silent = true;

    VMC vmc(generalParams, vmcParams, systemObjects, parParams, dmcParams.n_w, silent);
    ASGD asgd(&vmc, minimizerParams, parParams, "/tmp/");
    DMC dmc(generalParams, dmcParams, systemObjects, parParams, silent);
        
    vmc.set_error_estimator(new Blocking(vmcParams.n_c, parParams));

    arma::wall_clock a;
    a.tic();
    asgd.minimize();
    vmc.run_method();
    vmc.dump_subsamples();
    dmc.initFromVMC(&vmc);
    dmc.run_method();

    if (parParams.is_master) std::cout << "Time: " <<setprecision(3) << fixed << a.toc()/60 << std::endl;


    if (parParams.is_master) cout << "~.* QMC fin *.~" << endl;

#ifdef MPI_ON
    MPI_Finalize();
#endif

    return 0;
}
