#pragma once

#include "structs.h"

#include "Sampler/Sampler.h"

namespace QMC2
{


inline void initMPI(struct ParParams & parParams, int argc, char ** argv){
#ifdef MPI_ON
    int node, n_nodes;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &n_nodes);
    MPI_Comm_rank(MPI_COMM_WORLD, &node);

    parParams.n_nodes = n_nodes;
    parParams.node = node;
    parParams.parallel = (parParams.n_nodes > 1);
    parParams.is_master = (parParams.node == 0);

    int nodeSum = 0;
    MPI_Reduce(&node, &nodeSum, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    if (parParams.is_master)
    {
        if (nodeSum != ((n_nodes - 1)*n_nodes)/2)
        {
            std::cerr << "Error in MPI init: " << n_nodes <<  " ranks not properly distributed." << std::endl;
            MPI_Finalize();
            exit(1);
        }

        else
        {
            std::cout << "MPI load successful." << std::endl;
        }
    }

#else
    parParams.parallel = false;
    parParams.node = 0;
    parParams.n_nodes = 1;
    parParams.is_master = true;
#endif

    Sampler::setParParams(parParams);

}

inline void scaleWithProcs(struct ParParams & parParams,
                           struct GeneralParams & generalParams,
                           struct MinimizerParams & minimizerParams,
                           struct VMCparams & vmcParams,
                           struct DMCparams & dmcParams){

    generalParams.random_seed -= parParams.node;
    minimizerParams.n_c_SGD /= parParams.n_nodes;

    vmcParams.n_c /= parParams.n_nodes;
    dmcParams.n_w /= parParams.n_nodes;

    if (minimizerParams.n_c_SGD == 0) {
        minimizerParams.n_c_SGD = 1;
    }

    if (vmcParams.n_c == 0) {
        vmcParams.n_c = 1;
    }

    if (dmcParams.n_w == 0) {
        dmcParams.n_w = 1;
    }

}

}
