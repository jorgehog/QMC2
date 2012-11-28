/* 
 * File:   Blocking.cpp
 * Author: jorgmeister
 * 
 * Created on October 29, 2012, 3:58 PM
 */

#include "../../QMCheaders.h"

Blocking::Blocking(int n_c, ParParams & pp,
        std::string filename,
        std::string path,
        int n_b,
        int maxb,
        int minb,
        bool rerun)
: ErrorEstimator(n_c, filename, path, pp.parallel, pp.node, pp.n_nodes, rerun) {
    n_block_samples = n_b;
    max_block_size = maxb;
    min_block_size = minb;


    if (rerun) {
        if (minb < n_nodes) {
            std::cout << "minimum block size must at least match the number of nodes used." << std::endl;
            exit(1);
        }

        node_comm_scatter_data();

        init_file();

    }


}

Blocking::Blocking(int n_c, std::string filename, std::string path, int n_b, int maxb, int minb)
: ErrorEstimator(n_c, filename, path, false, 0, 1, false) {
    n_block_samples = n_b;
    max_block_size = maxb;
    min_block_size = minb;

}

double Blocking::estimate_error() {

    if (!rerun) {
        return 0;
    }

    int block_size, block_step_length;
    double error, var, mean;

    block_step_length = (max_block_size - min_block_size) / n_block_samples;


    get_initial_error();

    for (int j = 0; j < n_block_samples; j++) {
        block_size = (min_block_size + j * block_step_length) / n_nodes;

        block_data(block_size, var, mean);

        var = combine_variance(var, mean);

        if (is_master) {
            error = sqrt(var / (((n_nodes * n_c) / block_size) - 1.0));
            file << block_size << "\t" << error << std::endl;
        }
    }

    return error;
}

void Blocking::block_data(int block_size, double &var, double &mean) {
    using namespace arma;

    double block_mean;


    mean = 0;
    double mean2 = 0;

    int n_b = n_c / block_size;
    for (int j = 0; j < n_b; j++) {
        block_mean = sum(data(span(j*block_size, (j + 1) * block_size - 1))) / block_size;

        mean += block_mean;
        mean2 += block_mean*block_mean;
    }

    mean2 /= (n_b);
    mean /= (n_b);

    var = mean2 - mean*mean;
}

void Blocking::get_initial_error() {
    double var = arma::var(data);
    double mean = arma::mean(data);
    var = combine_variance(var, mean);
    if (is_master) std::cout << "Initial stddev: " << sqrt(var / (n_nodes * n_c - 1)) << std::endl;
}