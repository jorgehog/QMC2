/* 
 * File:   Blocking.cpp
 * Author: jorgmeister
 * 
 * Created on October 29, 2012, 3:58 PM
 */

#include "Blocking.h"

#include "../../structs.h"
#include <iomanip>


Blocking::Blocking(int n_c, ParParams & pp,
        std::string filename,
        std::string path,
        int n_b,
        int maxb,
        int minb,
        bool rerun)
: ErrorEstimator(n_c, filename, path, pp.parallel, pp.node, pp.n_nodes, rerun) {
    n_block_samples = n_b;
    max_block_size = maxb / n_nodes;
    min_block_size = minb / n_nodes;


    if (rerun) {

        if (parallel) {
            node_comm_scatter_data();
        }

        this->n_c = data.n_elem;

        if (is_master) {
            if (max_block_size > this->n_c/2) {
                std::cout << "invalid local max block size " << max_block_size << std::endl;
                std::cout << "max block size must not be greater than " << n_nodes * this->n_c/2 << std::endl;
                exit(1);
            }

            if (min_block_size < 1) {
                std::cout << "invalid local min block size " << min_block_size << std::endl;
                std::cout << "min block size must not be lower than n_nodes=" << n_nodes << std::endl;
                exit(1);
            }

            if (n_block_samples > max_block_size - min_block_size) {
                std::cout << "invalid amount of block samples " << n_block_samples << std::endl;
                std::cout << "block samples must be lower or equal " << max_block_size - min_block_size << std::endl;
                exit(1);
            }
        }

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

    int block_size, n;
    double error, var, mean;

    error = 0;
    get_initial_error();

    arma::Row<int> block_sizes = arma::zeros<arma::Row<int> >(n_block_samples);

    get_unique_blocks(block_sizes, n);
    
    for (int j = 0; j < n; j++) {
        block_size = block_sizes(j);

        block_data(block_size, var, mean);

        var = combine_variance(var, mean);

        if (is_master) {
            error = sqrt(var / ((n_nodes * n_c) / block_size - 1.0));
            
            file << block_size * n_nodes << "\t" << error << std::endl;
            if (j % 9 == 0) {
                std::cout << "\rBlocking progress: " << (double) (j + 1) / n * 100 << "%";
                std::cout.flush();
            }
        }
    }
    
    if (is_master) std::cout << " Done." << std::endl;
    return error;
}

void Blocking::block_data(int block_size, double &var, double &mean) {
    using namespace arma;

    double block_mean;


    mean = 0;
    double mean2 = 0;

    int n_b = n_c / block_size;
    for (int j = 0; j < n_b; j++) {
        block_mean = arma::mean(data(span(j*block_size, (j + 1) * block_size - 1)));

        mean += block_mean;
        mean2 += block_mean*block_mean;
    }

    mean /= n_b;
    
    var = mean2/(n_b-1) - n_b*mean*mean/(n_b-1);
}

void Blocking::get_initial_error() {
    using namespace std;
    
    double var = arma::var(data);
    double mean = arma::mean(data);
    var = combine_variance(var, mean);
    mean = combine_mean(mean, data.n_elem);
    if (is_master) cout <<"Initial mean: "<< setprecision(6) << fixed  << mean << endl;
    if (is_master) cout << "Initial stddev: " << sqrt(var / (n_nodes * n_c - 1)) << endl;
}

void Blocking::get_unique_blocks(arma::Row<int>& block_sizes, int& n) {

    int block_step_length = (max_block_size - min_block_size) / (n_block_samples-1);

    for (int j = 0; j < n_block_samples; j++) {
        block_sizes(j) = min_block_size + j * block_step_length;
    }
    block_sizes = arma::unique(block_sizes);
    
    n = block_sizes.n_elem;
}
