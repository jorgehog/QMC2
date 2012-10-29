/* 
 * File:   Blocking.cpp
 * Author: jorgmeister
 * 
 * Created on October 29, 2012, 3:58 PM
 */

#include "../../QMCheaders.h"

Blocking::Blocking(int n_c,
        std::string filename,
        std::string path,
        bool parallel,
        int my_rank,
        int num_procs)
: ErrorEstimator(n_c, filename, path, parallel, my_rank, num_procs) {

}

double Blocking::estimate_error() {
    int block_size, block_step_length;
    double error;

    block_step_length = (max_block_size - min_block_size) / n_block_samples;

    for (int i = 0; i < n_block_samples; i++) {
        block_size = min_block_size + i*block_step_length;
        error = block_data(block_size);

        file << block_size << "\t" << error << std::endl;

    }

    data.clear();
    
    return error;
}

double Blocking::block_data(int block_size) {
    using namespace arma;

    int n_b;
    double block_mean;
    double mean;
    double mean2;

    n_b = n_c / block_size;

    mean = 0;
    mean2 = 0;
    for (int i = 0; i < n_b; i++) {
        block_mean = sum(data(span(i*block_size, i * block_size + block_size)));
        mean += block_mean;
        mean2 += block_mean*block_mean;
    }

    mean2 /= block_size;
    mean /= block_size;

    return sqrt((mean2 - mean * mean) / ((n_c / block_size) - 1.0));
}

