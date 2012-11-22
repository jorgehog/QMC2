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
        int n_nodes,
        int n_b,
        int maxb,
        int minb,
        bool rerun)
: ErrorEstimator(n_c, filename, path, parallel, my_rank, n_nodes, rerun) {
    //    int step = 1;
    n_block_samples = n_b;
    max_block_size = maxb;
    min_block_size = minb;
    //    n_block_samples = (max_block_size - min_block_size) / step;
}

double Blocking::estimate_error() {
    using namespace std;

    int block_size, block_step_length;
    double error;

    block_step_length = (max_block_size - min_block_size) / n_block_samples;

    cout << "Initial stddev: " << sqrt(var(data) / (data.n_elem - 1)) << endl;

    for (int j = 0; j < n_block_samples; j++) {
        block_size = min_block_size + j*block_step_length;

        error = block_data(block_size);
        file << block_size << "\t" << error << std::endl;
    }

    finalize();

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

    for (int j = 0; j < n_b; j++) {
        block_mean = sum(data(span(j*block_size, (j + 1) * block_size - 1))) / block_size;

        mean += block_mean;
        mean2 += block_mean*block_mean;
    }

    mean2 /= (n_b);
    mean /= (n_b);

    return sqrt((mean2 - mean * mean) / ((n_c / block_size) - 1.0));
}

