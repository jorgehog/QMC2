/* 
 * File:   ErrorEstimator.cpp
 * Author: jorgmeister
 * 
 * Created on October 29, 2012, 3:58 PM
 */

#include "../QMCheaders.h"

ErrorEstimator::ErrorEstimator() {

}

ErrorEstimator::ErrorEstimator(int n_c,
        std::string filename,
        std::string path,
        bool parallel,
        int my_rank, int num_procs) {
    this->n_c = n_c;
    i = 0;
    data = arma::zeros<arma::rowvec > (n_c);

    this->my_rank = my_rank;
    this->num_procs = num_procs;
    this->parallel = parallel;

    this->filename = filename;
    this->path = path;

    if (parallel) {
        filename = filename + boost::lexical_cast<std::string > (my_rank);
    }

    this->file.open(((path + filename) + ".dat").c_str());
    
    this->do_output = true;
}





