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
        int my_rank, int num_procs, std::string* infile) {

    this->n_c = n_c;
    i = 0;

    if (infile != NULL) {
        
        bool success = data.load((path + *infile));
        if (!success) exit(1);
        
        this->n_c = data.n_rows;

        to_file = false;
        do_output = false;

    } else {
        data = arma::zeros<arma::rowvec > (n_c); //devide by numprocs?
        this->do_output = true;
        this->to_file = true;
    }

    this->my_rank = my_rank;
    this->num_procs = num_procs;
    this->parallel = parallel;

    this->filename = filename;
    this->path = path;

    if (parallel) {
        filename = filename + boost::lexical_cast<std::string > (my_rank);
    }

    this->file.open(((path + filename) + ".dat").c_str());

}

void ErrorEstimator::finalize() {
    if (to_file) data.save(path + (filename + "_RAWDATA.rowvec"));
    data.clear();
    file.close();
}


