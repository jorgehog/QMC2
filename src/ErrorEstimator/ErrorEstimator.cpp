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
        int node, int n_nodes, bool rerun) {

    this->n_c = n_c;
    i = 0;

    this->rerun = rerun;
    if (rerun) {

        bool success = data.load(path + (filename + "_RAWDATA.arma"));
        if (!success) exit(1);

        this->n_c = data.n_cols;
        data_to_file = false;

    } else {
        data = arma::zeros<arma::rowvec > (n_c);
    }

    this->node = node;
    this->n_nodes = n_nodes;
    this->parallel = parallel;
    this->is_master = (node == 0);

    this->filename = filename;
    this->path = path;

    if (parallel) {
        filename = filename + boost::lexical_cast<std::string > (node);
    }

    //will be overwritten (true) if a file is opened.
    output_to_file = false;

}

void ErrorEstimator::init_file() {
    output_to_file = true;
    if (node == 0) this->file.open(((path + filename) + ".dat").c_str());
}

void ErrorEstimator::finalize() {
    if (data_to_file) data.save(path + (filename + "_RAWDATA.arma"));
    if (output_to_file) file.close();

    data.clear();

}

void ErrorEstimator::node_comm() {
#ifdef MPI_ON

    int n = data.n_elem;

    if (is_master) {
        data.resize(n * n_nodes);
        MPI_Gather(MPI_IN_PLACE, n, MPI_DOUBLE, data.memptr(), n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    } else {
        MPI_Gather(data.memptr(), n, MPI_DOUBLE, NULL, 0, MPI_DATATYPE_NULL, 0, MPI_COMM_WORLD);
        data.clear();
    }

#endif
}
