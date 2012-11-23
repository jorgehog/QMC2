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

    using namespace arma;

    rowvec tmp_data;
    rowvec master_data;


    if (node == 0) {
        master_data = zeros<rowvec > (n * n_nodes);
    } else {
        tmp_data = data;
    }





    cout << mean(data(span(0, n - 1))) << endl;


    MPI_Gather(tmp_data.memptr(), n, MPI_DOUBLE, master_data.memptr(), n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (node == 0) cout << mean(data) << " " << data.n_elem << endl;

    if (node != 0) {
        data.clear();
    } else {
        data = master_data;
    }

    tmp_data.clear();
    master_data.clear();


#endif
}
