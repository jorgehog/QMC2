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

    this->node = node;
    this->n_nodes = n_nodes;
    this->parallel = parallel;
    this->is_master = (node == 0);

    this->filename = filename;
    this->path = path;

    this->rerun = rerun;
    if (rerun) {

        if (is_master) {

            bool success = data.load(path + (filename + "_RAWDATA.arma"));
            if (!success) exit(1);

            this->n_c = data.n_elem / n_nodes;

        }

        data_to_file = false;

    } else {
        data = arma::zeros<arma::rowvec > (n_c);
        data_to_file = true;
    }



    if (parallel) {
        filename = filename + boost::lexical_cast<std::string > (node);
    }

    //will be overwritten (true) if a file is opened.
    output_to_file = false;

}

void ErrorEstimator::init_file() {
    output_to_file = true;
    if (is_master) this->file.open(((path + filename) + ".dat").c_str());
}

void ErrorEstimator::finalize() {

    if (data_to_file) {
        node_comm_gather_data();
    }

    if (data_to_file && is_master) {
        data.save(path + (filename + "_RAWDATA.arma"));
    }

    if (output_to_file && is_master) file.close();

    data.clear();

}

void ErrorEstimator::node_comm_gather_data() {
#ifdef MPI_ON
    if (parallel) {

        int n = data.n_elem;

        if (is_master) {
            data.resize(n * n_nodes);
            MPI_Gather(MPI_IN_PLACE, n, MPI_DOUBLE, data.memptr(), n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        } else {
            MPI_Gather(data.memptr(), n, MPI_DOUBLE, NULL, 0, MPI_DATATYPE_NULL, 0, MPI_COMM_WORLD);
        }
    }
#endif
}

void ErrorEstimator::node_comm_scatter_data() {
#ifdef MPI_ON
    if (parallel) {
        using namespace arma;

        MPI_Bcast(&(this->n_c), 1, MPI_INT, 0, MPI_COMM_WORLD);

        rowvec raw_data;

        if (is_master) {

            raw_data = zeros<rowvec > (this->n_c);
            MPI_Scatter(data.memptr(), this->n_c, MPI_DOUBLE, raw_data.memptr(), this->n_c, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            data.clear();
            data = raw_data;
            raw_data.clear();
        } else {
            data = zeros<rowvec > (this->n_c);
            MPI_Scatter(new double(), 1, MPI_DOUBLE, data.memptr(), this->n_c, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        }
    }
#endif
}

double ErrorEstimator::combine_variance(double var, double mean) {

    int n = data.n_elem;

#ifdef MPI_ON

    if (parallel) {

        double combined_mean;

        MPI_Allreduce(&mean, &combined_mean, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        combined_mean /= n_nodes;

        double combined_var = n * (mean - combined_mean)*(mean - combined_mean) + (n - 1) * var;

        if (is_master) {
            MPI_Reduce(MPI_IN_PLACE, &combined_var, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
            var = combined_var / (n * n_nodes - 1);
        } else {
            MPI_Reduce(&combined_var, new double(), 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        }

    }

#endif

    return var;
}