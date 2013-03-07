/* 
 * File:   Distribution.cpp
 * Author: jorgehog
 *
 * Created on 3. sept 2012, 13:17
 */

#include "../../QMCheaders.h"

Distribution::Distribution(ParParams & pp, std::string path, std::string name)
: OutputHandler("", path, pp.parallel, pp.node, pp.n_nodes) {
    this->name = name;
}

void Distribution::dump() {
    qmc->save_distribution();
}

void Distribution::finalize() {
    generate_distribution();
}

void Distribution::generate_distribution() {

    using namespace arma;

    //scrap out all the over-allocated space (DMC)
    qmc->dist.resize(qmc->last_inserted, qmc->dim);

    int n_tot;
    int n = qmc->dist.n_rows;

    int x_i;
    int y_i;

    double x, y;

#ifdef MPI_ON
    MPI_Allreduce(&n, &n_tot, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#endif

    colvec R = arma::sqrt(sum(qmc->dist % qmc->dist, 1));

    double stretch = 3;
    double mean_r = qmc->error_estimator->combine_mean(mean(R), n, n_tot);

    double bin_edge = stretch*mean_r;

    int N = 200;

    umat distribution = zeros<umat > (N, N);

    double dr = 2 * bin_edge / (N - 1);

    for (int ni = 0; ni < n; ni++) {

        x = qmc->dist(ni, 0);
        y = qmc->dist(ni, 1);

        x_i = (N / 2 + int(x / dr))*(x > 0) + (N / 2 + int(x / dr) - 1)*(x <= 0);
        y_i = (N / 2 + int(y / dr))*(y > 0) + (N / 2 + int(y / dr) - 1)*(y <= 0);

        if (x_i < 0 || x_i >= N) {
            continue;
        } else if (y_i < 0 || y_i >= N) {
            continue;
        }

        distribution(x_i, y_i)++;

    }

#ifdef MPI_ON
    if (node == 0) {
        MPI_Reduce(MPI_IN_PLACE, distribution.memptr(), N*N, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    } else {
        MPI_Reduce(distribution.memptr(), MPI_DATATYPE_NULL, N*N, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    }
#endif

    qmc->dist.reset();

    if (node == 0) {
        cout << n_tot << endl;
        mat normalized_dist = conv_to< mat>::from(distribution);
        normalized_dist /= (accu(normalized_dist) * dr * dr);

        s << path << "walker_positions/dist_out_" << name << ".arma";
        normalized_dist.save(s.str());
        s.str(std::string());

        normalized_dist.reset();

    }


    distribution.reset();

}

