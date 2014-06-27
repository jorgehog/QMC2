#include "Distribution.h"

#include "../../structs.h"

#include "../../QMC/QMC.h"
#include "../../ErrorEstimator/ErrorEstimator.h"


using namespace QMC2;

Distribution::Distribution(ParParams & pp, std::string path, std::string name, bool silent)
: OutputHandler("", path, pp.parallel, pp.node, pp.n_nodes) {
    this->name = name;
    this->silent = silent;
}


void Distribution::detect_deadlock(const arma::mat& dist, int n_p, int dim, int n) {

    deadlock_x = dist(0, 0);
    deadlock_y = dist(0, 1);

    if (dim == 3) {
        deadlock_z = dist(0, 2);
    }
    locked = true;

    for (int i = n_p; i < n; i += n_p) {

        if (!is_deadlocked(dist, dim, i)) {
            locked = false;
            return;
        }

    }

    if ((node == 0) && (!silent)) std::cout << "Detected deadlock at " << deadlock_x << "  " << deadlock_y << "  " << deadlock_z << "  " << std::endl;
}

void Distribution::generate_distribution3D(arma::mat& dist,
        int n_p,
        double bin_edge,
        int N,
        bool rerun) {

    int x_i, y_i, z_i, r_i, n, n_tot;
    double x, y, z, r, dr, dr_R, stretch, mean_r;

    int dim = 3;
    
    using namespace arma;

    ucube distribution(N, N, N);
    uvec radial_dist = zeros<uvec > (N);
    mat tot_dist;

    vec R = arma::sqrt(sum(dist % dist, 1));

    n = dist.n_rows;

#ifdef MPI_ON
    ivec n_list = zeros<ivec > (n_nodes);

    MPI_Allgather(&n, 1, MPI_INT, n_list.memptr(), 1, MPI_INT, MPI_COMM_WORLD);

    n_tot = accu(n_list);
#else
    n_tot = n;
#endif

    detect_deadlock(dist, n_p, dim, n);

    //On fly calculation during QMC initialized by a binedge 0.
    if (bin_edge == 0) {
        stretch = 3;

        //calculate the bin edge and size based on the mean radius
        if (!locked) {
            mean_r = ErrorEstimator::combine_mean(mean(R), n, n_tot);
        } else {

            mean_r = 0;
            int k = 0;

            for (int i = 0; i < n; i++) {
                if (is_deadlocked(dist, dim, i)) continue;

                mean_r += R(i);
                k++;
            }

            mean_r = ErrorEstimator::combine_mean(mean_r / k, n, n_tot);

        }
        bin_edge = stretch*mean_r;
        //        if (node==0) cout << "pre " << mean_r << endl;
    }

    dr = 2 * bin_edge / (N - 1);

    for (int ni = 0; ni < n; ni++) {

        if (is_deadlocked(dist, dim, ni)) continue;

        x = dist(ni, 0);
        y = dist(ni, 1);
        z = dist(ni, 2);

        x_i = (N / 2 + int(x / dr))*(x > 0) + (N / 2 + int(x / dr) - 1)*(x <= 0);
        y_i = (N / 2 + int(y / dr))*(y > 0) + (N / 2 + int(y / dr) - 1)*(y <= 0);
        z_i = (N / 2 + int(z / dr))*(z > 0) + (N / 2 + int(z / dr) - 1)*(z <= 0);

        if (x_i < 0 || x_i >= N) {
            continue;
        } else if (y_i < 0 || y_i >= N) {
            continue;
        } else if (z_i < 0 || z_i >= N) {
            continue;
        }

        distribution(x_i, y_i, z_i)++;

    }

    dr_R = dr / 2;
    for (int ni = 0; ni < n; ni++) {

        if (is_deadlocked(dist, dim, ni)) continue;

        r = R(ni);

        r_i = r / dr_R;

        if (r_i >= N) {
            continue;
        }

        radial_dist(r_i)++;

    }




#ifdef MPI_ON
    if (node == 0) {

        MPI_Reduce(MPI_IN_PLACE, distribution.memptr(), N * N*N, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(MPI_IN_PLACE, radial_dist.memptr(), N, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);

        if (!rerun) {

            ivec displs = zeros<ivec > (n_nodes);

            double displ_sum = 0;
            for (int j = 0; j < n_nodes; j++) {
                displs(j) = displ_sum;
                displ_sum += n_list(j);
            }

            tot_dist = zeros<mat > (n_tot, dim);

            for (int i = 0; i < dim; i++) {
                MPI_Gatherv(dist.colptr(i), n, MPI_DOUBLE,
                        tot_dist.colptr(i), n_list.memptr(), displs.memptr(),
                        MPI_DOUBLE, 0, MPI_COMM_WORLD);
            }

            //            cout << as_scalar(mean(arma::sqrt(sum(tot_dist % tot_dist, 1)))) << endl;

            displs.reset();

        }

    } else {
        MPI_Reduce(distribution.memptr(), NULL, N * N*N, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(radial_dist.memptr(), NULL, N, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);

        if (!rerun) {
            for (int i = 0; i < dim; i++) {
                MPI_Gatherv(dist.colptr(i), n, MPI_DOUBLE,
                        NULL, NULL, NULL,
                        NULL, 0, MPI_COMM_WORLD);
            }
        }

        radial_dist.reset();
        distribution.reset();
    }
#else
    tot_dist(dist.memptr(), n, dim, false, true);
#endif

    n_list.reset();
    R.reset();

    if (node == 0) {

        if (!silent) cout << "3D Distribution calculated using " << n_tot << " samples." << endl;

        cube normalized_dist = conv_to< cube>::from(distribution);
        vec normalized_radd = conv_to< vec>::from(radial_dist);

        vec radial_axis = linspace(0, bin_edge, N);

        normalized_dist *= n_p / (accu(normalized_dist) * dr * dr * dr);


        //project out a symmetric axis and normalize (skip singularity)
        //        normalized_radd(span(1, N - 1)) /= radial_axis(span(1, N - 1));
        //
        //        normalized_radd(span(1, N - 1)) /= radial_axis(span(1, N - 1));
        //
        //        normalized_radd(0) = normalized_radd(1);

        //        cout << normalized_radd.max()/(2*datum::pi*accu(normalized_radd)*dr_R) << endl;
        normalized_radd *= n_p / (accu(normalized_radd) * dr_R);


        s << path << "walker_positions/dist_out_" << name + suffix << "_edge" << bin_edge << ".arma3D";
        normalized_dist.save(s.str());
        normalized_dist.reset();
        s.str(std::string());

        if (!rerun) {
            s << path << "walker_positions/dist_rawdata_" << name + suffix << ".arma";
            tot_dist.save(s.str());
            tot_dist.reset();
            s.str(std::string());
        }

        s << path << "walker_positions/radial_out_" << name + suffix << "_edge" << bin_edge << ".arma";
        normalized_radd.save(s.str());
        normalized_radd.reset();
        s.str(std::string());

        radial_axis.reset();
        radial_dist.reset();
        distribution.reset();

    }

}

void Distribution::generate_distribution2D(arma::mat & dist,
        int n_p,
        double bin_edge,
        int N,
        bool rerun) {

    using namespace arma;

    int dim = 2;
    
    int x_i, y_i, r_i, n, n_tot;
    double x, y, r, dr, dr_R, stretch, mean_r;



    umat distribution = zeros<umat > (N, N);
    uvec radial_dist = zeros<uvec > (N);
    mat tot_dist;

    vec R = arma::sqrt(sum(dist % dist, 1));
    //    vec R = arma::sqrt(sum(dist % dist, 1));

    n = dist.n_rows;

#ifdef MPI_ON
    ivec n_list = zeros<ivec > (n_nodes);

    MPI_Allgather(&n, 1, MPI_INT, n_list.memptr(), 1, MPI_INT, MPI_COMM_WORLD);

    n_tot = accu(n_list);
#else
    n_tot = n;
#endif

    detect_deadlock(dist, n_p, dim, n);

    //On fly calculation during QMC initialized by a binedge 0.
    if (bin_edge == 0) {
        stretch = 3;

        //calculate the bin edge and size based on the mean radius
        if (!locked) {
            mean_r = ErrorEstimator::combine_mean(mean(R), n, n_tot);
        } else {

            mean_r = 0;
            int k = 0;

            for (int i = 0; i < n; i++) {
                if (is_deadlocked(dist, dim, i)) continue;

                mean_r += R(i);
                k++;
            }

            mean_r = ErrorEstimator::combine_mean(mean_r / k, n, n_tot);

        }
        bin_edge = stretch*mean_r;
        //        if (node==0) cout << "pre " << mean_r << endl;
    }

    dr = 2 * bin_edge / (N - 1);

    for (int ni = 0; ni < n; ni++) {

        if (is_deadlocked(dist, dim, ni)) continue;

        x = dist(ni, 0);
        y = dist(ni, 1);

        x_i = (N / 2 + int(x / dr))*(x > 0) + (N / 2 + int(x / dr) - 1)*(x <= 0);
        y_i = (N / 2 + int(y / dr))*(y > 0) + (N / 2 + int(y / dr) - 1)*(y <= 0);

        if (x_i < 0 || x_i >= N) {
            continue;
        } else if (y_i < 0 || y_i >= N) {
            continue;
        }

        distribution(x_i, y_i)++;

    }

    dr_R = dr / 2;
    for (int ni = 0; ni < n; ni++) {

        if (is_deadlocked(dist, dim, ni)) continue;

        r = R(ni);

        r_i = r / dr_R;

        if (r_i >= N) {
            continue;
        }

        radial_dist(r_i)++;

    }

#ifdef MPI_ON
    if (node == 0) {

        MPI_Reduce(MPI_IN_PLACE, distribution.memptr(), N*N, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(MPI_IN_PLACE, radial_dist.memptr(), N, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);

        if (!rerun) {

            ivec displs = zeros<ivec > (n_nodes);

            double displ_sum = 0;
            for (int j = 0; j < n_nodes; j++) {
                displs(j) = displ_sum;
                displ_sum += n_list(j);
            }

            tot_dist = zeros<mat > (n_tot, dim);

            for (int i = 0; i < dim; i++) {
                MPI_Gatherv(dist.colptr(i), n, MPI_DOUBLE,
                        tot_dist.colptr(i), n_list.memptr(), displs.memptr(),
                        MPI_DOUBLE, 0, MPI_COMM_WORLD);
            }

            //            cout << as_scalar(mean(arma::sqrt(sum(tot_dist % tot_dist, 1)))) << endl;

            displs.reset();

        }

    } else {
        MPI_Reduce(distribution.memptr(), NULL, N*N, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(radial_dist.memptr(), NULL, N, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);

        radial_dist.reset();
        distribution.reset();

        if (!rerun) {
            for (int i = 0; i < dim; i++) {
                MPI_Gatherv(dist.colptr(i), n, MPI_DOUBLE,
                        NULL, NULL, NULL,
                        NULL, 0, MPI_COMM_WORLD);
            }
        }

    }
#else
    tot_dist(dist.memptr(), n, dim, false, true);
#endif

    n_list.reset();
    R.reset();

    if (node == 0) {

        if (!silent) cout << "Distribution calculated using " << n_tot << " samples." << endl;

        mat normalized_dist = conv_to< mat>::from(distribution);
        vec normalized_radd = conv_to< vec>::from(radial_dist);

        vec radial_axis = linspace(0, bin_edge, N);

        normalized_dist *= n_p / (accu(normalized_dist) * dr * dr);

        //        project out a symmetric axis and normalize(skip singularity)
        normalized_radd(span(1, N - 1)) /= radial_axis(span(1, N - 1));

        normalized_radd(0) = normalized_radd(1);

        //        cout << normalized_radd.max()/(2*datum::pi*accu(normalized_radd)*dr_R) << endl;
        normalized_radd *= n_p / (accu(normalized_radd) * dr_R);

        if (!rerun) {
            s << path << "walker_positions/dist_rawdata_" << name + suffix<< ".arma";
            tot_dist.save(s.str());
            tot_dist.reset();
            s.str(std::string());
        }

        s << path << "walker_positions/dist_out_" << name + suffix<< "_edge" << bin_edge << ".arma";
        normalized_dist.save(s.str());
        normalized_dist.reset();
        s.str(std::string());

        s << path << "walker_positions/radial_out_" << name + suffix<< "_edge" << bin_edge << ".arma";
        normalized_radd.save(s.str());
        normalized_radd.reset();
        s.str(std::string());

        radial_axis.reset();
        radial_dist.reset();
        distribution.reset();

    }

}

void Distribution::rerun(int n_p, int N, double bin_edge) {

    using namespace arma;

    mat dist;
    int n, dim;

    if (node == 0) {
        s << path << "walker_positions/dist_rawdata_" << name + suffix<< ".arma";
        dist.load(s.str());
        s.str(std::string());
        n = dist.n_rows / n_nodes;
        dim = dist.n_cols;

        if (n == 0) {
            cout << "Not enough data to scatter." << endl;
            exit(1);
        }
    }

#ifdef MPI_ON

    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&dim, 1, MPI_INT, 0, MPI_COMM_WORLD);

    mat raw_data;

    if (node == 0) {

        raw_data = zeros<mat > (n, dim);

        for (int i = 0; i < dim; i++) {
            MPI_Scatter(dist.colptr(i), n, MPI_DOUBLE, raw_data.colptr(i), n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        }

        dist.reset();
        dist = raw_data;
        raw_data.reset();

    } else {

        dist = zeros<mat > (n, dim);

        for (int i = 0; i < dim; i++) {
            MPI_Scatter(NULL, 0, MPI_DATATYPE_NULL, dist.colptr(i), n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        }
    }

#endif

    if (dim == 3) {
        generate_distribution3D(dist, n_p, bin_edge, N, true);
    } else {
        generate_distribution2D(dist, n_p, bin_edge, N, true);
    }

    dist.reset();

}

