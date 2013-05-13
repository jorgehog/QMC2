/* 
 * File:   Walker.cpp
 * Author: jorgehog
 * 
 * Created on 30. mars 2012, 16:50
 */

#include "../QMCheaders.h"

Walker::Walker(int n_p, int dim, bool alive) {
    using namespace arma;

    this->dim = dim;
    this->n_p = n_p;
    this->n2 = ceil(n_p / 2.0);

    is_murdered = !alive;

    r = zeros<mat > (n_p, dim);
    r_rel = zeros<mat > (n_p, n_p);
    qforce = zeros<mat > (n_p, dim);
    jast_grad = zeros<mat > (n_p, dim);
    spatial_grad = zeros<mat > (n_p, dim);

    r2 = zeros(1, n_p);
    abs_r = zeros(1, n_p);

    lapl_sum = 0;
    spatial_ratio = 0;
    E = 0;
    inv = zeros<mat > (n2, n_p);

    phi = zeros<mat > (n_p, n2);
    dell_phi = arma::field<mat > (n_p, 1);
    for (int i = 0; i < n_p; i++) {
        dell_phi(i) = zeros<mat > (dim, n2);
    }

    dJ = zeros<cube > (n_p, n_p, dim);

}

void Walker::calc_r_i2(int i) {

    double r2i = 0;

    for (int j = 0; j < dim; j++) {
        r2i += r(i, j) * r(i, j);
    }

    r2(i) = r2i;

}

void Walker::make_rel_matrix() {
    int i, j;

    for (i = 0; i < n_p - 1; i++) {
        for (j = i + 1; j < n_p; j++) {
            r_rel(i, j) = r_rel(j, i) = calc_r_rel(i, j);
        }
    }

}

double Walker::calc_r_rel(int i, int j) const {
    double r_ij, tmp;

    r_ij = 0;
    for (int k = 0; k < dim; k++) {
        tmp = (r(i, k) - r(j, k));
        r_ij += tmp*tmp;
    }
    r_ij = sqrt(r_ij);

    return r_ij;
}

void Walker::print(std::string header) const {
    using namespace std;

    cout << endl;
    cout << "---- ---- ---- " << header << " ---- ---- ----" << endl;

    cout << "r\n" << this->r << endl;
    cout << "r_rel\n" << this->r_rel << endl;
    cout << "r2\n" << this->r2 << endl;
    cout << "r\n" << this->abs_r << endl;
    cout << "S grad\n" << this->spatial_grad << endl;

    cout << "J grad\n" << this->jast_grad << endl;
    cout << "J derivs\n" << this->dJ << endl;

    cout << "Qforce\n" << this->qforce << endl;
    cout << "inv\n" << this->inv << endl;

    cout << "Phi\n" << this->phi << endl;
    cout << "DellPhi\n" << this->dell_phi << endl;

    cout << "Lapl_sum|S ratio|E\t" << this->lapl_sum << "|";
    cout << this->spatial_ratio << "|";
    cout << this->E << endl;

    cout << "n_p|dim|n2\t" << this->n_p << "|";
    cout << this->dim << "|";
    cout << this->n2 << endl;

    cout << "Walker is currently ";
    if (is_murdered) {
        cout << "dead." << endl;
    } else {
        cout << "alive." << endl;
    }

    cout << "---- ---- ---- ---- ---- ---- ----" << endl;
    cout << endl;
}

void Walker::send_soul(int dest) {
#ifdef MPI_ON

    MPI_Send(&E, 1, MPI_DOUBLE, dest, 0, MPI_COMM_WORLD);
    MPI_Send(r.memptr(), r.n_elem, MPI_DOUBLE, dest, 1, MPI_COMM_WORLD);
    MPI_Send(r_rel.memptr(), r_rel.n_elem, MPI_DOUBLE, dest, 2, MPI_COMM_WORLD);
    MPI_Send(r2.memptr(), r2.n_elem, MPI_DOUBLE, dest, 3, MPI_COMM_WORLD);
    MPI_Send(phi.memptr(), phi.n_elem, MPI_DOUBLE, dest, 4, MPI_COMM_WORLD);
    MPI_Send(spatial_grad.memptr(), spatial_grad.n_elem, MPI_DOUBLE, dest, 5, MPI_COMM_WORLD);

    //if jastrow
    MPI_Send(dJ.memptr(), dJ.n_elem, MPI_DOUBLE, dest, 6, MPI_COMM_WORLD);
    MPI_Send(jast_grad.memptr(), jast_grad.n_elem, MPI_DOUBLE, dest, 7, MPI_COMM_WORLD);

    //if IS
    MPI_Send(qforce.memptr(), qforce.n_elem, MPI_DOUBLE, dest, 8, MPI_COMM_WORLD);

    //if fermions
    MPI_Send(inv.memptr(), inv.n_elem, MPI_DOUBLE, dest, 9, MPI_COMM_WORLD);

    for (int i = 0; i < n_p; i++) {
        MPI_Send(dell_phi(i).memptr(), dell_phi(i).n_elem, MPI_DOUBLE, dest, i + 10, MPI_COMM_WORLD);
    }

    kill();

#endif
}

void Walker::recv_soul(int root) {
#ifdef MPI_ON

    MPI_Recv(&E, 1, MPI_DOUBLE, root, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(r.memptr(), r.n_elem, MPI_DOUBLE, root, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(r_rel.memptr(), r_rel.n_elem, MPI_DOUBLE, root, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(r2.memptr(), r2.n_elem, MPI_DOUBLE, root, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(phi.memptr(), phi.n_elem, MPI_DOUBLE, root, 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(spatial_grad.memptr(), spatial_grad.n_elem, MPI_DOUBLE, root, 5, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    //if jastrow
    MPI_Recv(dJ.memptr(), dJ.n_elem, MPI_DOUBLE, root, 6, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(jast_grad.memptr(), jast_grad.n_elem, MPI_DOUBLE, root, 7, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    //if IS
    MPI_Recv(qforce.memptr(), qforce.n_elem, MPI_DOUBLE, root, 8, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    //if fermions
    MPI_Recv(inv.memptr(), inv.n_elem, MPI_DOUBLE, root, 9, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    for (int i = 0; i < n_p; i++) {
        MPI_Recv(dell_phi(i).memptr(), dell_phi(i).n_elem, MPI_DOUBLE, root, i + 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }


    ressurect();

#endif
}