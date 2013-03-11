/* 
 * File:   Distribution.h
 * Author: jorgehog
 *
 * Created on 3. sept 2012, 13:17
 */

#ifndef DISTRIBUTION_H
#define	DISTRIBUTION_H

class Distribution : public OutputHandler {
public:
    Distribution(ParParams &, std::string path, std::string name);

    void dump();
    void finalize();
    void rerun(int n_p, int N, double bin_edge);

private:

    int dim;

    std::string name;
    void generate_distribution(arma::mat & dist,
            int n_p,
            double bin_edge = 0,
            int N = 400,
            bool rerun = false);

    void post_pointer_init() {
        this->dim = qmc->dim;
    }

};


#endif	/* DISTRIBUTION_H */