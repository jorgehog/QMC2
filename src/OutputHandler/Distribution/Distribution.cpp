/* 
 * File:   Distribution.cpp
 * Author: jorgehog
 *
 * Created on 3. sept 2012, 13:17
 */

#include "../../QMCheaders.h"

Distribution::Distribution(std::string filename,
        std::string path,
        bool parallel,
        int my_rank,
        int num_procs)
: OutputHandler(filename, path, parallel, my_rank, num_procs) {

}

void Distribution::dump() {

    if ((qmc->cycle > qmc->n_c / 2) && (qmc->cycle % 100 == 0)) {
        for (int i = 0; i < qmc->n_p; i++) {
            for (int j = 0; j < qmc->dim; j++) {
                if (j == qmc->dim - 1) {
                    file << qmc->original_walker->r(i, j);
                } else {
                    file << qmc->original_walker->r(i, j) << " ";
                }
            }
            file << endl;
        }
    }

}