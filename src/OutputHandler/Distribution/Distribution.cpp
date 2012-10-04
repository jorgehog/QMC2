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
    this->is_vmc = true;
}

void Distribution::dump() {

    if ((vmc->cycle > vmc->n_c / 2) && (vmc->cycle % 200 == 0)) {
        for (int i = 0; i < vmc->n_p; i++) {
            for (int j = 0; j < vmc->dim; j++) {
                if (j == vmc->dim - 1) {
                    file << vmc->original_walker->r(i, j);
                } else {
                    file << vmc->original_walker->r(i, j) << " ";
                }
            }
            file << endl;
        }
    }

}