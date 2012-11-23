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
        int node,
        int n_nodes)
: OutputHandler(filename, path, parallel, node, n_nodes) {
    i = 0;
    this->is_vmc = true;
}

void Distribution::dump() {

    if ((vmc->cycle > vmc->n_c / 2) && (vmc->cycle % 100 == 0)) {
        s << path << "walker_positions/" << filename << node <<"_"<< i <<".arma";
        vmc->original_walker->r.save(s.str());
        s.str(std::string());
        i++;
    }

}