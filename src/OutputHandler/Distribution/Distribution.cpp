/* 
 * File:   Distribution.cpp
 * Author: jorgehog
 *
 * Created on 3. sept 2012, 13:17
 */

#include "../../QMCheaders.h"

Distribution::Distribution(ParParams & pp, std::string filename,
        std::string path)
: OutputHandler(filename, path, pp.parallel, pp.node, pp.n_nodes) {
    i = 0;
    this->is_vmc = true;
}

void Distribution::dump() {

    if ((vmc->cycle >= vmc->n_c / 2) && (vmc->cycle % 100 == 0)) {
        
        //s: unique file name for the current snapshot of the diffusing walker
        s << path << "walker_positions/" << filename << node <<"_"<< i <<".arma";
        
        //saves the position matrix to file
        vmc->original_walker->r.save(s.str());
        
        //clears the stringstream and iterates identifyer 'i'
        s.str(std::string());
        i++;
    }
}