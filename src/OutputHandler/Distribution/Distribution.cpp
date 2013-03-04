/* 
 * File:   Distribution.cpp
 * Author: jorgehog
 *
 * Created on 3. sept 2012, 13:17
 */

#include "../../QMCheaders.h"

Distribution::Distribution(ParParams & pp, std::string path)
: OutputHandler("", path, pp.parallel, pp.node, pp.n_nodes) {

}

void Distribution::dump() {

    qmc->save_distribution();

}

void Distribution::finalize() {
    s << path << "walker_positions/dist_out_" << qmc->name << node << ".arma";
    qmc->dist.save(s.str());
    qmc->dist.reset();
    s.str(std::string());
}