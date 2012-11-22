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

    Distribution(std::string filename = "dist_out",
            std::string path = "./",
            bool parallel = false,
            int my_rank = 0,
            int n_nodes = 1
            );

    virtual void dump();

};


#endif	/* DISTRIBUTION_H */