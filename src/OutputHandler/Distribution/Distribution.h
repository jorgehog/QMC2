/* 
 * File:   Distribution.h
 * Author: jorgehog
 *
 * Created on 3. sept 2012, 13:17
 */

#ifndef DISTRIBUTION_H
#define	DISTRIBUTION_H

class Distribution : public OutputHandler {
protected:
    int i;
public:
    Distribution(ParParams &, std::string filename = "dist_out",
            std::string path = "./");

    virtual void dump();

};


#endif	/* DISTRIBUTION_H */