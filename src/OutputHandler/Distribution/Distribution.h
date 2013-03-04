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
    Distribution(ParParams &, std::string path);

    void dump();
    void finalize();
    
};


#endif	/* DISTRIBUTION_H */