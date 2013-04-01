/* 
 * File:   function.h
 * Author: jorgehog
 *
 * Created on 27. juni 2012, 13:17
 */

#ifndef BASISFUNCTIONS_H
#define	BASISFUNCTIONS_H

/*! \brief The Superclass shell for orbital basis functions.
 * 
 * Each single particle orbital has it's own implementation as a subclass of this class.
 * A set of orbitals can then be loaded into the Orbitals basis_function vectors.
 * An orbitalGenerator script is supplied to autogenerate CPP files using this class.
 * @see Orbitals::basis_functions
 * @see Orbitals::dell_basis_functions
 * @see Orbitals::lapl_basis_functions
 */

class BasisFunctions {
public:
    BasisFunctions();
    
    //! The method representing the orbitals functional expression.
    /*!
     * @param walker The Walker whose position the orbital is to be evaluted at.
     * @param i The particle to be evaluated (Single particle function).
     */
    virtual double eval(const Walker* walker, int i) = 0;
};

#endif	/* BASISFUNCTIONS_H */

