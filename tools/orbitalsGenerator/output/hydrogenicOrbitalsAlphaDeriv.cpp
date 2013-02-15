double hydrogenicOrbitals::get_dell_alpha_phi(const Walker* walker, int qnum, int i){
    
    double dphi;
    
    if (qnum == 0) {    
    
        //-Z*r
        
        dphi = -Z*walker->get_r_i(i);
        
    } else if (qnum == 1) {    
    
        //-Z*r*(k*r - 4)/(2*(k*r - 2))
        
        dphi = -Z*walker->get_r_i(i)*((*k)*walker->get_r_i(i) - 4)/(2*((*k)*walker->get_r_i(i) - 2));
        
    } else if (qnum == 2) {    
    
        //-Z*r/2
        
        dphi = -Z*walker->get_r_i(i)/2;
        
    } else if (qnum == 3) {    
    
        //-Z*r/2
        
        dphi = -Z*walker->get_r_i(i)/2;
        
    } else if (qnum == 4) {    
    
        //-Z*r/2
        
        dphi = -Z*walker->get_r_i(i)/2;
        
    }
    
    return dphi;
    
}