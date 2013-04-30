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
        
    } else if (qnum == 5) {    
    
        //-Z*r*(2*k^2*r^2 - 30*k*r + 81)/(3*(2*k^2*r^2 - 18*k*r + 27))
        
        dphi = -Z*walker->get_r_i(i)*(2*(*k2)*walker->get_r_i2(i) - 30*(*k)*walker->get_r_i(i) + 81)/(3*(2*(*k2)*walker->get_r_i2(i) - 18*(*k)*walker->get_r_i(i) + 27));
        
    } else if (qnum == 6) {    
    
        //-Z*r*(k*r - 9)/(3*(k*r - 6))
        
        dphi = -Z*walker->get_r_i(i)*((*k)*walker->get_r_i(i) - 9)/(3*((*k)*walker->get_r_i(i) - 6));
        
    } else if (qnum == 7) {    
    
        //-Z*r*(k*r - 9)/(3*(k*r - 6))
        
        dphi = -Z*walker->get_r_i(i)*((*k)*walker->get_r_i(i) - 9)/(3*((*k)*walker->get_r_i(i) - 6));
        
    } else if (qnum == 8) {    
    
        //-Z*r*(k*r - 9)/(3*(k*r - 6))
        
        dphi = -Z*walker->get_r_i(i)*((*k)*walker->get_r_i(i) - 9)/(3*((*k)*walker->get_r_i(i) - 6));
        
    } else if (qnum == 9) {    
    
        //-Z*r/3
        
        dphi = -Z*walker->get_r_i(i)/3;
        
    } else if (qnum == 10) {    
    
        //-Z*r/3
        
        dphi = -Z*walker->get_r_i(i)/3;
        
    } else if (qnum == 11) {    
    
        //-Z*r/3
        
        dphi = -Z*walker->get_r_i(i)/3;
        
    } else if (qnum == 12) {    
    
        x = walker->r(i, 0);
        y = walker->r(i, 1);
    
        x2 = x*x;
        y2 = y*y;
    
        //-Z*r*(x - y)*(x + y)/(3*(x^2 - y^2))
        
        dphi = -Z*walker->get_r_i(i)*(x - y)*(x + y)/(3*(x2 - y2));
        
    } else if (qnum == 13) {    
    
        //-Z*r/3
        
        dphi = -Z*walker->get_r_i(i)/3;
        
    }
    
    return dphi;
    
}