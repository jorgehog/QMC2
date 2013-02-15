basis_functions[0] = new hydrogenic_0(k, k2, r22d, r2d, exp_factor_n1);
basis_functions[1] = new hydrogenic_1(k, k2, r22d, r2d, exp_factor_n2);
basis_functions[2] = new hydrogenic_2(k, k2, r22d, r2d, exp_factor_n2);
basis_functions[3] = new hydrogenic_3(k, k2, r22d, r2d, exp_factor_n2);
basis_functions[4] = new hydrogenic_4(k, k2, r22d, r2d, exp_factor_n2);

dell_basis_functions[0][0] = new dell_hydrogenic_0_x(k, k2, r22d, r2d, exp_factor_n1);
dell_basis_functions[1][0] = new dell_hydrogenic_0_y(k, k2, r22d, r2d, exp_factor_n1);
dell_basis_functions[2][0] = new dell_hydrogenic_0_z(k, k2, r22d, r2d, exp_factor_n1);
dell_basis_functions[0][1] = new dell_hydrogenic_1_x(k, k2, r22d, r2d, exp_factor_n2);
dell_basis_functions[1][1] = new dell_hydrogenic_1_y(k, k2, r22d, r2d, exp_factor_n2);
dell_basis_functions[2][1] = new dell_hydrogenic_1_z(k, k2, r22d, r2d, exp_factor_n2);
dell_basis_functions[0][2] = new dell_hydrogenic_2_x(k, k2, r22d, r2d, exp_factor_n2);
dell_basis_functions[1][2] = new dell_hydrogenic_2_y(k, k2, r22d, r2d, exp_factor_n2);
dell_basis_functions[2][2] = new dell_hydrogenic_2_z(k, k2, r22d, r2d, exp_factor_n2);
dell_basis_functions[0][3] = new dell_hydrogenic_3_x(k, k2, r22d, r2d, exp_factor_n2);
dell_basis_functions[1][3] = new dell_hydrogenic_3_y(k, k2, r22d, r2d, exp_factor_n2);
dell_basis_functions[2][3] = new dell_hydrogenic_3_z(k, k2, r22d, r2d, exp_factor_n2);
dell_basis_functions[0][4] = new dell_hydrogenic_4_x(k, k2, r22d, r2d, exp_factor_n2);
dell_basis_functions[1][4] = new dell_hydrogenic_4_y(k, k2, r22d, r2d, exp_factor_n2);
dell_basis_functions[2][4] = new dell_hydrogenic_4_z(k, k2, r22d, r2d, exp_factor_n2);

lapl_basis_functions[0] = new lapl_hydrogenic_0(k, k2, r22d, r2d, exp_factor_n1);
lapl_basis_functions[1] = new lapl_hydrogenic_1(k, k2, r22d, r2d, exp_factor_n2);
lapl_basis_functions[2] = new lapl_hydrogenic_2(k, k2, r22d, r2d, exp_factor_n2);
lapl_basis_functions[3] = new lapl_hydrogenic_3(k, k2, r22d, r2d, exp_factor_n2);
lapl_basis_functions[4] = new lapl_hydrogenic_4(k, k2, r22d, r2d, exp_factor_n2);
