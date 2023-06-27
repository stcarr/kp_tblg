function [t] = dft_interlayer_coupling(vector, theta_row, theta_col, a)
    % old function call goes like: 
    % dft_interlayer_coupling(orbit_row, orbit_col, vector, theta_row, theta_col, a)
    
    r = sqrt(vector(:,1).^2 + vector(:,2).^2);
    ac = atan2(vector(:,2), vector(:,1));

    % I think we do not need A/B suborbital info since
    % theta_row/theta_col is the LOCAL angle to NN bond...
    % i.e. just treat everything like an A orbital...

    % theta21 (angle to bond on sheet 1)
    % In our convention the NN in-plane bond is along the x-dir
    theta21 = ac - theta_row;% - pi/6;
    %{
    if (orbit_row == 2) % B orbital
        theta21 = theta21 + pi/6;
    else  % A orbital
        theta21 = theta21 - pi/6;
    end
    %}
    % theta12 (angle to bond on sheet 2)
    theta12 = (ac+pi) - theta_col;% - pi/6;
    %{
    if (orbit_col == 2) % B orbital
        theta12 = theta12 + pi/6;
    else % A orbital
        theta12 = theta12 - pi/6;
    end
    %}
    rs = r/a; % a = lattice param
    z_eq = 3.35; % Equilbrium AB height from our DFT Calculations
    z_eps = (abs(vector(:,3)) - z_eq)/z_eq;
    z_eps_sq = z_eps.*z_eps;
    % New coefficients (with compression dependence fitting)
    lambda_0 =  0.310 - 1.882*z_eps + 7.741*z_eps_sq;
    xi_0     =  1.750 - 1.618*z_eps + 1.848*z_eps_sq;
    kappa_0  =  1.990 + 1.007*z_eps + 2.427*z_eps_sq;

    lambda_3 = -0.068 + 0.399*z_eps - 1.739*z_eps_sq;
    xi_3     =  3.286 - 0.914*z_eps + 12.011*z_eps_sq;
    x_3      =  0.500 + 0.322*z_eps + 0.908*z_eps_sq;

    lambda_6 = -0.008 + 0.046*z_eps - 0.183*z_eps_sq;
    xi_6     =  2.272 - 0.721*z_eps - 4.414*z_eps_sq;
    x_6      =  1.217 + 0.027*z_eps - 0.658*z_eps_sq;
    kappa_6  =  1.562 - 0.371*z_eps - 0.134*z_eps_sq;

    % Old coefficients (From: S. Fang, Phys. Rev. B 93, 235153 - Published 27 June 2016)
    %{
    lambda_0 =  0.3155;
    xi_0     =  1.7543;
    kappa_0  =  2.0010;

    lambda_3 = -0.0688;
    xi_3     =  3.4692;
    x_3      =  0.5212;

    lambda_6 = -0.0083;
    xi_6     =  2.8764;
    x_6      =  1.5206;
    kappa_6  =  1.5731;
    %}

    V0 = lambda_0 .*         exp(-xi_0.*(rs    ).*(rs    )) .* cos(kappa_0.*rs);
    V3 = lambda_3 .* rs.*rs .* exp(-xi_3.*(rs-x_3).*(rs-x_3));                       
    V6 = lambda_6 .*         exp(-xi_6.*(rs-x_6).*(rs-x_6)) .* sin(kappa_6.*rs);

    t = V0 + V3.*(cos(3*theta12)+cos(3*theta21)) + V6.*(cos(6*theta12)+cos(6*theta21));
    %t = V0; % no angular component
end