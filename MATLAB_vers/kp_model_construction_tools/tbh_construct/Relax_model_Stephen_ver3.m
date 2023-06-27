function [atom_pos_shift] = Relax_model_Stephen_ver3(rot_theta,sc_b1,sc_b2,sc_t1,sc_t2,top_layer,atom_pos_list)

% Conventions:
% For the rigid twist:
% bottom layer is at z=0 plane, 
% while top layer at z=3.35A plane.
% The top layer is twisted/rotated (counter clockwise, CCW) by rot_theta (radian) with respect to
% the bottom layer.

% what this code does: given a list of atomic positions, return the
% displacement vector (3d) for each one of them.

% input parameters:

% rot_theta : twist angle of the top layer, positive number in radian (not degree)
% one should expect the input angle to be from 0.5 degree to about 6
% degree.

% sc_b1, sc_b2 : supercell k space reciprocal lattice vectors, 1 by 2 vector, in
% units of 1/Angstrom or 1/a (a: graphene lattice constant ~2.46A)

% sc_t1, sc_t2 : supercell (real space) lattice vectors, 1 by 2 vector, in
% units of Angstrom or a (a: graphene lattice constant ~2.46A)

% top_layer: a single boolean value, true for the top layer, false for the bottom layer.

% atom_pos_list: a N by 2 matrix that is the list of unrelaxed atomic
% position (to be used to compute for the displacement vector for each one of them.)
% they are all top layer if the top_layer is true, or all bottom layer
% positions if otherwise.
% the unit of these positions will be the same as sc_t1/sc_t2, A or a.

% atom_pos_shift: the output, a N by 3 matrix for the corresponding x,y,z
% displacements for each atoms given, with the displacement distance in
% unit of A (Angstrom).

% also note that don't add z=3.35A to the top layer. I will deal with that
% somewhere.
% further note: careful for the dc term for the z shifts. I guess you want
% the minimal shift to be zero for the top layer, not something negative.

sizetmp=size(atom_pos_list);
atom_pos_shift=zeros(sizetmp(1),3);

% your code for interpolation/fft summation.
u_sign = 1;

if top_layer
    % code for computing with atoms in the TOP layers 
    u_sign = -1;    
else
    % code for computing with atoms in the BOTTOM layers 
end

% filepath is where the relaxation data is located!
%filepath = './sweep_coeffs2/';
%filepath = './0p1_30p0_k25_08-24-2018/';

filepath = './relax_data/zfull_0p1_30p0_k25_09-25-2018/';

% load the fourier coeff data
theta_list = dlmread([filepath 'thetas.txt']);
coeffs_x = dlmread([filepath 'coeffs_x.txt']);
coeffs_y = dlmread([filepath 'coeffs_y.txt']);
coeffs_z = dlmread([filepath 'coeffs_z.txt']);

% get the distance between theta samples
d_theta = theta_list(2) - theta_list(1);

% convert to degrees
tar_theta = rot_theta*180/pi;
%tar_theta = 0.4; %hardcoded for symm checking

% interpolate coefficients from nearest values in theta_list
[interp_diff, nearest_theta] = min(abs(tar_theta - theta_list));
if (tar_theta - theta_list(nearest_theta) < 0 )
   nearest_theta = nearest_theta - 1;
   interp_diff = -interp_diff+d_theta;
end
interp_diff = interp_diff/d_theta;

coeffs_1 = [coeffs_x(:,nearest_theta) coeffs_y(:,nearest_theta) coeffs_z(:,nearest_theta)];
coeffs_2 = [coeffs_x(:,nearest_theta+1) coeffs_y(:,nearest_theta+1) coeffs_z(:,nearest_theta+1)];
coeffs = coeffs_1*(1-interp_diff) + coeffs_2*interp_diff;

% Number of total atoms
N = length(atom_pos_list);

% set up some temporary vectors for storing displacements (might need to
% rotate them into a different basis at end)
disp_x = zeros(N,1);
disp_y = zeros(N,1);
disp_z = zeros(N,1);

% Used for finding rotated directions of the reciprocal lattice vectors
rot = [cosd(60) -sind(60);
        sind(60)  cosd(60)];
    
% idx keeps track of where we are in the coeffs array
idx = 1;

% max_k = 2 for testing (up to [2,1] component)!
% if max_k <  5, the Z component is not reliable at small angles...
max_k = 25;%25;
%max_k = 2;

% manually clean up the symmetry!
for i = 1:max_k %skips i = 0, the (0,0) constant term
    
   % first do the (m,n) to (m,m-n) symmetry) 
    
   pivot = i/2;

   if mod(i,2) == 0
       sep = 1.0;
   else
       sep = 0.5;
   end
   
   while (pivot + sep) < i
       tar_j_p = pivot + sep;
       tar_j_m = pivot - sep;
       
       n_hat = [-1/2,sqrt(3)/2]; %mirror-plane direction for u_x,u_y symm
       
       p_idx = 2 + (i-1)*i/2 + tar_j_p;
       m_idx = 2 + (i-1)*i/2 + tar_j_m;
       
       u_A = coeffs(m_idx,1:2);
       u_B = coeffs(p_idx,1:2);
       
       u_B = (u_B + 2*dot(u_A,n_hat)*n_hat - u_A) / 2;
       u_A = 2*dot(u_B,n_hat)*n_hat - u_B;
       
       coeffs(p_idx,1) = u_B(1);
       coeffs(p_idx,2) = u_B(2);
       coeffs(p_idx,3) = (coeffs(p_idx,3) + coeffs(m_idx,3)) / 2; %u_z =  u'_z

       coeffs(m_idx,1) = u_A(1);
       coeffs(m_idx,2) = u_A(2);
       coeffs(m_idx,3) = coeffs(p_idx,3);

       sep = sep+1;
   end
   
   
   % now the (m,m/2) symm
   if mod(i,2) == 0
      tar_idx = 2 + (i-1)*i/2 + i/2;
      coeffs(tar_idx,2) = ( coeffs(tar_idx,2) - sqrt(3)*coeffs(tar_idx,1)) / 2; % u_y = -sqrt(3) u_x
      coeffs(tar_idx,1) = -coeffs(tar_idx,2)/sqrt(3);  
   end
   
   % and the (m,0) symm %need to fix!
   zero_idx = 2 + (i-1)*i/2;
   coeffs(zero_idx,1) = 0.0; % u = (0,u_y)
   
end

% loop over fourier component
for i = 0:max_k
    for j = 0:max((i-1),0)
        r_max = 2;
        if (i == 0)
            % dont do any rotations for the [0, 0] component
           r_max = 0; 
        end
        % loop for 3 symm. rotated directions
        for r_idx = 0:r_max
            
            k1 = rot^(r_idx)*(i*sc_b1 + j*sc_b2);
            
            c = rot^(r_idx)*[coeffs(idx,1:2)]';
    
            % coeffs are for exp(1i*k*r), but we can use their trig
            % representations by multiplying by 2 and doing just 3 r_idx's
            disp_x = disp_x + c(1)         * 2*sin(k1(1)*atom_pos_list(:,1) + k1(2)*atom_pos_list(:,2));
            disp_y = disp_y + c(2)         * 2*sin(k1(1)*atom_pos_list(:,1) + k1(2)*atom_pos_list(:,2));
            disp_z = disp_z + coeffs(idx,3)* 2*cos(k1(1)*atom_pos_list(:,1) + k1(2)*atom_pos_list(:,2));

        end
        %fprintf("[%d, %d] coeffs(%d,:) = [%f, %f %f] \n",i,j,idx,coeffs(idx,1),coeffs(idx,2),coeffs(idx,3));
        %fprintf("[%d, %d] coeffs(%d) = %f (%f deg) \n",i,j,idx,norm(coeffs(idx,1:2)),angle(coeffs(idx,1) + 1j*coeffs(idx,2))*180/pi);

        idx = idx+1;

    end
end

% correct u_z to be at 0 at the AB location
disp_z = disp_z - max(max(disp_z));

% change sign depending on bot vs top layer
disp_x = u_sign*disp_x;
disp_y = u_sign*disp_y;
disp_z = u_sign*disp_z;

% compute how far sc_b1 is off the x-axis
% so we can rotate our basis into it 
% (coeffs assume that b1 is parallel to x)
b1_t = atan2(sc_b1(2),sc_b1(1));
b1_rot = [cos(b1_t) -sin(b1_t);
            sin(b1_t)  cos(b1_t)]; 

% save final results in correct basis!
atom_pos_shift(:,1) = b1_rot(1,1)*disp_x + b1_rot(1,2)*disp_y;
atom_pos_shift(:,2) = b1_rot(2,1)*disp_x + b1_rot(2,2)*disp_y;
atom_pos_shift(:,3) = disp_z;

end

