y_pos = 30;

sc_t1 = 100*[1 0]';
sc_t2 = 100*[1/2 sqrt(3)/2]';

rot_theta = (pi/180)*0.5;

A = [sc_t1 sc_t2];
G = 2*pi*inv(A);
sc_b1 = G(1,:);
sc_b2 = G(2,:);

sc_centerline = (sc_t1 + sc_t2);
fix_theta = angle(sc_centerline(1) + 1j*sc_centerline(2));
R = [cos(fix_theta) sin(fix_theta); -sin(fix_theta) cos(fix_theta)];
rot_centerline = R*sc_centerline;
rot_b1 = R*sc_b1(1:2)';
rot_b2 = R*sc_b2(1:2)';

sym_test_pos(1,:) = [rot_centerline(1)/2+6 y_pos];
sym_test_pos(2,:) = [rot_centerline(1)/2+6 -y_pos];
sym_test_pos(3,:) = [rot_centerline(1)/3 y_pos];
sym_test_pos(4,:) = [rot_centerline(1)/3 -y_pos];
sym_test_pos(5,:) = [2*rot_centerline(1)/3 y_pos];
sym_test_pos(6,:) = [2*rot_centerline(1)/3 -y_pos];
sym_test_pos(7,:) = [rot_centerline(1)/3 0]; % AB stacking
sym_test_pos(8,:) = [0 0]; % AA stacking

clf
hold on
plot(sym_test_pos(:,1), sym_test_pos(:,2),'.r');
rot_t1 = R*sc_t1;
rot_t2 = R*sc_t2;
plot([0 rot_t1(1)],[0 rot_t1(2)],'--k')
plot([0 rot_t2(1)],[0 rot_t2(2)],'--k')
plot(rot_t2(1)+[0 rot_t1(1)],rot_t2(2)+[0 rot_t1(2)],'--k')
plot(rot_t1(1)+[0 rot_t2(1)],rot_t1(2)+[0 rot_t2(2)],'--k')
axis equal

sym_test_disp = Relax_model_Stephen_ver3(rot_theta,rot_b1,rot_b2,sc_t1,sc_t2,0,sym_test_pos);

for l = 1:3
    error_table(l,1) = sym_test_disp(2*l - 1,1) + sym_test_disp(2*l,1);
    error_table(l,2) = sym_test_disp(2*l - 1,2) - sym_test_disp(2*l,2);
    error_table(l,3) = sym_test_disp(2*l - 1,3) - sym_test_disp(2*l,3);
end
error_table
%max(abs(sym_test_disp(:,3)))*2

%%

% filepath is where the relaxation data is located!
%filepath = './sweep_coeffs2/';
filepath = './0p1_30p0_k25_08-24-2018/';
% load the fourier coeff data
theta_list = dlmread([filepath 'thetas.txt']);
coeffs_x = dlmread([filepath 'coeffs_x.txt']);
coeffs_y = dlmread([filepath 'coeffs_y.txt']);
coeffs_z = dlmread([filepath 'coeffs_z.txt']);

% get the distance between theta samples
d_theta = theta_list(2) - theta_list(1);

% convert to degrees
tar_theta = rot_theta*180/pi;

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

b1_t = atan2(rot_b1(2),rot_b1(1));
b1_rot = [cos(b1_t) -sin(b1_t);
            sin(b1_t)  cos(b1_t)]; 
coeffs_xh = coeffs(:,1);
coeffs_yh = coeffs(:,2);
coeffs(:,1) = b1_rot(1,1)*coeffs_xh + b1_rot(1,2)*coeffs_yh;
coeffs(:,2) = b1_rot(2,1)*coeffs_xh + b1_rot(2,2)*coeffs_yh;

max_k = 25;
%max_k = 2;
idx = 1;

atom_pos_list = sym_test_pos;
N = length(atom_pos_list);

% set up some temporary vectors for storing displacements (might need to
% rotate them into a different basis at end)
disp_x = zeros(N,1);
disp_y = zeros(N,1);
disp_z = zeros(N,1);

rot = [cosd(60) -sind(60);
        sind(60)  cosd(60)];
 
% manually clean up the z-component symmetry!
for i = 0:max_k
   pivot = i/2;

   if mod(i,2) == 0
       sep = 1.0;
   else
       sep = 0.5;
   end
   
   while (pivot + sep) < i
       tar_j_p = pivot + sep;
       tar_j_m = pivot - sep;
       
       p_idx = 1 + (i-1)*i/2 + tar_j_p+1;
       m_idx = 1 + (i-1)*i/2 + tar_j_m+1;
       coeffs(p_idx,3) = (coeffs(p_idx,3) + coeffs(m_idx,3)) / 2;
       coeffs(m_idx,3) = coeffs(p_idx,3);
       sep = sep+1;
   end
   
end
clear fourier_comps
% loop over fourier componenet
for i = 0:max_k
    for j = 0:max((i-1),0)
        r_max = 2;
        if (i == 0)
            % dont do any rotations for the [0, 0] component
           r_max = 0; 
        end
        % loop for 3 symm. rotated directions
        for r_idx = 0:r_max
            
            k1 = rot^(r_idx)*(i*rot_b1 + j*rot_b2);
            fourier_comps(idx,:) = [i j];
            
            c = rot^(r_idx)*[coeffs(idx,1:2)]';
    
            % coeffs are for the exp(1i*k*r), but we can use their trig
            % representations by multiplying by 2 and doing just 3 r_idx's
            disp_x = disp_x + c(1)         * 2*sin(k1(1)*atom_pos_list(:,1) + k1(2)*atom_pos_list(:,2));
            disp_y = disp_y + c(2)         * 2*sin(k1(1)*atom_pos_list(:,1) + k1(2)*atom_pos_list(:,2));
            disp_z = disp_z + coeffs(idx,3)* 2*cos(k1(1)*atom_pos_list(:,1) + k1(2)*atom_pos_list(:,2));

        end
        
        idx = idx+1;

    end
end

sym_test_disp = [disp_x disp_y disp_z];
for l = 1:3
    error_table(l,1) = sym_test_disp(2*l - 1,1) + sym_test_disp(2*l,1);
    error_table(l,2) = sym_test_disp(2*l - 1,2) - sym_test_disp(2*l,2);
    error_table(l,3) = sym_test_disp(2*l - 1,3) - sym_test_disp(2*l,3);
end
error_table
%max(abs(sym_test_disp(:,3)))*2