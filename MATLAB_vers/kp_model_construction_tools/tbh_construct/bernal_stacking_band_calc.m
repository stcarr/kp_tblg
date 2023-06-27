clear all;
lattice_a=1.42*sqrt(3);

layer_d=[0,0,3.35];

% lattice_a*
a1=lattice_a*[sqrt(3)/2,-1/2]';
a2=lattice_a*[sqrt(3)/2,+1/2]';

A = [a1 a2];
G = 2*pi*inv(A)';

b2 = G(:,1)+G(:,2);
b1 = G(:,1);

K_pt = 1/3*(2*G(:,1)+G(:,2));

b_shift = (a1+a2)/3;

pos_list = zeros(3,4);

% B orbital positions
pos_list(1:2,2) = b_shift;
pos_list(1:2,4) = b_shift;

% bernal stacking
pos_list(:,3) = pos_list(:,3) + [b_shift;layer_d(3)];
pos_list(:,4) = pos_list(:,4) + [b_shift;layer_d(3)];

sc_grid = 5;

% graphene TBH from DFT with strain corrections
t_list = 1.2*[-2.822, 0*0.254, -0.180];
r_list = lattice_a/sqrt(3)*[1,sqrt(3),2];

k_samp = 20;

dmesh = linspace(0,1,3*k_samp+1);
dmesh = dmesh(1:end-1);
[mesh_x,mesh_y] = meshgrid(dmesh,dmesh);

k_list(1,:) = b1(1)*mesh_x(:) + b2(1)*mesh_y(:);
k_list(2,:) = b1(2)*mesh_x(:) + b2(2)*mesh_y(:);

for k_idx = 1:size(k_list,2)

    kh = k_list(:,k_idx);
    H = zeros(4,4);
            
    % interlayer terms
    for i = -sc_grid:sc_grid
        for j = -sc_grid:sc_grid
            R = i*a1 + j*a2;
            for l = 1%:2
                for o_from = 1:2
                    idx_f = (l-1)*2+o_from;
                    pos_f = pos_list(:,idx_f);
                    for o_to = 1:2
                        idx_t = (l)*2+o_to;
                        pos_t = pos_list(:,idx_t);

                        dr = [R;0] + pos_t - pos_f;
                        t_h = dft_interlayer_coupling(dr', 0, 0, lattice_a);

                        phase_h = exp(1j*dot(kh,dr(1:2)));
                        H(idx_t,idx_f) = H(idx_t,idx_f) + t_h*phase_h;

                    end

                end
            end
        end
    end
    H = H + H'; 
            
            % monolayer terms
    for i = -sc_grid:sc_grid
        for j = -sc_grid:sc_grid
            R = i*a1 + j*a2;
            for l = 1:2
                for o_from = 1:2
                    idx_f = (l-1)*2+o_from;
                    pos_f = pos_list(1:2,idx_f);
                    for o_to = 1:2
                        idx_t = (l-1)*2+o_to;
                        pos_t = pos_list(1:2,idx_t);

                        dr = R + pos_t - pos_f;
                        rl = norm(dr);

                        t_h = 0;
                        for idx = 1:3
                            if abs(rl - r_list(idx)) < 1e-4
                               t_h = t_list(idx);
                            end
                        end

                        phase_h = exp(1j*dot(kh,dr(1:2)));
                        H(idx_t,idx_f) = H(idx_t,idx_f) + t_h*phase_h;

                    end

                end
            end
                       

        end
    end

    bands(k_idx,:) = sort(real(eigs(H,4)));
    
end

%%
E_f = 0*-0.9144; %  to move Fermi energy back to 0 eV
[dos_sweep, ~, E_list, ~] = interp_kp_dos([1], {bands'-E_f}, {k_list'});

%%

clf

subplot(2,1,1)
plot(bands,'k')
%axis([80 82 -2 2])
title('Bands (BZ sampling)')
ylabel('Energy (eV)')
xlabel('k index')

dos = dos_sweep{1};

subplot(2,1,2)
plot(E_list,dos,'k')
%axis([-0.5 0.5 0 1.5])
title('DOS')

xlabel('Energy (eV)')
ylabel('DoS (states per (eV nm^2))')
