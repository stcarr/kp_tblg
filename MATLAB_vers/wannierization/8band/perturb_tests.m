clear all;
addpath('../..')
addpath('../projection_tools/')
addpath('../../../data')

% folder where we saved the wannier projection
tar_folder = 'run_folder_8band';
tar_theta = 1.1;
numk = 16;

% file of the output of step 5 (before REWAN)
h_new_filename = [tar_folder,'/eightWANnew_H_' strrep(num2str(tar_theta,'%.2f'),'.','p') '_' num2str(numk) 'k.mat'];

% perturbation strength, in eV
delta_i = (1e-3)*30; % 30 meV

% we do one at a time, second comment has values from paper figures
disp_str        = 0;           %(1e-3)*150;
sub_str_sym     = delta_i;%delta_i;      %(1e-3)*30;
sub_str_asym    = delta_i;%delta_i;            %(1e-3)*30;

% get the band structure of the perturbed k-dot-p theory
[kp_bands, scale_axis] = perturb_band_kp_calc(tar_theta,31,disp_str, sub_str_sym, sub_str_asym);

% get the perturbed reduced Hamiltonian (first order in pert. theory)
[perturb_filename] = eightWAN_perturb_calc(h_new_filename,tar_theta,tar_folder,numk,disp_str, sub_str_sym, sub_str_asym);

% NOTE!! We do not need to pass perturb_filename through the REWAN function
% instead we will just use the rewan gauge transformation U and add it
% ontop of the perturbed hamiltonian.
% This avoids issues with the REWAN process failing due to the perturbation

%% get realspace H perturbation
% this is for viewing the size of the hopping terms
% introduced by the perturbation(s).

postproc_filename = [tar_folder '/hmat_8band_' strrep(num2str(tar_theta,'%.2f'),'.','p') '_rewan_postproc.mat'];

% we can get the pertubation before or after REWAN
% after REWAN:
[H_perturb, hex_R, grid_to_index] = Generate_8band_TBH_perturb_rewan(perturb_filename,postproc_filename);
% before REWAN:
[H_perturb_pre, hex_R, grid_to_index] = Generate_8band_TBH_perturb(perturb_filename);


%% plot H_perturb 

    grid_size = size(grid_to_index,1);
    H_map = zeros(grid_size*9);
    H_pos_x = zeros(grid_size*9);
    H_pos_y = zeros(grid_size*9);
    
    %H_tar = H_perturb_pre; % pre-REWAN H
    H_tar = H_perturb; % post-REWAN H

    % plots the Hamiltonian

    %splot(3) = subplot(1,7,3);
    clf
    hold on

    ax_c = [31.5 21.5];
    
    MMM = 20;

    for ind1 = 1:2*MMM+1
        for ind2 = 1:2*MMM+1
            idx_here = grid_to_index(ind1,ind2);
            i_start = 9*(ind1-1);
            j_start = 9*(ind2-1);
            if (idx_here ~= 0)
                H_map(i_start+[1:8],j_start+[1:8]) = fliplr(abs(H_tar(:,:,idx_here)));
            end
            for idx = 1:9
                H_pos_x(i_start+[1:9],j_start+idx) = ind1 + (ind2-1)*0.5 + linspace(0,1,9);
                H_pos_y(i_start+idx,j_start+[1:9]) = ind2 + linspace(0,1,9);
            end
            if (abs(ind1-MMM) + abs(ind2-MMM) < 15)
                plot3(ind1 + (ind2-1)*0.5+[0 1]-ax_c(1),[ind2 ind2]-ax_c(2),[5 5],'-w')
                plot3(ind1+(ind2-1)*0.5+[0 0]-ax_c(1),[ind2 ind2+1]-ax_c(2),[5 5],'-w')
            end


        end
    end
    
    H_map = (H_map/delta_i);

    %c_min = -10;
    %c_max = -1;
    c_min = 0;
    c_max = max(abs(H_map(:)))*1;

    caxis([c_min c_max])
    nc = size(colormap,1);
    c_axis = linspace(c_min,c_max,nc);

    %scatter(H_pos_x(:),H_pos_y(:))
    surf(H_pos_x - ax_c(1), H_pos_y - ax_c(2),(H_map),'EdgeAlpha',.4)
    %imagesc(H_pos_x(:),H_pos_y(:),interp1(c_axis,colormap,log(H_map(:))) ) 

    axis equal
    ax_w = 3;
    axis([-ax_w ax_w -ax_w ax_w -inf inf])
    view(2)
    %cbar = colorbar('north');
    colorbar
    colormap hot
    %set(cbar,'pos',[.41 .1 .215 .02])
    title('$\theta = 1.1^\circ$, Sublattice $|\Delta H^{ij}_R|$')
    xlabel('$x$ ($\lambda_\theta$)')
    ylabel('$y$ ($\lambda_\theta$)')




%% get bands
% for comparing band-structure effects

% can't use the symmetrized from of the perturbed Hamiltonian!
[~, scale_axis, bands_perturb] = Generate_8band_TBH_sym_ver4(perturb_filename);

% save the symmetrized unperturbed Hamiltonian for (optional) comparison.
[bands_orig_sym, scale_axis, bands_orig] = Generate_8band_TBH_sym_ver4(h_new_filename);
%% check bands
clf
hold on

% unperturbed original 8-band Hamiltonian
%plot(scale_axis,bands_orig_sym,'--k')

% perturbed kp band structure
plot(scale_axis,kp_bands,'k')

% perturbed 8-band Hamiltonian band structure (should agree at low energy)
plot(scale_axis,bands_perturb,'r')

ax_m = 0.25;
axis([0 1 -ax_m ax_m])


%max(abs(Heff_perturb(:)))

