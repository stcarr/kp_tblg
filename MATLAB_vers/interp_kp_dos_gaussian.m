function [dos_sweep, idos_sweep, E_list, half_filling_hole_E] = interp_kp_dos_gaussian(theta_list, sweep_vals, sweep_kpts)

    tar_theta_list = theta_list;

    % number of extra bands to include
    b_size = 4;


    for t_idx = 1:length(tar_theta_list)
        fprintf("Starting theta %d/%d: ",t_idx,length(tar_theta_list));
        tic
        
        tri_idx = 1;

        all_kpts1 = sweep_kpts{t_idx};
        allbands1 = sweep_vals{t_idx};

        kpoints = all_kpts1(:,[2 1]);
        eig_vals = allbands1';

        nk = sqrt(size(eig_vals,1));
        nb = size(eig_vals,2);
        max_E = 0.3;
        dE = max_E/6000;
        
        E_list = [-max_E:dE:max_E];
        dos = zeros(length(E_list),1);
        
        sig = .005;

        for tar_b = (nb/2)-b_size:(nb/2+1)+b_size
            
            for i = 1:nk
                for j = 1:nk
                    Eh = eig_vals((i-1)*nk + j,tar_b);
                    dos = dos + exp(-(E_list'-Eh).^2/(2*sig^2))/(sqrt(2*pi*sig^2))/nk^2;
                end
            end
            
        end

        %dos_sweep{t_idx} = dos;

        idos = zeros(size(dos));
        for x = 1:length(E_list)
            if (x > 1)
                idos(x) = trapz(E_list(1:x),dos(1:x));
            end
        end

        tot_bands = 2*(b_size+1);
        tot_bands = 4*tot_bands; % 2 for valley, 2 for spin
        
        alpha = 2.47;
        sc_alpha = alpha/(2*sind(tar_theta_list(t_idx)/2));
        sc_area = sc_alpha^2*sind(60)*1e-2; %area in nm^2
        n0 = 1/sc_area;
        
        idos_rescale = tot_bands/idos(end);
        dos_rescale = idos_rescale*n0

        idos(:) = idos_rescale*(idos(:) - 0.5*idos(end));
        [val, idx] = min(abs(idos - (-2)));

        idos_sweep{t_idx} = idos;
        dos_sweep{t_idx} = dos_rescale*dos;
        half_filling_hole_E(t_idx) = E_list(idx);
    
        toc
    end
    
end