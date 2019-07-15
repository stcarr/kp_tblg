%clear all
orig_dir = pwd;
addpath('../..')
addpath('../projection_tools/')
addpath('../../../data')

tar_folder = 'run_folder_8band';

% theta_list: angles to do projection
% just one:
theta_list = 1.10;

% multiple:
%theta_list = 0.6:.05:1.3;

% start from the first step
start_step = 1;

% generate the WF from the rewan process
rewan_wf = 1;

% calculate and plot band structure at finish or not
%calc_bands = 0; % no longer used, use the final *_rewan_rs.dat file instead


for idx = 1:size(theta_list,2)
    tar_theta = theta_list(idx);

    numk_cs = 15;
    numxy_cs = 31;

    numk = 16;
    numxy = 15;
    
    cutoff_scale = 4;
    
    mkdir(tar_folder);

    % STEP 1: do CS projections
    if (start_step <= 1) % step 1
        eight_cs_filename = eightWANcs(tar_theta,tar_folder,numk_cs,numxy_cs);
        %five_cs_top_filename = fiveWANcs(tar_theta,tar_folder,numk_cs,numxy_cs,0);
        %five_cs_bot_filename = fiveWANcs(tar_theta,tar_folder,numk_cs,numxy_cs,1);
        save([tar_folder '/wanFilenames_' strrep(num2str(tar_theta,'%.2f'),'.','p') '.mat'], ...
                'eight_cs_filename');
    else
        load([tar_folder '/wanFilenames_' strrep(num2str(tar_theta,'%.2f'),'.','p') '.mat']);
    end
    % STEPS 2 & 3: do realistic projection (WF_0)
    if (start_step <= 3) % step 2 & 3
        trio_cs_top_filename = trioWANcs(tar_theta,tar_folder,numk_cs,numxy_cs,0);
        trio_cs_bot_filename = trioWANcs(tar_theta,tar_folder,numk_cs,numxy_cs,1);
        [eight_wf_filename,eight_h_filename] = eightWANfull_trios_together(trio_cs_top_filename,trio_cs_bot_filename,eight_cs_filename,tar_theta,tar_folder,numk,numxy);
        %[five_top_wf_filename, five_top_h_filename] = fiveWANfull(trio_cs_top_filename,five_cs_top_filename,tar_theta,tar_folder,numk,numxy,0);
        %[five_bot_wf_filename, five_bot_h_filename] = fiveWANfull(trio_cs_bot_filename,five_cs_bot_filename,tar_theta,tar_folder,numk,numxy,1);
        save([tar_folder '/wanFilenames_' strrep(num2str(tar_theta,'%.2f'),'.','p') '.mat'], ...
                'eight_cs_filename', ...
                'trio_cs_top_filename','trio_cs_bot_filename',...
                'eight_wf_filename','eight_h_filename');
    end
    % STEPS 4 & 5: do unweighted realistic projection (WF_1)
    if (start_step <= 5) % step 4 & 5
        [eight_wf_new_filename, eight_h_new_filename] = eightWANfull_fromWF(eight_wf_filename,tar_theta,tar_folder,numk);
        %[five_top_wf_new_filename, five_top_h_new_filename] = fiveWANfull_fromWF(five_top_wf_filename,tar_theta,tar_folder,numk,numxy,0);
        %[five_bot_wf_new_filename, five_bot_h_new_filename] = fiveWANfull_fromWF(five_bot_wf_filename,tar_theta,tar_folder,numk,numxy,1);
        [hmat_filename] = generate_hmat_datafile(eight_h_new_filename, tar_folder, tar_theta);
        
        save([tar_folder '/wanFilenames_' strrep(num2str(tar_theta,'%.2f'),'.','p') '.mat'], ...
            'eight_cs_filename', ...
            'trio_cs_top_filename','trio_cs_bot_filename',...
            'eight_wf_filename','eight_h_filename', ...
            'eight_wf_new_filename','eight_h_new_filename',...
            'hmat_filename');
    end
    
    % REWAN step:
    if (start_step <= 6)
        
        % define the sigma_k smearing used for fixing the Gamma-point
        % leakage of the flat-bands near the magic angle.
        if (tar_theta > 1.09)
            k_radius = 0.25;
        else
            k_radius = (1.1/tar_theta)*0.25;
        end
        
        % rewannier, calculate updated Unitary transformation
        [rewan_filename, postproc_filename] = rewannier_aapm(hmat_filename,tar_folder,tar_theta,numk,k_radius,0);
        
        % calc rewan wavefunctions, if turned on
        if (rewan_wf == 1)
            eight_wf_rewan_filename = eightWAN_WF_modifyU(eight_h_new_filename,eight_wf_new_filename,postproc_filename,tar_theta,tar_folder,numk);
        end

        % symmetrize final 8-band hamiltonain
        rewan_rs_filename = real_space_sym(rewan_filename,tar_folder,tar_theta,1);
        
        save([tar_folder '/wanFilenames_' strrep(num2str(tar_theta,'%.2f'),'.','p') '.mat'], ...
            'eight_cs_filename', ...
            'trio_cs_top_filename','trio_cs_bot_filename',...
            'eight_wf_filename','eight_h_filename', ...
            'eight_wf_new_filename','eight_h_new_filename',...
            'hmat_filename', ...
            'rewan_filename','postproc_filename', ...
            'rewan_rs_filename');
        
    end

    
    %old bandstructure code, use the rewan_rs_filename instead    
    %{
    if (calc_bands == 1)
        bandcalcTBHsym_eightband(eight_h_new_filename,eight_wf_new_filename,tar_theta,cutoff_scale);    
    end
    %}
end

