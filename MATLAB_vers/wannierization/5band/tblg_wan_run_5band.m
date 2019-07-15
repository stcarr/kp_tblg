%clear all
orig_dir = pwd;
addpath('../..')
addpath('../projection_tools/')
addpath('../../../data')

tar_folder = 'run_folder_5band';

% theta_list: angles to do projection
% just one:
theta_list = 1.10;

% multiple:
%theta_list = 0.85:.05:1.3;

% start from the first step
start_step = 1;

% generate the WF from the rewan process
rewan_wf = 1;

clear frames
for idx = 1:size(theta_list,2)
    tar_theta = theta_list(idx);

    % k-space sampling for chiral symmetric model
    numk_cs = 15;
    % real-space grid size (per moire-cell) for CS model
    numxy_cs = 31;
    
    
    % k-space sampling for realistic model
    numk = 16;
    % real-space grid size (per moire-cell) for realistic model
    numxy = 15;
    
    % truncation length, in moire lengths
    cutoff_scale = 4;
    
    mkdir(tar_folder);

    % STEP 1: do CS projections
    if (start_step <= 1) % step 1
        %eight_cs_filename = eightWANcs(tar_theta,tar_folder,numk_cs,numxy_cs);
        five_cs_top_filename = fiveWANcs(tar_theta,tar_folder,numk_cs,numxy_cs,0);
        five_cs_bot_filename = fiveWANcs(tar_theta,tar_folder,numk_cs,numxy_cs,1);
        save([tar_folder '/wanFilenames_' strrep(num2str(tar_theta,'%.2f'),'.','p') '.mat'], ...
                'five_cs_top_filename', 'five_cs_bot_filename');
    else
        load([tar_folder '/wanFilenames_' strrep(num2str(tar_theta,'%.2f'),'.','p') '.mat']);
    end
    
    % STEPS 2 & 3: do realistic projection (WF_0)
    if (start_step <= 3) % step 2 & 3
        trio_cs_top_filename = trioWANcs(tar_theta,tar_folder,numk_cs,numxy_cs,0);
        trio_cs_bot_filename = trioWANcs(tar_theta,tar_folder,numk_cs,numxy_cs,1);
        %[eight_wf_filename,eight_h_filename] = eightWANfull_trios_together(trio_cs_top_filename,trio_cs_bot_filename,eight_cs_filename,tar_theta,tar_folder,numk,numxy);
        [five_top_wf_filename, five_top_h_filename] = fiveWANfull(trio_cs_top_filename,five_cs_top_filename,tar_theta,tar_folder,numk,numxy,0);
        [five_bot_wf_filename, five_bot_h_filename] = fiveWANfull(trio_cs_bot_filename,five_cs_bot_filename,tar_theta,tar_folder,numk,numxy,1);
        save([tar_folder '/wanFilenames_' strrep(num2str(tar_theta,'%.2f'),'.','p') '.mat'], ...
                'five_cs_top_filename', 'five_cs_bot_filename', ...
                'trio_cs_top_filename','trio_cs_bot_filename',...
                'five_top_wf_filename','five_top_h_filename',...
                'five_bot_wf_filename','five_bot_h_filename');
    end
    
    % STEPS 4 & 5: do unweighted realistic projection (WF_1)
    if (start_step <= 5) % step 4 & 5
        %[eight_wf_new_filename, eight_h_new_filename] = eightWANfull_fromWF(eight_wf_filename,tar_theta,tar_folder,numk,numxy);
        [five_top_wf_new_filename, five_top_h_new_filename] = fiveWANfull_fromWF(five_top_wf_filename,tar_theta,tar_folder,numk,numxy,0);
        [five_bot_wf_new_filename, five_bot_h_new_filename] = fiveWANfull_fromWF(five_bot_wf_filename,tar_theta,tar_folder,numk,numxy,1);
        
        [hmat_top_filename] = generate_hmat_datafile_5band(five_top_h_new_filename, tar_folder, tar_theta, cutoff_scale, 0);
        [hmat_bot_filename] = generate_hmat_datafile_5band(five_bot_h_new_filename, tar_folder, tar_theta, cutoff_scale, 1);

        save([tar_folder '/wanFilenames_' strrep(num2str(tar_theta,'%.2f'),'.','p') '.mat'], ...
                'five_cs_top_filename', 'five_cs_bot_filename', ...
                'trio_cs_top_filename','trio_cs_bot_filename',...
                'five_top_wf_filename','five_top_h_filename',...
                'five_bot_wf_filename','five_bot_h_filename',...
                'five_top_wf_new_filename','five_top_h_new_filename',...
                'five_bot_wf_new_filename','five_bot_h_new_filename',...
                'hmat_top_filename','hmat_bot_filename');
    end
    
    % REWAN step:
    if (start_step <= 6)
        
        % define the sigma_k smearing used for fixing the Gamma-point
        % leakage of the flat-bands near the magic angle.
        if (tar_theta > 1.09)
            k_radius = .4;
        else
            k_radius = (1.1/tar_theta)*0.4;
        end
                        
        % rewannier, calculate updated Unitary transformation
        [rewan_top_filename, postproc_top_filename]= rewannier_aapm_5band(hmat_top_filename,tar_folder,tar_theta,numk,k_radius,0,0);
        [rewan_bot_filename, postproc_bot_filename]= rewannier_aapm_5band(hmat_bot_filename,tar_folder,tar_theta,numk,k_radius,1,0);
        
        % calc rewan wavefunctions, if turned on
        if (rewan_wf == 1)
            five_top_wf_rewan_filename = fiveWAN_WF_modifyU(five_top_h_new_filename,five_top_wf_new_filename,postproc_top_filename,tar_theta,tar_folder,0,numk);
            five_bot_wf_rewan_filename = fiveWAN_WF_modifyU(five_top_h_new_filename,five_top_wf_new_filename,postproc_bot_filename,tar_theta,tar_folder,1,numk);

        end

        % symmetrize final 5-band hamiltonains
        rewan_top_rs_filename = real_space_sym_5band(rewan_top_filename,tar_folder,tar_theta,0,1);
        rewan_bot_rs_filename = real_space_sym_5band(rewan_bot_filename,tar_folder,tar_theta,1,1);

        save([tar_folder '/wanFilenames_' strrep(num2str(tar_theta,'%.2f'),'.','p') '.mat'], ...
            'five_cs_top_filename', 'five_cs_bot_filename', ...
            'trio_cs_top_filename','trio_cs_bot_filename',...
            'five_top_wf_filename','five_top_h_filename',...
            'five_bot_wf_filename','five_bot_h_filename',...
            'five_top_wf_new_filename','five_top_h_new_filename',...
            'five_bot_wf_new_filename','five_bot_h_new_filename',...
            'hmat_top_filename','hmat_bot_filename', ...
            'rewan_top_filename','rewan_bot_filename', ...
            'postproc_top_filename','postproc_bot_filename',...
            'rewan_top_rs_filename','rewan_bot_rs_filename');
        
    end
      
    
end

