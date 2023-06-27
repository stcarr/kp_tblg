%clear all;

% set path / environment variables
home_dir = '/home/stc/devspace/codes/shiang_kp_relaxed_tblg/kp_construction_tools';

%base_dir = '/home/stc/devspace/codes/shiang_kp_relaxed_tblg/kp_runs/';
%run_name = 'first_pass';
%main_dir = [base_dir run_name '/'];
fprintf('Running kp_param_extract_job.m for %s \n',main_dir);

global_var_file = [main_dir '/global_vars.m'];

data_dir = [main_dir '/data/'];
tbh_data_dir = [data_dir 'tbh_data/'];
kp_data_dir = [data_dir 'kp_data/'];
bands_data_dir = [data_dir 'bands/'];

clearvars -except home_dir base_dir run_name main_dir global_var_file data_dir tbh_data_dir kp_data_dir bands_data_dir

% make the subfolders if they do not exist

cd(main_dir);
mkdir('data');

cd(data_dir);
mkdir('bands');
mkdir('kp_data');
mkdir('tbh_data');

% begin the job!

cd(home_dir);
cd kp_fit
%KP_construct_eff_theory   % NN only, no k+- dependence
%KP_construct_eff_theory_relax_ext
KP_construct_eff_theory_relax_miller
clearvars -except home_dir base_dir run_name main_dir global_var_file data_dir tbh_data_dir kp_data_dir bands_data_dir
fprintf('done with KP fitting \n');

%{
cd(home_dir);
cd kp_bands
%kp_band_calc_lambda_nn    % NN only, no k+- dependence
TwBLG_relax_KP_kdep_ver2
clearvars -except home_dir base_dir run_name main_dir global_var_file data_dir tbh_data_dir kp_data_dir bands_data_dir
fprintf('done with KP bands \n');
%}

cd(home_dir);

% done
