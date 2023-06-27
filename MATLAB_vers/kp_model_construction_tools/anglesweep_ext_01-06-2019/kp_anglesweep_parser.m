
clear all;

sc_m_list = [7:1:30,32:2:40,45:5:80,100,120,150,180];
%sc_m_list = [7:1:35,40:5:100];


%tar_folder = 'no_relax/data/kp_data/';
tar_folder = 'dft_full_relax_larger_cutoff/kp_data/';
%tar_folder = 'dft_flat_relax/data/kp_data/';

for m_idx = 1:length(sc_m_list)
   sc_m = sc_m_list(m_idx);
   sc_n = sc_m-1;
   filename = [tar_folder 'TwBLG_EffKP_' num2str(sc_m) '_' num2str(sc_n) '_dft-inter.mat'];
   load(filename)
   thetas(m_idx) = acos((sc_m^2+sc_n^2+4*sc_m*sc_n)/2/(sc_m^2+sc_m*sc_n+sc_n^2));
   inter_kp(m_idx,:,:,:,:) = All_Eff_inter_ext;
   intra_bot_kp(m_idx,:,:,:,:) = All_Eff_intra_bot_ext;
   intra_top_kp(m_idx,:,:,:,:) = All_Eff_intra_top_ext;
   inter_shells = All_Eff_inter_shell_indices;
   intra_shells = All_Eff_intra_shell_indices;
   
end
%%
save('dft_full_relax_data_02-04-2019','thetas','inter_kp','intra_top_kp','intra_bot_kp','inter_shells','intra_shells');
