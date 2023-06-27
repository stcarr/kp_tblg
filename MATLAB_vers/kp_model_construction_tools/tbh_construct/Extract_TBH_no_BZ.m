run(global_var_file);

for indm=1:length(all_sc_m)
    
    % M > N
    % preparation step
    % M > N
    sc_m=all_sc_m(indm);
    sc_n=sc_m-1;
    
    tic
    graphene_setup_supercell;
    fprintf("done with graphene_setup_supercell.m : ");
    toc
    
    tic
    graphene_setup_supercell_2;
    fprintf("done with graphene_setup_supercell_2.m :");
    toc
    
    tic
    build_supercell_ham;
    fprintf("done with build_supercell_ham.m :");
    toc
    
    know=sc_gamma;
    
    tic
    graphene_setup_supercell_relax_ver2;
    fprintf("done with graphene_setup_supercell_relax_ver2.m : ");
    toc
    
    all_rot_angle(indm)=rot_theta/pi*180;
    
    filename_here_tbh = ['TwBLG_TBH-dat_',num2str(sc_m),'_',num2str(sc_n),'_' interlayer_model '-inter'];
    fprintf(['saving to file: ' filename_here_tbh '.mat : ']);
    
    save([tbh_data_dir '/' filename_here_tbh]);

end