
run(global_var_file);

for indm=1:length(all_sc_m)
    
%     % M > N
%     % preparation step
%     % M > N
      sc_m=all_sc_m(indm);
      sc_n=sc_m-1;
      [sc_m,sc_n]
%     
%     fft_relax_onoff_factor=1.0;
%     graphene_setup_supercell;
%     graphene_setup_supercell_2;
%     build_supercell_ham;
%     
%     know=sc_gamma;
%     
%     graphene_setup_supercell_relax_ver2;
%     
%     all_rot_angle(indm)=rot_theta/pi*180;
%     
%     
%     graphene_onsite_e=-0.7876;
%     
%     Hmat_gamma = ham_supercell_N8_Koshino_relax_height( know/lattice_a, total_num, ,inplaneGAB_top_relax, inplaneGAABB_bottom_relax,inplaneGAABB_top_relax,inter_relax ,graphene_onsite_e);
%     
%     % input variables:
%     %
%     % in Shiang's code, normally I use the length scale lattice_a = 1 !!!
%     
    
    tbh_filename_here = ['TwBLG_TBH-dat_',num2str(sc_m),'_',num2str(sc_n),'_' interlayer_model '-inter'];
    bands_filename_here = ['TwBLG_tbh-bands_',num2str(sc_m),'_',num2str(sc_n),'_' interlayer_model '-inter'];
    
    fprintf(['loading file: ' tbh_filename_here '.mat \n']);
    load([tbh_data_dir '/' tbh_filename_here]);    

    %load(['TwBLG_EffKP_relax_symtest1_',num2str(sc_m),'_',num2str(sc_n),'_xy']);
    
    %'done loading'
    
    %knum=20;
    % scan_klist=[[sc_gamma,0];[sc_mpoint,0];[sc_kpoint,0];[sc_gamma,0]];
    
    scan_klist=[[sc_kpoint,0];[sc_gamma,0];[sc_mpoint,0];[sc_kpoint,0]];
    
    [ all_kpts, scale_axis] = generate_k_line( knum, scan_klist );
    knum_tot=size(all_kpts);
    knum_tot=knum_tot(1);
    
    allbands=zeros(tot_num,knum_tot);
    graphene_onsite_e=-0.7876;
    
    %parfor
    for indk=1:knum_tot
        tic
        %fprintf([num2str(100*indk/knum_tot,'%.0f') '%% done with bandstructure k sampling \n']);

        know=all_kpts(indk,1:2);
        
        Hmat = ham_supercell_N8_Koshino_relax_height( know/lattice_a, total_num, inplaneGAB_bottom_relax,inplaneGAB_top_relax, inplaneGAABB_bottom_relax,inplaneGAABB_top_relax,inter_relax ,graphene_onsite_e);
        
        
        [bands]=eig(Hmat);
        
        allbands(:,indk)=sort(real(bands),'ascend');
 	  
        fprintf([num2str(100*indk/knum_tot,'%.0f') '%% done with bandstructure k sampling \n']);
        toc 
end
    
    fprintf(['saving file: ' bands_filename_here '.mat \n']);
    save([bands_data_dir '/' bands_filename_here],'scale_axis','all_kpts','allbands');    
    
end

%%

% plot(all_sc_m,all_rot_angle)



