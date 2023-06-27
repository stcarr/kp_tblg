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
    
    tic
    save([tbh_data_dir '/' filename_here_tbh],'-v7.3');
    toc
    
    graphene_onsite_e=0.0;
    
    Hmat_gamma = ham_supercell_N8_Koshino_relax_height_sparse( know/lattice_a, total_num, inplaneGAB_bottom_relax,inplaneGAB_top_relax, inplaneGAABB_bottom_relax,inplaneGAABB_top_relax,inter_relax ,graphene_onsite_e);
    
    % input variables:
    %
    % in Shiang's code, normally I use the length scale lattice_a = 1 !!!
    
    BZknum = 4;
    %BZknum = 1;
    tmpkk=linspace(0,1,BZknum+1);
    tmpkk=tmpkk(1:BZknum);

    [allk1,allk2]=meshgrid(tmpkk,tmpkk);
    allkxy(:,1)=sc_b1(1)*allk1(:)+sc_b2(1)*allk2(:);
    allkxy(:,2)=sc_b1(2)*allk1(:)+sc_b2(2)*allk2(:);
   
    BZ_allham={};
    indcnt=1;
    for indbz=1:(BZknum^2)

        tic
        know=sc_gamma;
        know(1)=allkxy(indbz,1);
        know(2)=allkxy(indbz,2);

        Hmat_tmp = ham_supercell_N8_Koshino_relax_height_sparse( know/lattice_a, total_num, inplaneGAB_bottom_relax,inplaneGAB_top_relax, inplaneGAABB_bottom_relax,inplaneGAABB_top_relax,inter_relax ,graphene_onsite_e);

        BZ_allham{indcnt}=Hmat_tmp;
        indcnt=indcnt+1;
        fprintf([num2str(100*indbz/BZknum^2,'%.0f') '%% done with k sampling: ']);
        toc


    end


    filename_here = ['TwBLG_BZscan_',num2str(sc_m),'_',num2str(sc_n),'_' interlayer_model '-inter'];
    fprintf(['saving to file: ' filename_here '.mat \n']);
    save([tbh_data_dir '/' filename_here],'-v7.3');
    
    %'done'
end

%%

% plot(all_sc_m,all_rot_angle)


