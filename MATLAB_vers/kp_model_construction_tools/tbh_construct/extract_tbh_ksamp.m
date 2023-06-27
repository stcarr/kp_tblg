function [] = extract_tbh_ksamp(m,knum)
   

    % M > N
    sc_m=m;
    sc_n=sc_m-1;

    fft_relax_onoff_factor=1.0;
    graphene_setup_supercell;
    graphene_setup_supercell_2;
    build_supercell_ham;

    know=sc_gamma;

    graphene_setup_supercell_relax_ver2;

    rot_angle = rot_theta/pi*180;


    graphene_onsite_e=-0.7876;

    Hmat_gamma = ham_supercell_N8_Koshino_relax_height_sparse( know/lattice_a, total_num, inplaneGAB_bottom_relax,inplaneGAB_top_relax, inplaneGAABB_bottom_relax,inplaneGAABB_top_relax,inter_relax ,graphene_onsite_e);

    % input variables:
    %
    % in Shiang's code, normally I use the length scale lattice_a = 1 !!!

    BZknum=knum;
    tmpkk=linspace(0,1,BZknum+1);
    tmpkk=tmpkk(1:BZknum);

    [allk1,allk2]=meshgrid(tmpkk,tmpkk);
    allkxy(:,1)=sc_b1(1)*allk1(:)+sc_b2(1)*allk2(:);
    allkxy(:,2)=sc_b1(2)*allk1(:)+sc_b2(2)*allk2(:);

    BZ_allham={};
    indcnt=1;
    for indbz=1:(BZknum^2)
        fprintf([num2str(100*indbz/BZknum^2,'%.0f') '%% done with k sampling \n']);
        know=sc_gamma;
        know(1)=allkxy(indbz,1);
        know(2)=allkxy(indbz,2);

        Hmat_tmp = ham_supercell_N8_Koshino_relax_height_sparse( know/lattice_a, total_num, inplaneGAB_bottom_relax,inplaneGAB_top_relax, inplaneGAABB_bottom_relax,inplaneGAABB_top_relax,inter_relax ,graphene_onsite_e);

        BZ_allham{indcnt}=Hmat_tmp;
        indcnt=indcnt+1;

    end
    
    save( ['reduced_dat_BZscan_',num2str(sc_m),'_',num2str(sc_n)], 'rot_angle', ...
        'sc_b1','sc_b2','pc_vec_b1','pc_vec_b2','pc_vec_rb1','pc_vec_rb2', ...
        'BZ_allham','allkxy','sc_all_points_shift','tot_num', ...
        'bot_A_atoms','bot_B_atoms','top_A_atoms','top_B_atoms');

end

