clear all;

% 16:35 sc_m
% all_sc_m=13;
all_sc_m=30;
all_rot_angle=zeros(length(all_sc_m),1);


for indm=1:length(all_sc_m)
    
    % M > N
    % preparation step
    % M > N
    sc_m=all_sc_m(indm);
    sc_n=sc_m-1;
    
    fft_relax_onoff_factor=1.0;
    graphene_setup_supercell;
    graphene_setup_supercell_2;
    build_supercell_ham;
    
    know=sc_gamma;
    
    graphene_setup_supercell_relax_ver3;
    
    all_rot_angle(indm)=rot_theta/pi*180;
    
    
    graphene_onsite_e=-0.7876;
    
    Hmat_gamma = ham_supercell_N8_Koshino_relax_height_sparse( know/lattice_a, total_num, inplaneGAB_bottom_relax,inplaneGAB_top_relax, inplaneGAABB_bottom_relax,inplaneGAABB_top_relax,inter_relax ,graphene_onsite_e);
    
    % input variables:
    %
    % in Shiang's code, normally I use the length scale lattice_a = 1 !!!
    
    BZknum=10;
    tmpkk=linspace(0,1,BZknum+1);
    tmpkk=tmpkk(1:BZknum);

    [allk1,allk2]=meshgrid(tmpkk,tmpkk);
    allkxy(:,1)=sc_b1(1)*allk1(:)+sc_b2(1)*allk2(:);
    allkxy(:,2)=sc_b1(2)*allk1(:)+sc_b2(2)*allk2(:);
   
    BZ_allham={}
    indcnt=1;
    for indbz=1:(BZknum^2)
        indbz
        know=sc_gamma;
        know(1)=allkxy(indbz,1);
        know(2)=allkxy(indbz,2);

        Hmat_tmp = ham_supercell_N8_Koshino_relax_height_sparse( know/lattice_a, total_num, inplaneGAB_bottom_relax,inplaneGAB_top_relax, inplaneGAABB_bottom_relax,inplaneGAABB_top_relax,inter_relax ,graphene_onsite_e);

        BZ_allham{indcnt}=Hmat_tmp;
        indcnt=indcnt+1;

    end



 
    save(['TwBLG_EffKP_relax_sym_BZscanb_',num2str(sc_m),'_',num2str(sc_n),'_xy'],'-v7.3');
    
    'done'
end

%%

% plot(all_sc_m,all_rot_angle)


