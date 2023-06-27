
run(global_var_file);

for indm=1:length(all_sc_m)
   
    % M > N
    % preparation step
    % M > N
    sc_m=all_sc_m(indm);
    sc_n=sc_m-1;

    [sc_m,sc_n]
 
    tbh_filename_here = ['TwBLG_BZscan_',num2str(sc_m),'_',num2str(sc_n),'_' interlayer_model '-inter'];
    kp_filename_here = ['TwBLG_EffKP_',num2str(sc_m),'_',num2str(sc_n),'_' interlayer_model '-inter'];
    
    fprintf(['loading file: ' tbh_filename_here '.mat \n']);
    load([tbh_data_dir '/' tbh_filename_here]);
    
    %load(['TwBLG_EffKP_relax_symtest1_',num2str(sc_m),'_',num2str(sc_n),'_xy']);
    
    %savename=['TwBLG_EffKP_symparms1_',num2str(sc_m),'_',num2str(sc_n),'_xyz'];
    
    % note the sc_bi are in units of 1/lattice_a
    
    hex_cut=5.51;
    hex_cut=hex_cut*sqrt(dot(sc_b1,sc_b1));
    
    hex_shift=(-sc_b1+sc_b2)/3;
    
    hex_M=40;
    
    hex_table=zeros(2*hex_M+1);
    
    hex_index=0;
    hex_coor=0;
    
    ind=1;
    for ind1=(-hex_M):hex_M
        for ind2=(-hex_M):hex_M
            vec=sc_b1*ind1+sc_b2*ind2+hex_shift;
            
            if sqrt(dot(vec,vec))<hex_cut
                %hex_index(ind,1:2)=[ind1,ind2];
                hex_table(ind1+hex_M+1,ind2+hex_M+1)=ind;
                hex_coor(ind,1:2)=vec(1:2);
                
                ind=ind+1;
            end
            
        end
    end
    
    
    
    
    % find all sym qs
    
    bot_K_point=(-pc_vec_b1+pc_vec_b2)/3
    top_K_point=(-pc_vec_rb1+pc_vec_rb2)/3
    hex_shift=(-sc_b1+sc_b2)/3;
    
    all_qs(1,:)=hex_shift;
    all_qs(2,:)=hex_shift+sc_b2;
    all_qs(3,:)=hex_shift-sc_b2;
    all_qs(4,:)=hex_shift-2*sc_b2;
    all_qs(5,:)=hex_shift-sc_b1;
    all_qs(6,:)=hex_shift-sc_b1-sc_b2;
    all_qs(7,:)=hex_shift-sc_b1-2*sc_b2;
    all_qs(8,:)=hex_shift+sc_b1;
    all_qs(9,:)=hex_shift+sc_b1+sc_b2;
    all_qs(10,:)=hex_shift+sc_b1-sc_b2;
    all_qs(11,:)=hex_shift+2*sc_b1;
    all_qs(12,:)=hex_shift+2*sc_b1+sc_b2;
    
    all_bot_qs=all_qs;
    all_bot_qs(:,1)=all_bot_qs(:,1)+bot_K_point(1);
    all_bot_qs(:,2)=all_bot_qs(:,2)+bot_K_point(2);
    all_top_qs=-all_qs;
    all_top_qs(:,1)=all_top_qs(:,1)+top_K_point(1);
    all_top_qs(:,2)=all_top_qs(:,2)+top_K_point(2);
    
    % in plane terms
    
    
    all_given_qs(1,:)=sc_b1;
    all_given_qs(2,:)=sc_b1+sc_b2;
    all_given_qs(3,:)=sc_b2;
    all_given_qs(4,:)=-sc_b1;
    all_given_qs(5,:)=-sc_b1-sc_b2;
    all_given_qs(6,:)=-sc_b2;
    
    all_given_qs(7,:)=2*sc_b1+sc_b2;
    all_given_qs(8,:)=2*sc_b2+sc_b1;
    all_given_qs(9,:)=-sc_b1+sc_b2;
    all_given_qs(10,:)=-2*sc_b1-sc_b2;
    all_given_qs(11,:)=-2*sc_b2-sc_b1;
    all_given_qs(12,:)=sc_b1-sc_b2;
    
    numq_intra=size(all_given_qs);
    numq_intra=numq_intra(1);
    numqs=12;
    
    q_resolution=1E-6;
    
    for indqq=1:numq_intra
        given_q=all_given_qs(indqq,:);
        
        clear Eff_intra_bot;
        clear Eff_intra_top;
        
        indc=0;
        for indq1=1:numqs
            qtmpbot1=all_qs(indq1,:);
            for indq2=1:numqs
                qtmpbot2=all_qs(indq2,:);
                
                qtmpdiff=qtmpbot2-qtmpbot1-given_q;
                if sqrt(dot(qtmpdiff,qtmpdiff))<q_resolution
                    indc=indc+1;
                    
                    bkpt_test1=qtmpbot1+bot_K_point;
                    bkpt_test2=qtmpbot2+bot_K_point;
                    
                    wf_set1_b1=zeros(tot_num,2);
                    
                    tmp_k=bkpt_test1(1:2);
                    tmp_k=reshape(tmp_k,2,1);
                    tmp_pos=sc_all_points_shift(bot_A_atoms,:);
                    wf_set1_b1(bot_A_atoms,1)=exp(i*tmp_pos*tmp_k);
                    tmp_pos=sc_all_points_shift(bot_B_atoms,:);
                    wf_set1_b1(bot_B_atoms,2)=exp(i*tmp_pos*tmp_k);
                    wf_set1_b1=wf_set1_b1/sqrt(tot_num/4);
                    
                    wf_set1_b2=zeros(tot_num,2);
                    
                    tmp_k=bkpt_test2(1:2);
                    tmp_k=reshape(tmp_k,2,1);
                    tmp_pos=sc_all_points_shift(bot_A_atoms,:);
                    wf_set1_b2(bot_A_atoms,1)=exp(i*tmp_pos*tmp_k);
                    tmp_pos=sc_all_points_shift(bot_B_atoms,:);
                    wf_set1_b2(bot_B_atoms,2)=exp(i*tmp_pos*tmp_k);
                    wf_set1_b2=wf_set1_b2/sqrt(tot_num/4);
                    
                    Eff_intra_bot(:,:,indc)=wf_set1_b2'*Hmat_gamma*wf_set1_b1;
                    
                end
                
                
            end
        end
        
        Eff_intra_bot_avg=sum(Eff_intra_bot,3)/indc;
        All_Eff_intra_bot(:,:,indqq)=Eff_intra_bot_avg;
        
        indc=0;
        for indq1=1:numqs
            qtmptop1=-all_qs(indq1,:);
            for indq2=1:numqs
                qtmptop2=-all_qs(indq2,:);
                
                qtmpdiff=qtmptop2-qtmptop1-given_q;
                if sqrt(dot(qtmpdiff,qtmpdiff))<q_resolution
                    indc=indc+1;
                    
                    tkpt_test1=qtmptop1+top_K_point;
                    tkpt_test2=qtmptop2+top_K_point;
                    
                    wf_set1_t1=zeros(tot_num,2);
                    
                    tmp_k=tkpt_test1(1:2);
                    tmp_k=reshape(tmp_k,2,1);
                    tmp_pos=sc_all_points_shift(top_A_atoms,:);
                    wf_set1_t1(top_A_atoms,1)=exp(i*tmp_pos*tmp_k);
                    tmp_pos=sc_all_points_shift(top_B_atoms,:);
                    wf_set1_t1(top_B_atoms,2)=exp(i*tmp_pos*tmp_k);
                    wf_set1_t1=wf_set1_t1/sqrt(tot_num/4);
                    
                    wf_set1_t2=zeros(tot_num,2);
                    
                    tmp_k=tkpt_test2(1:2);
                    tmp_k=reshape(tmp_k,2,1);
                    tmp_pos=sc_all_points_shift(top_A_atoms,:);
                    wf_set1_t2(top_A_atoms,1)=exp(i*tmp_pos*tmp_k);
                    tmp_pos=sc_all_points_shift(top_B_atoms,:);
                    wf_set1_t2(top_B_atoms,2)=exp(i*tmp_pos*tmp_k);
                    wf_set1_t2=wf_set1_t2/sqrt(tot_num/4);
                    
                    Eff_intra_top(:,:,indc)=wf_set1_t2'*Hmat_gamma*wf_set1_t1;
                    
                end
                
                
            end
        end
        
        Eff_intra_top_avg=sum(Eff_intra_top,3)/indc;
        All_Eff_intra_top(:,:,indqq)=Eff_intra_top_avg;
        
        
    end
    
    
    
    % inter plane terms
    
    q_resolution=1E-6;
    
    interall_given_qs(1,:)=(-2*hex_shift+sc_b2);
    interall_given_qs(2,:)=-hex_shift-(hex_shift+sc_b1);
    interall_given_qs(3,:)=(-hex_shift-sc_b1)-(hex_shift-sc_b2);
    
    interall_given_qs(4,:)=-2*hex_shift;
    interall_given_qs(5,:)=-2*hex_shift+2*sc_b2;
    interall_given_qs(6,:)=-2*hex_shift-2*sc_b1;
    
    interall_given_qs(7,:)=(-2*hex_shift+sc_b2)+sc_b1;
    interall_given_qs(8,:)=(-2*hex_shift+sc_b2)+sc_b1+sc_b2;
    
    interall_given_qs(9,:)=(-2*hex_shift+sc_b2)-sc_b1+sc_b2;
    interall_given_qs(10,:)=(-2*hex_shift+sc_b2)-2*sc_b1;
    
    interall_given_qs(11,:)=(-2*hex_shift+sc_b2)-(sc_b1+sc_b2)-sc_b2;
    interall_given_qs(12,:)=(-2*hex_shift+sc_b2)-(sc_b1+sc_b2)*2;
    
    
    
    numq_inter=size(interall_given_qs);
    numq_inter=numq_inter(1);
    numqs=12;
    
    q_resolution=1E-6;
    
    for indqq=1:numq_inter
        given_q=interall_given_qs(indqq,:);
        
        clear Eff_inter;
        
        
        indc=0;
        for indq1=1:numqs
            qtmptop=-all_qs(indq1,:);
            for indq2=1:numqs
                qtmpbot=all_qs(indq2,:);
                
                qtmpdiff=qtmptop-qtmpbot-given_q;
                if sqrt(dot(qtmpdiff,qtmpdiff))<q_resolution
                    indc=indc+1;
                    
                    bkpt_test=qtmpbot+bot_K_point;
                    tkpt_test=qtmptop+top_K_point;
                    
                    wf_set1_b=zeros(tot_num,2);
                    
                    tmp_k=bkpt_test(1:2);
                    tmp_k=reshape(tmp_k,2,1);
                    tmp_pos=sc_all_points_shift(bot_A_atoms,:);
                    wf_set1_b(bot_A_atoms,1)=exp(i*tmp_pos*tmp_k);
                    tmp_pos=sc_all_points_shift(bot_B_atoms,:);
                    wf_set1_b(bot_B_atoms,2)=exp(i*tmp_pos*tmp_k);
                    wf_set1_b=wf_set1_b/sqrt(tot_num/4);
                    
                    wf_set1_t=zeros(tot_num,2);
                    
                    tmp_k=tkpt_test(1:2);
                    tmp_k=reshape(tmp_k,2,1);
                    tmp_pos=sc_all_points_shift(top_A_atoms,:);
                    wf_set1_t(top_A_atoms,1)=exp(i*tmp_pos*tmp_k);
                    tmp_pos=sc_all_points_shift(top_B_atoms,:);
                    wf_set1_t(top_B_atoms,2)=exp(i*tmp_pos*tmp_k);
                    wf_set1_t=wf_set1_t/sqrt(tot_num/4);
                    
                    Eff_inter(:,:,indc)=wf_set1_t'*Hmat_gamma*wf_set1_b;
                    
                end
                
                
            end
        end
        indc
        
        Eff_inter_avg=sum(Eff_inter,3)/indc;
        All_Eff_inter(:,:,indqq)=Eff_inter_avg;
        
        
    end
    
    
    
    
    Intra_qs=all_given_qs/lattice_a;
    Inter_qs=interall_given_qs/lattice_a;
    
    fprintf(['saving file: ' kp_filename_here '.mat \n']);
    save([kp_data_dir '/' kp_filename_here],'Intra_qs','All_Eff_intra_top','All_Eff_intra_bot','Inter_qs','All_Eff_inter');
    
    
end




