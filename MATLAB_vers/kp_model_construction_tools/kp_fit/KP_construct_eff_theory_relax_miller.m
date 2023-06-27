run(global_var_file);

for indm=1:length(all_sc_m)
    sc_m=all_sc_m(indm);
    sc_n=sc_m-1;
    
    intra_shell_max = 10;
    inter_shell_max = 10;

    % here the length scale adopts a=2.46A = 1

    tbh_filename_here = ['TwBLG_BZscan_',num2str(sc_m),'_',num2str(sc_n),'_' interlayer_model '-inter'];
    kp_filename_here = ['TwBLG_EffKP_',num2str(sc_m),'_',num2str(sc_n),'_' interlayer_model '-inter'];

    fprintf(['loading file: ' tbh_filename_here '.mat \n']);
    load([tbh_data_dir '/' tbh_filename_here]);

    
    % note the sc_bi are in units of 1/lattice_a

    hex_cut=3.11;
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

    num_hex=ind-1;


    % find all sym qs

    bot_K_point=(-pc_vec_b1+pc_vec_b2)/3;
    top_K_point=(-pc_vec_rb1+pc_vec_rb2)/3;
    hex_shift=(-sc_b1+sc_b2)/3;
    %
    % all_qs(1,:)=hex_shift;
    % all_qs(2,:)=hex_shift+sc_b2;
    % all_qs(3,:)=hex_shift-sc_b2;
    % all_qs(4,:)=hex_shift-2*sc_b2;
    % all_qs(5,:)=hex_shift-sc_b1;
    % all_qs(6,:)=hex_shift-sc_b1-sc_b2;
    % all_qs(7,:)=hex_shift-sc_b1-2*sc_b2;
    % all_qs(8,:)=hex_shift+sc_b1;
    % all_qs(9,:)=hex_shift+sc_b1+sc_b2;
    % all_qs(10,:)=hex_shift+sc_b1-sc_b2;
    % all_qs(11,:)=hex_shift+2*sc_b1;
    % all_qs(12,:)=hex_shift+2*sc_b1+sc_b2;

    all_bot_qs=hex_coor;
    all_bot_qs(:,1)=all_bot_qs(:,1)+bot_K_point(1);
    all_bot_qs(:,2)=all_bot_qs(:,2)+bot_K_point(2);
    all_top_qs=-hex_coor;
    all_top_qs(:,1)=all_top_qs(:,1)+top_K_point(1);
    all_top_qs(:,2)=all_top_qs(:,2)+top_K_point(2);

    % INTRA CONSTRUCTION
    %%{
    
    % intra layer terms

    % fast hexagonal miller index construction
    
    grid_max = 15;
    shell_max = intra_shell_max;

    shell_max = shell_max + 1; %origin counts as a shell, so +1
    grid_base = [-grid_max:grid_max];

    [a1_idx, a2_idx] = meshgrid(grid_base,grid_base);

    idx = 1;
    for a1_r = 0:grid_max
        a2_min = ceil(a1_r/2);
        for a2_r = -(a2_min:a1_r)
            a3_r = -(a1_r + a2_r);
            r_shells(idx) = a1_r^2 + a2_r^2 + a3_r^2;
            idx = idx + 1;
        end
    end

    % uncomment to get all shells found within grid search
    %shell_max = length(r_shells);

    r_counts = zeros(shell_max,1);

    for r_idx = 1:shell_max
       idx_shells{r_idx} = [];
    end

    a3_idx = -(a1_idx + a2_idx);
    idx = 1;
    
    intra_idx = 1;
    shell_to_INTRAall = -ones(shell_max,24);
    miller_b1 = sc_b1;
    miller_b2 = sc_b1 + sc_b2;
    
    for r_idx = 2:shell_max % skip constant
        for x = 1:size(a1_idx,1)
            for y = 1:size(a2_idx,2)

            r_val = a1_idx(x,y)^2 + a2_idx(x,y)^2 + a3_idx(x,y)^2;

                if r_val == r_shells(r_idx) % if the radius here matches this shell's radius

                    idx_here = idx_shells{r_idx};
                    idx_here(r_counts(r_idx)+1,:) = [a1_idx(x,y) a2_idx(x,y) a3_idx(x,y)];
                    idx_shells{r_idx} = idx_here;
                    %fprintf('[%d, %d, %d], shell %d \n',idx_here(r_counts(r_idx)+1,1),idx_here(r_counts(r_idx)+1,2),idx_here(r_counts(r_idx)+1,3),r_idx);

                    INTRAall_given_qs(intra_idx,:)= miller_b1*idx_here(r_counts(r_idx)+1,1) + miller_b2*idx_here(r_counts(r_idx)+1,2);
                    INTRAall_to_shell(intra_idx,:) = [r_idx r_counts(r_idx)+1]; % save both indices
                    shell_to_INTRAall(r_idx, r_counts(r_idx)+1) = intra_idx;

                    r_counts(r_idx) = r_counts(r_idx) + 1; % keeps track of index within shell
                    intra_idx = intra_idx + 1; % keeps track of total number of intra-layer momenta

                end

            end

        end
    end
    

    % old hardcoded intra momenta
    %{
    INTRAall_given_qs(1,:)=sc_b1;
    INTRAall_given_qs(2,:)=sc_b1+sc_b2;
    INTRAall_given_qs(3,:)=sc_b2;
    INTRAall_given_qs(4,:)=-sc_b1;
    INTRAall_given_qs(5,:)=-sc_b1-sc_b2;
    INTRAall_given_qs(6,:)=-sc_b2;

    INTRAall_given_qs(7,:)=2*sc_b1+sc_b2;
    INTRAall_given_qs(8,:)=2*sc_b2+sc_b1;
    INTRAall_given_qs(9,:)=-sc_b1+sc_b2;
    INTRAall_given_qs(10,:)=-2*sc_b1-sc_b2;
    INTRAall_given_qs(11,:)=-2*sc_b2-sc_b1;
    INTRAall_given_qs(12,:)=sc_b1-sc_b2;
    %}

    %scatter(INTRAall_given_qs(:,1),INTRAall_given_qs(:,2))

    numq_intra=size(INTRAall_given_qs);
    numq_intra=numq_intra(1);
    %numqs_inter=12;


    All_intra_bot_list={};
    All_intra_top_list={};

    q_resolution=1E-6;

    for indqq=1:numq_intra
        given_q=INTRAall_given_qs(indqq,1:2);

        %clear Eff_inter;


        tmp_list=[];
        indc=0;
        for indq1=1:num_hex
            qtmptop=hex_coor(indq1,:);
            for indq2=1:num_hex
                qtmpbot=hex_coor(indq2,:);

                qtmpdiff=qtmptop-qtmpbot-given_q;
                if sqrt(dot(qtmpdiff,qtmpdiff))<q_resolution
                    indc=indc+1;

                    tmp_list(indc,1:2)=[indq1,indq2]; %,qtmptop(1)-qtmpbot(1),qtmptop(2)-qtmpbot(2)];


                end


            end
        end
        indc
        All_intra_bot_list{indqq}=tmp_list;

        tmp_list=[];
        indc=0;
        for indq1=1:num_hex
            qtmptop=-hex_coor(indq1,:);
            for indq2=1:num_hex
                qtmpbot=-hex_coor(indq2,:);

                qtmpdiff=qtmptop-qtmpbot-given_q;
                if sqrt(dot(qtmpdiff,qtmpdiff))<q_resolution
                    indc=indc+1;

                    tmp_list(indc,1:2)=[indq1,indq2]; %,qtmptop(1)-qtmpbot(1),qtmptop(2)-qtmpbot(2)];


                end


            end
        end
        indc
        All_intra_top_list{indqq}=tmp_list;

    end
    
    %


    QQ_collect_list_bot={};
    QQ_all_collect_kk_bot={};
    QQ_collect_list_top={};
    QQ_all_collect_kk_top={};

    for indqq=1:numq_intra
        indqq

        %given_q=INTRAall_given_qs(indqq,1:2);


        list_now=All_intra_bot_list{indqq}; % index pairs that have correct dk
        sizetmp=size(list_now,1);

        indcc=1;

        collect_list_bot=zeros(2,2,sizetmp*length(BZ_allham));
        all_collect_kk_bot=zeros(sizetmp*length(BZ_allham),2);


        for indk=1:length(BZ_allham)
            Hmatnow=BZ_allham{indk};
            k_scnow=allkxy(indk,1:2);

            for indl=1:sizetmp
                tmpnn=list_now(indl,1:2); % pair of hex_coor that have correct dk

                kkbot1=all_bot_qs(tmpnn(2),1:2); % momenta of first in pair
                kkbot2=all_bot_qs(tmpnn(1),1:2); % momenta of second in pair


                % construct the wavefunction for first momenta
                wf_set1_b1=zeros(tot_num,2);

                tmp_k=kkbot1(1:2);
                tmp_k=reshape(tmp_k,2,1);
                tmp_pos=sc_all_points_shift(bot_A_atoms,:);
                wf_set1_b1(bot_A_atoms,1)=exp(i*tmp_pos*tmp_k);
                tmp_pos=sc_all_points_shift(bot_B_atoms,:);
                wf_set1_b1(bot_B_atoms,2)=exp(i*tmp_pos*tmp_k);
                wf_set1_b1=wf_set1_b1/sqrt(tot_num/4);

                % construct the wavefunction for second momenta
                wf_set1_b2=zeros(tot_num,2);

                tmp_k=kkbot2(1:2);
                tmp_k=reshape(tmp_k,2,1);
                tmp_pos=sc_all_points_shift(bot_A_atoms,:);
                wf_set1_b2(bot_A_atoms,1)=exp(i*tmp_pos*tmp_k);
                tmp_pos=sc_all_points_shift(bot_B_atoms,:);
                wf_set1_b2(bot_B_atoms,2)=exp(i*tmp_pos*tmp_k);
                wf_set1_b2=wf_set1_b2/sqrt(tot_num/4);

                % overlap integral
                Eff_intra_bot_tmp=wf_set1_b2'*Hmatnow*wf_set1_b1;


                % relative momenta (hex_coor is before layer momentum shift)
                kkbot_rel=hex_coor(tmpnn(2),1:2);


                all_collect_kk_bot(indcc,:)=kkbot_rel(1:2)+k_scnow(1:2);
                collect_list_bot(:,:,indcc)=Eff_intra_bot_tmp(:,:);

                indcc=indcc+1;
            end
        end


        QQ_collect_list_bot{indqq}=collect_list_bot;
        QQ_all_collect_kk_bot{indqq}=all_collect_kk_bot;


        % same thing as above, but for the top layer!
        list_now=All_intra_top_list{indqq};
        sizetmp=size(list_now,1);

        indcc=1;

        collect_list_top=zeros(2,2,sizetmp*length(BZ_allham));
        all_collect_kk_top=zeros(sizetmp*length(BZ_allham),2);


        for indk=1:length(BZ_allham)
            Hmatnow=BZ_allham{indk};
            k_scnow=allkxy(indk,1:2);

            for indl=1:sizetmp
                tmpnn=list_now(indl,1:2);

                kktop1=all_top_qs(tmpnn(2),1:2);
                kktop2=all_top_qs(tmpnn(1),1:2);



                wf_set1_t1=zeros(tot_num,2);

                tmp_k=kktop1(1:2);
                tmp_k=reshape(tmp_k,2,1);
                tmp_pos=sc_all_points_shift(top_A_atoms,:);
                wf_set1_t1(top_A_atoms,1)=exp(i*tmp_pos*tmp_k);
                tmp_pos=sc_all_points_shift(top_B_atoms,:);
                wf_set1_t1(top_B_atoms,2)=exp(i*tmp_pos*tmp_k);
                wf_set1_t1=wf_set1_t1/sqrt(tot_num/4);

                wf_set1_t2=zeros(tot_num,2);

                tmp_k=kktop2(1:2);
                tmp_k=reshape(tmp_k,2,1);
                tmp_pos=sc_all_points_shift(top_A_atoms,:);
                wf_set1_t2(top_A_atoms,1)=exp(i*tmp_pos*tmp_k);
                tmp_pos=sc_all_points_shift(top_B_atoms,:);
                wf_set1_t2(top_B_atoms,2)=exp(i*tmp_pos*tmp_k);
                wf_set1_t2=wf_set1_t2/sqrt(tot_num/4);

                Eff_intra_top_tmp=wf_set1_t2'*Hmatnow*wf_set1_t1;


                kktop_rel=-hex_coor(tmpnn(2),1:2);


                all_collect_kk_top(indcc,:)=kktop_rel(1:2)+k_scnow(1:2);
                collect_list_top(:,:,indcc)=Eff_intra_top_tmp(:,:);

                indcc=indcc+1;
            end
        end

        QQ_collect_list_top{indqq}=collect_list_top;
        QQ_all_collect_kk_top{indqq}=all_collect_kk_top;



    end


    % intra terms analysis:


    all_intraQ_bot_CCopt=zeros(3,2,2,numq_intra);
    all_intraQ_top_CCopt=zeros(3,2,2,numq_intra);

    for indqq=1:numq_intra
        indqq

        collect_list=QQ_collect_list_bot{indqq};
        all_collect_kk=QQ_all_collect_kk_bot{indqq};

        All_CCopt=zeros(3,2,2);

        for ind1=1:2
            for ind2=1:2


                collect_kk=sqrt(all_collect_kk(:,1).^2+all_collect_kk(:,2).^2);
                sel_collect_kk=(collect_kk<=0.3); % only keep momenta with small total (first BZ)

                % k_\pm = kx \pm iky
                kk_plus=all_collect_kk(sel_collect_kk,1)+i*all_collect_kk(sel_collect_kk,2);
                kk_minus=all_collect_kk(sel_collect_kk,1)-i*all_collect_kk(sel_collect_kk,2);


                sel_data=squeeze(collect_list(ind1,ind2,sel_collect_kk));


                CC_mat=zeros(3);
                CC_mat(1,1)=sum(conj(kk_plus).*kk_plus);
                CC_mat(1,2)=sum(conj(kk_plus).*kk_minus);
                CC_mat(1,3)=sum(conj(kk_plus));

                CC_mat(2,1)=sum(conj(kk_minus).*kk_plus);
                CC_mat(2,2)=sum(conj(kk_minus).*kk_minus);
                CC_mat(2,3)=sum(conj(kk_minus));

                CC_mat(3,1)=sum(kk_plus);
                CC_mat(3,2)=sum(kk_minus);
                CC_mat(3,3)=length(kk_minus);

                DDvec=zeros(3,1);
                DDvec(1)=sum(conj(kk_plus).*sel_data);
                DDvec(2)=sum(conj(kk_minus).*sel_data);
                DDvec(3)=sum(sel_data);

                CCoptimal=inv(CC_mat)*DDvec;
                All_CCopt(:,ind1,ind2)=CCoptimal(:);

            end


        end

        all_intraQ_bot_CCopt(:,:,:,indqq)=All_CCopt(:,:,:);




        collect_list=QQ_collect_list_top{indqq};
        all_collect_kk=QQ_all_collect_kk_top{indqq};

        All_CCopt=zeros(3,2,2);

        for ind1=1:2
            for ind2=1:2


                collect_kk=sqrt(all_collect_kk(:,1).^2+all_collect_kk(:,2).^2);
                sel_collect_kk=(collect_kk<=0.3);

                % k_\pm = kx \pm iky
                kk_plus=all_collect_kk(sel_collect_kk,1)+i*all_collect_kk(sel_collect_kk,2);
                kk_minus=all_collect_kk(sel_collect_kk,1)-i*all_collect_kk(sel_collect_kk,2);


                sel_data=squeeze(collect_list(ind1,ind2,sel_collect_kk));


                CC_mat=zeros(3);
                CC_mat(1,1)=sum(conj(kk_plus).*kk_plus);
                CC_mat(1,2)=sum(conj(kk_plus).*kk_minus);
                CC_mat(1,3)=sum(conj(kk_plus));

                CC_mat(2,1)=sum(conj(kk_minus).*kk_plus);
                CC_mat(2,2)=sum(conj(kk_minus).*kk_minus);
                CC_mat(2,3)=sum(conj(kk_minus));

                CC_mat(3,1)=sum(kk_plus);
                CC_mat(3,2)=sum(kk_minus);
                CC_mat(3,3)=length(kk_minus);

                DDvec=zeros(3,1);
                DDvec(1)=sum(conj(kk_plus).*sel_data);
                DDvec(2)=sum(conj(kk_minus).*sel_data);
                DDvec(3)=sum(sel_data);

                CCoptimal=inv(CC_mat)*DDvec;
                All_CCopt(:,ind1,ind2)=CCoptimal(:);

            end


        end

        all_intraQ_top_CCopt(:,:,:,indqq)=All_CCopt(:,:,:);





    end
   
    % enforce symmetry on intra-layer terms
    
    % INTRAall_given_qs(intra index,:) = [k1, k2]
    % INTRAall_to_shell(intra index,:) = [shell index, internal index]
    % shell_to_INTRAall(shell index, internal index) = intra index;
      
    
    % basic variables
    sigx = [0,1;1,0];
    sigy = [0,-1j;1j,0];
    sigz = [1,0;0,-1];
    rot_mat=diag([exp(-1j*2*pi/3),exp(1j*2*pi/3)]);

    
    bot_pre_symm = all_intraQ_bot_CCopt;
    top_pre_symm = all_intraQ_top_CCopt;
    
    bot_post_symm = all_intraQ_bot_CCopt;
    top_post_symm = all_intraQ_top_CCopt;   
    
    % enforce symmetry
    for symm_cond = 1:3 % go over all symmetry conditions (onsite, mirror, rotational)

        for r_idx = 1:shell_max % each shell has its own conditions

            idx_arr = idx_shells{r_idx};

            % evaluate plane-wave with Hamiltonian
            % t_arr = ... [ done above ]

            for c = 1:size(idx_arr,1)
                idx_here = idx_arr(c,:);

                mirror_tar = -idx_here([2,1,3]); % mirror = swap front two values
                rot1_tar = idx_here([3,1,2]); % rot = all even permutations
                rot2_tar = idx_here([2,3,1]);

                mirror_zero = sum((idx_arr - mirror_tar).^2,2);
                rot1_zero = sum((idx_arr - rot1_tar).^2,2);
                rot2_zero = sum((idx_arr - rot2_tar).^2,2);

                [~,mirror_c] = min(mirror_zero);
                [~,rot1_c] = min(rot1_zero);
                [~,rot2_c] = min(rot2_zero);

                % get the INTRA_list index
                intra_h = shell_to_INTRAall(r_idx, c);
                intra_m = shell_to_INTRAall(r_idx, mirror_c);
                intra_r1 = shell_to_INTRAall(r_idx, rot1_c);
                intra_r2 = shell_to_INTRAall(r_idx, rot2_c);

                M_hb  = squeeze(bot_pre_symm(3,:,:,intra_h));
                M_mb  = squeeze(bot_pre_symm(3,:,:,intra_m));
                M_r1b = squeeze(bot_pre_symm(3,:,:,intra_r1));
                M_r2b = squeeze(bot_pre_symm(3,:,:,intra_r2));

                M_ht  = squeeze(top_pre_symm(3,:,:,intra_h));
                M_mt  = squeeze(top_pre_symm(3,:,:,intra_m));
                M_r1t = squeeze(top_pre_symm(3,:,:,intra_r1));
                M_r2t = squeeze(top_pre_symm(3,:,:,intra_r2));                    

                % use symmetry conditions to calculate average term
                % onsite symm.
                if (symm_cond == 1)
                    post_symm_b = (M_hb+sigx*conj(M_hb)*(sigx))/2;
                    post_symm_t = (M_ht+sigx*conj(M_ht)*(sigx))/2;

                end
                % mirror symm.
                if (symm_cond == 2)
                    post_symm_b =(M_hb+sigx*(M_mt)*sigx)/2; % C2 rotation
                    post_symm_t =(M_ht+sigx*(M_mb)*sigx)/2; % so invert layer index!
                end

                % rotational symm.
                if (symm_cond == 3)
                    post_symm_b = (M_hb+(rot_mat')*M_r1b*rot_mat+(rot_mat)*M_r2b*(rot_mat'))/3;
                    post_symm_t = (M_ht+(rot_mat')*M_r1t*rot_mat+(rot_mat)*M_r2t*(rot_mat'))/3;
                end

                    bot_post_symm(3,:,:,intra_h) = post_symm_b;
                    top_post_symm(3,:,:,intra_h) = post_symm_t;

            end

        end
        
        bot_pre_symm = bot_post_symm;
        top_pre_symm = top_post_symm;
        
    end
    
    All_Eff_intra_shell_indices = idx_shells; % save for later

    
    

    % all_intraQ_bot_CCopt
    % all_intraQ_top_CCopt
    % order from fitting: kplus, kminus, const term

    %scatter3(all_collect_kk(:,1),all_collect_kk(:,2),squeeze(collect_list(ind1,ind2,:)));



    %%}
    % END INTRA CONSTRUCTION
    
    % inter plane terms


    q_resolution=1E-6;

    
    %hex_shift=(-sc_b1+sc_b2)/3;

    % OLD HARD CODED INTERLAYER MOMENTA
    %{
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
    %}
    %numqs=12;

    % interq miller construct
    a_rot = 120;
    
    r_mat = [cosd(a_rot) sind(a_rot); -sind(a_rot) cosd(a_rot)];
    
    a1 = (-2*hex_shift(1:2)+sc_b2(1:2))'; % first INTERQ
    a2 = r_mat*a1;
    a3 = r_mat*a2;

    grid_max = 10;
    shell_max = inter_shell_max;

    %shell_max = shell_max + 1; %origin counts as a shell, so +1
    grid_base = [-grid_max:grid_max];

    [a1_idx, a2_idx] = meshgrid(grid_base,grid_base);

    idx = 1;
    for a1_r = 1:grid_max
        a2_min = 0;%ceil(a1_r/2);
        a2_max = grid_max;
        for a2_r = -(a2_min:a2_max)
            a3_r = -(a1_r + a2_r) + 1;
            r_shells_list(idx) = a1_r^2 + a2_r^2 + a3_r^2;
            idx = idx + 1;
        end
    end

    r_shells = unique(sort(r_shells_list)); % just find the unique shells, in order

    % uncomment to get all shells found within grid search
    %shell_max = length(r_shells);

    r_counts = zeros(shell_max,1);

    clear idx_shells;
    for r_idx = 1:shell_max
       %pos_shells{r_idx} = []; 
       idx_shells{r_idx} = [];
    end

    a3_idx = -(a1_idx + a2_idx)+1;

    %inter_idx = 1;
    shell_to_interall = -ones(shell_max,24);

    
    for x = 1:size(a1_idx,1)
        for y = 1:size(a2_idx,2)

            pos_h = a1*a1_idx(x,y) + a2*a2_idx(x,y) + a3*a3_idx(x,y);
            pos_x(x,y) = pos_h(1);
            pos_y(x,y) = pos_h(2);
            r_vals(x,y) = a1_idx(x,y)^2 + a2_idx(x,y)^2 + a3_idx(x,y)^2;

            for r_idx = 1:shell_max
                if r_vals(x,y) == r_shells(r_idx)

                    %pos_here = pos_shells{r_idx};
                    %pos_here(r_counts(r_idx)+1,:) = [pos_x(x,y) pos_y(x,y)];
                    %pos_shells{r_idx} = pos_here;

                    
                    idx_here = idx_shells{r_idx};
                    idx_here(r_counts(r_idx)+1,:) = [a1_idx(x,y) a2_idx(x,y) a3_idx(x,y)];
                    idx_shells{r_idx} = idx_here;
                    %fprintf('[%d, %d, %d], shell %d \n',idx_here(r_counts(r_idx)+1,1),idx_here(r_counts(r_idx)+1,2),idx_here(r_counts(r_idx)+1,3),r_idx);

                    %interall_given_qs(inter_idx,:)= pos_h(1:2);
                    %interall_to_shell(inter_idx,:) = [r_idx r_counts(r_idx)+1]; % save both shell number and shell index
                    %shell_to_interall(r_idx, r_counts(r_idx)+1) = inter_idx;

                    r_counts(r_idx) = r_counts(r_idx) + 1; % keeps track of index within shell
                    %inter_idx = inter_idx + 1; % keeps track of total number of inter-layer momenta
                    

                end

                %index_list(x,y) = idx;
                %pos_list(idx,:) = pos_h;
                %r_list(idx) = r_vals(x,y);
                %idx = idx+1;
            end


        end
    end
    
    inter_idx = 1;
    for r = 1:shell_max
        shell_here = idx_shells{r};
        for idx = 1:size(shell_here,1) % go through every momenta in this shell
    
            indices_here = shell_here(idx,:);
            
            interall_given_qs(inter_idx,:)= a1*indices_here(1) + a2*indices_here(2) + a3*indices_here(3);
            interall_to_shell(inter_idx,:) = [r, idx]; % save both shell number and shell index
            shell_to_interall(r, idx) = inter_idx;

            inter_idx = inter_idx + 1; % keeps track of total number of inter-layer momenta
        end
    end

    numq_inter=size(interall_given_qs,1);
    
    All_Eff_inter_shell_indices = idx_shells; % save for later

        
% end of Miller interq construction
    
    
    All_inter_list={};

    for indqq=1:numq_inter
        given_q=interall_given_qs(indqq,1:2);

        clear Eff_inter;


        tmp_list=[];
        indc=0;
        for indq1=1:num_hex
            qtmptop=-hex_coor(indq1,:);
            for indq2=1:num_hex
                qtmpbot=hex_coor(indq2,:);

                qtmpdiff=qtmptop-qtmpbot-given_q;
                if sqrt(dot(qtmpdiff,qtmpdiff))<q_resolution
                    indc=indc+1;

                    tmp_list(indc,1:2)=[indq1,indq2]; %,qtmptop(1)-qtmpbot(1),qtmptop(2)-qtmpbot(2)];


                end


            end
        end
        %indc
        All_inter_list{indqq}=tmp_list;

    end



    %{
    figure(3)
    scatter(interall_given_qs(:,1),interall_given_qs(:,2))
    hold on;
    scatter([0],[0],'r')
    hold off;
    axis equal;
    %}

    %


    all_interQ_CCopt=zeros(3,2,2,numq_inter);

    for indqq=1:numq_inter

        indqq

        list_now=All_inter_list{indqq};
        sizetmp=size(list_now,1);

        indcc=1;

        collect_list=zeros(2,2,sizetmp*length(BZ_allham));
        all_collect_kk=zeros(sizetmp*length(BZ_allham),2);

        for indk=1:length(BZ_allham)
            Hmatnow=BZ_allham{indk};
            k_scnow=allkxy(indk,1:2);

            for indl=1:sizetmp
                tmpnn=list_now(indl,1:2);

                kkbot=all_bot_qs(tmpnn(2),1:2);
                kktop=all_top_qs(tmpnn(1),1:2);


                bkpt_test=kkbot;
                tkpt_test=kktop;

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

                Eff_inter_tmp=wf_set1_t'*Hmatnow*wf_set1_b;


                kkbot_rel=hex_coor(tmpnn(2),1:2);


                all_collect_kk(indcc,:)=kkbot_rel(1:2)+k_scnow(1:2);
                collect_list(:,:,indcc)=Eff_inter_tmp(:,:);

                indcc=indcc+1;
            end






        end

        All_CCopt=zeros(3,2,2);

        for ind1=1:2
            for ind2=1:2


                collect_kk=sqrt(all_collect_kk(:,1).^2+all_collect_kk(:,2).^2);
                sel_collect_kk=(collect_kk<=0.1);

                % k_\pm = kx \pm iky
                kk_plus=all_collect_kk(sel_collect_kk,1)+i*all_collect_kk(sel_collect_kk,2);
                kk_minus=all_collect_kk(sel_collect_kk,1)-i*all_collect_kk(sel_collect_kk,2);


                sel_data=squeeze(collect_list(ind1,ind2,sel_collect_kk));


                CC_mat=zeros(3);
                CC_mat(1,1)=sum(conj(kk_plus).*kk_plus);
                CC_mat(1,2)=sum(conj(kk_plus).*kk_minus);
                CC_mat(1,3)=sum(conj(kk_plus));

                CC_mat(2,1)=sum(conj(kk_minus).*kk_plus);
                CC_mat(2,2)=sum(conj(kk_minus).*kk_minus);
                CC_mat(2,3)=sum(conj(kk_minus));

                CC_mat(3,1)=sum(kk_plus);
                CC_mat(3,2)=sum(kk_minus);
                CC_mat(3,3)=length(kk_minus);

                DDvec=zeros(3,1);
                DDvec(1)=sum(conj(kk_plus).*sel_data);
                DDvec(2)=sum(conj(kk_minus).*sel_data);
                DDvec(3)=sum(sel_data);
                
                CCoptimal=inv(CC_mat)*DDvec;
                All_CCopt(:,ind1,ind2)=CCoptimal(:);

            end


        end

        all_interQ_CCopt(:,:,:,indqq)=All_CCopt(:,:,:);

    end
    
    % enforce symmetry on inter-layer terms
    
    % interall_given_qs(intra index,:) = [k1, k2]nt
    % interall_to_shell(intra index,:) = [shell index, internal index]
    % shell_to_interall(shell index, internal index) = intra index;
      
    
    % basic variables
    sigx = [0,1;1,0];
    sigy = [0,-1j;1j,0];
    sigz = [1,0;0,-1];
    rot_mat=diag([exp(-1j*2*pi/3),exp(1j*2*pi/3)]);

    
    inter_pre_symm = all_interQ_CCopt;
    inter_post_symm = all_interQ_CCopt;
    
    % enforce symmetry
    for symm_cond = 1:3 % go over symmetry conditions (NO ONSITE, but has mirror, rotational)

        for r_idx = 1:shell_max % each shell has its own conditions

            idx_arr = idx_shells{r_idx};

            % evaluate plane-wave with Hamiltonian
            % t_arr = ... [ done above ]

            for c = 1:size(idx_arr,1)
                idx_here = idx_arr(c,:);

                mirror_tar = idx_here([1,3,2]); % mirror = swap last two values
                rot1_tar = idx_here([3,1,2]); % rot = all even permutations
                rot2_tar = idx_here([2,3,1]);

                mirror_zero = sum((idx_arr - mirror_tar).^2,2);
                rot1_zero = sum((idx_arr - rot1_tar).^2,2);
                rot2_zero = sum((idx_arr - rot2_tar).^2,2);

                [~,mirror_c] = min(mirror_zero);
                [~,rot1_c] = min(rot1_zero);
                [~,rot2_c] = min(rot2_zero);

                % get the interall index
                inter_h = shell_to_interall(r_idx, c);
                inter_m = shell_to_interall(r_idx, mirror_c);
                inter_r1 = shell_to_interall(r_idx, rot1_c);
                inter_r2 = shell_to_interall(r_idx, rot2_c);

                M_h  = squeeze(inter_pre_symm(3,:,:,inter_h));
                M_m  = squeeze(inter_pre_symm(3,:,:,inter_m));
                M_r1 = squeeze(inter_pre_symm(3,:,:,inter_r1));
                M_r2 = squeeze(inter_pre_symm(3,:,:,inter_r2));                    

                % use symmetry conditions to calculate average term
                % onsite symm.
                if (symm_cond == 1)
                    post_symm_temp = (M_h+sigx*conj(M_h)*(sigx))/2;

                end
                % mirror symm.
                if (symm_cond == 2)
                    post_symm_temp =(M_h+sigx*(M_m')*sigx)/2; % C2 rotation
                end

                % rotational symm.
                if (symm_cond == 3)
                    post_symm_temp = (M_h+(rot_mat')*M_r2*rot_mat+(rot_mat)*M_r1*(rot_mat'))/3;
                end

                inter_post_symm(3,:,:,inter_h) = post_symm_temp;

            end

        end
        
        inter_pre_symm = inter_post_symm;
    end
    

    %interall_given_qs(:,:)

    %interall_given_qs(indqq,:)
    All_CCopt_const=squeeze(All_CCopt(3,:,:));
    %All_CCopt_const
    % angle(All_CCopt_const)/pi*180

    % q=1
    %
    %    0.1116 + 0.0000i   0.1116 + 0.0000i
    %    0.1116 + 0.0000i   0.1116 - 0.0000i

    % q=2
    %    0.1117 + 0.0000i  -0.0559 - 0.0967i
    %   -0.0559 + 0.0967i   0.1117 - 0.0000i

    % q=3
    %    0.1115 + 0.0000i  -0.0557 + 0.0966i
    %   -0.0557 - 0.0966i   0.1115 - 0.0000i


    All_CCopt_kplus=squeeze(All_CCopt(1,:,:));
    %All_CCopt_kplus
    % q=1
    %   -0.0002 + 0.0451i  -0.0002 + 0.0451i
    %   -0.0002 + 0.0451i  -0.0002 + 0.0451i

    % q=2
    %   -0.0395 - 0.0216i   0.0010 + 0.0450i
    %    0.0385 - 0.0234i  -0.0395 - 0.0216i

    % q=3
    %    0.0413 - 0.0244i   0.0005 + 0.0480i
    %   -0.0418 - 0.0236i   0.0413 - 0.0244i



    All_CCopt_kminus=squeeze(All_CCopt(2,:,:));
    %All_CCopt_kminus
    % q=1
    %   -0.0002 - 0.0451i  -0.0002 - 0.0451i
    %   -0.0002 - 0.0451i  -0.0002 - 0.0451i

    % q=2
    %   -0.0395 + 0.0216i   0.0385 + 0.0234i
    %    0.0010 - 0.0450i  -0.0395 + 0.0216i

    % q=3
    %    0.0413 + 0.0244i  -0.0418 + 0.0236i
    %    0.0005 - 0.0480i   0.0413 + 0.0244i


    www=exp(i*2*pi/3);




    % save data for the next phase:


    %All_Eff_intra_bot_ext=all_intraQ_bot_CCopt;
    %All_Eff_intra_top_ext=all_intraQ_top_CCopt;
    %All_Eff_inter_ext=all_interQ_CCopt;
    
    All_Eff_intra_bot_ext = bot_post_symm;
    All_Eff_intra_top_ext = top_post_symm;
    All_Eff_inter_ext = inter_post_symm;

    fprintf(['saving file: ' kp_filename_here '.mat \n']);
    kp_data_dir = '/home/stc/devspace/codes/shiang_kp_relaxed_tblg/kp_runs/0p5_TBH_validate/BZ_sweep/data/kp_data';
    save([kp_data_dir '/' kp_filename_here],'All_Eff_intra_bot_ext','All_Eff_intra_top_ext','All_Eff_inter_ext','All_Eff_intra_shell_indices','All_Eff_inter_shell_indices');
end







