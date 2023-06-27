run(global_var_file);

for indm=1:length(all_sc_m)
    sc_m=all_sc_m(indm);
    sc_n=sc_m-1;

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

    bot_K_point=(-pc_vec_b1+pc_vec_b2)/3
    top_K_point=(-pc_vec_rb1+pc_vec_rb2)/3
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


    qqmat=[sc_b1(1),sc_b2(1);sc_b1(2),sc_b2(2)];
    invqqmat=inv(qqmat);

    for inds=1:size(all_bot_qs,1)
        qqnow=all_bot_qs(inds,1:2)';

        qqconv=invqqmat*qqnow;

        qqconv


        qqnow=all_top_qs(inds,1:2)';

        qqconv2=invqqmat*qqnow;

        qqconv2
    end


% intra layer terms



    bot_intra_list=zeros(2,2,num_hex*length(BZ_allham));
    bot_all_intra_kk=zeros(num_hex*length(BZ_allham),2);

    top_intra_list=zeros(2,2,num_hex*length(BZ_allham));
    top_all_intra_kk=zeros(num_hex*length(BZ_allham),2);

    indcc=1;
    for indk=1:length(BZ_allham)
        Hmatnow=BZ_allham{indk};
        k_scnow=allkxy(indk,1:2);

        for indh=1:num_hex

            kkbot=all_bot_qs(indh,1:2);
            kktop=all_top_qs(indh,1:2);


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

            Eff_intra_bot=wf_set1_b'*Hmatnow*wf_set1_b;
            Eff_intra_top=wf_set1_t'*Hmatnow*wf_set1_t;

            bot_intra_list(:,:,indcc)=Eff_intra_bot(:,:);
            top_intra_list(:,:,indcc)=Eff_intra_top(:,:);


            kkbot_rel=hex_coor(indh,1:2);
            kktop_rel=-hex_coor(indh,1:2);


            bot_all_intra_kk(indcc,:)=kkbot_rel(1:2)+k_scnow(1:2);
            top_all_intra_kk(indcc,:)=kktop_rel(1:2)+k_scnow(1:2);


            indcc=indcc+1;
        end
    end

    %scatter3(bot_all_intra_kk(:,1),bot_all_intra_kk(:,2),abs(bot_intra_list(1,2,:)))

    % data_kx=top_all_intra_kk(:,1);
    % data_ky=top_all_intra_kk(:,2);
    % data_fit=squeeze(top_intra_list(2,1,:));

    data_kx=bot_all_intra_kk(:,1);
    data_ky=bot_all_intra_kk(:,2);
    data_fit=squeeze(bot_intra_list(2,1,:));

    data_kplus=data_kx+1i*data_ky;
    data_kminus=data_kx-1i*data_ky;

    CC_mat=zeros(3);
    CC_mat(1,1)=sum(conj(data_kplus).*data_kplus);
    CC_mat(1,2)=sum(conj(data_kplus).*data_kminus);
    CC_mat(1,3)=sum(conj(data_kplus));

    CC_mat(2,1)=sum(conj(data_kminus).*data_kplus);
    CC_mat(2,2)=sum(conj(data_kminus).*data_kminus);
    CC_mat(2,3)=sum(conj(data_kminus));

    CC_mat(3,1)=sum(data_kplus);
    CC_mat(3,2)=sum(data_kminus);
    CC_mat(3,3)=length(data_fit);

    DDvec(1,1)=sum(conj(data_kplus).*data_fit);
    DDvec(2,1)=sum(conj(data_kminus).*data_fit);
    DDvec(3,1)=sum(data_fit);


    CC_optimal=inv(CC_mat)*DDvec;

    CC_optimal


    % top (1,2):
    %   -0.0448 - 2.1198i
    %    0.0415 - 0.0225i
    %   -0.0010 + 0.0020i

    % top (2,1)
    %    0.0415 + 0.0225i
    %   -0.0448 + 2.1198i
    %   -0.0010 - 0.0020i

    % bot (1,2)
    %    0.0033 - 2.1203i
    %    0.0407 - 0.0239i
    %   -0.0013 + 0.0019i

    % bot (2,1)
    %    0.0407 + 0.0239i
    %    0.0033 + 2.1203i
    %   -0.0013 - 0.0019i


    data_optimal=CC_optimal(1)*data_kplus+CC_optimal(2)*data_kminus+CC_optimal(3);

    %scatter3(data_kx,data_ky,abs(data_fit-data_optimal))


    % intra layer terms


    % in plane terms


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
    
    INTRAall_given_qs


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


        list_now=All_intra_bot_list{indqq};
        sizetmp=size(list_now,1);

        indcc=1;

        collect_list_bot=zeros(2,2,sizetmp*length(BZ_allham));
        all_collect_kk_bot=zeros(sizetmp*length(BZ_allham),2);


        for indk=1:length(BZ_allham)
            Hmatnow=BZ_allham{indk};
            k_scnow=allkxy(indk,1:2);

            for indl=1:sizetmp
                tmpnn=list_now(indl,1:2);

                kkbot1=all_bot_qs(tmpnn(2),1:2);
                kkbot2=all_bot_qs(tmpnn(1),1:2);



                wf_set1_b1=zeros(tot_num,2);

                tmp_k=kkbot1(1:2);
                tmp_k=reshape(tmp_k,2,1);
                tmp_pos=sc_all_points_shift(bot_A_atoms,:);
                wf_set1_b1(bot_A_atoms,1)=exp(i*tmp_pos*tmp_k);
                tmp_pos=sc_all_points_shift(bot_B_atoms,:);
                wf_set1_b1(bot_B_atoms,2)=exp(i*tmp_pos*tmp_k);
                wf_set1_b1=wf_set1_b1/sqrt(tot_num/4);

                wf_set1_b2=zeros(tot_num,2);

                tmp_k=kkbot2(1:2);
                tmp_k=reshape(tmp_k,2,1);
                tmp_pos=sc_all_points_shift(bot_A_atoms,:);
                wf_set1_b2(bot_A_atoms,1)=exp(i*tmp_pos*tmp_k);
                tmp_pos=sc_all_points_shift(bot_B_atoms,:);
                wf_set1_b2(bot_B_atoms,2)=exp(i*tmp_pos*tmp_k);
                wf_set1_b2=wf_set1_b2/sqrt(tot_num/4);

                Eff_intra_bot_tmp=wf_set1_b2'*Hmatnow*wf_set1_b1;


                kkbot_rel=hex_coor(tmpnn(2),1:2);


                all_collect_kk_bot(indcc,:)=kkbot_rel(1:2)+k_scnow(1:2);
                collect_list_bot(:,:,indcc)=Eff_intra_bot_tmp(:,:);

                indcc=indcc+1;
            end
        end


        QQ_collect_list_bot{indqq}=collect_list_bot;
        QQ_all_collect_kk_bot{indqq}=all_collect_kk_bot;


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

    % all_intraQ_bot_CCopt
    % all_intraQ_top_CCopt
    % order from fitting: kplus, kminus, const term

    %scatter3(all_collect_kk(:,1),all_collect_kk(:,2),squeeze(collect_list(ind1,ind2,:)));




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
    %numqs=12;


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


    All_Eff_intra_bot_ext=all_intraQ_bot_CCopt;
    All_Eff_intra_top_ext=all_intraQ_top_CCopt;
    All_Eff_inter_ext=all_interQ_CCopt;
    
    kp_data_dir = '/home/stc/devspace/codes/shiang_kp_relaxed_tblg/kp_runs/0p5_TBH_validate/BZ_sweep/data/kp_data';
    fprintf(['saving file: ' kp_filename_here '.mat \n']);
    save([kp_data_dir '/' kp_filename_here],'All_Eff_intra_bot_ext','All_Eff_intra_top_ext','All_Eff_inter_ext');
end







