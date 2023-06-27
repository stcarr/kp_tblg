%run(global_var_file);

for indm=1:1%length(all_sc_m)
    
    %{
    sc_m = all_sc_m(1);
    sc_n = sc_m - 1;
    rot_theta=acos((sc_m^2+sc_n^2+4*sc_m*sc_n)/2/(sc_m^2+sc_m*sc_n+sc_n^2));%*180/pi;

    
    kp_filename_here = ['TwBLG_EffKP_',num2str(sc_m),'_',num2str(sc_n),'_' interlayer_model '-inter'];
    bands_filename_here = ['TwBLG_kp-bands_',num2str(sc_m),'_',num2str(sc_n),'_' interlayer_model '-inter'];

    fprintf(['loading file: ' kp_filename_here '.mat \n']);
    load([kp_data_dir '/' kp_filename_here]);
    %}
    All_Eff_inter=squeeze(All_Eff_inter_ext(3,:,:,:));
    All_Eff_intra_bot=squeeze(All_Eff_intra_bot_ext(3,:,:,:));
    All_Eff_intra_top=squeeze(All_Eff_intra_top_ext(3,:,:,:));


    All_Eff_inter_kplus=squeeze(All_Eff_inter_ext(1,:,:,:));
    All_Eff_intra_bot_kplus=squeeze(All_Eff_intra_bot_ext(1,:,:,:));
    All_Eff_intra_top_kplus=squeeze(All_Eff_intra_top_ext(1,:,:,:));

    All_Eff_inter_kminus=squeeze(All_Eff_inter_ext(2,:,:,:));
    All_Eff_intra_bot_kminus=squeeze(All_Eff_intra_bot_ext(2,:,:,:));
    All_Eff_intra_top_kminus=squeeze(All_Eff_intra_top_ext(2,:,:,:));
    


    %[All_Eff_inter,All_Eff_intra_bot,All_Eff_intra_top] = Fit_Model_parms(fname_parm,rot_theta);

    lattice_a=1.42*sqrt(3);

    unit_dim=2;
    % cutoff in momentum space, which would affect accuracies for the solution


    Dirac_v=5.2275;
    % Dirac_v=5.2;

    KD=4*pi/3/lattice_a;
    KTH=2*KD*sin(rot_theta/2);

    HEX_BLEN=KTH*sqrt(3);

    % basic variables
    sigx=[0,1;1,0];
    sigy=[0,-i;i,0];
    sigz=[1,0;0,-1];

    umat=expm(-i*rot_theta*sigz/2);
    rsigx=umat'*sigx*umat;
    rsigy=umat'*sigy*umat;
    rsigz=umat'*sigz*umat;

    hex_b1=HEX_BLEN*[sqrt(3)/2,-1/2];
    hex_b2=HEX_BLEN*[0,1];

    hex_shift=(-hex_b1+hex_b2)/3;



    % create hex table
    num_intra_qs=12;
    intra_qs(1,:)=[1,0];
    intra_qs(2,:)=[1,1];
    intra_qs(3,:)=[0,1];
    intra_qs(4,:)=[-1,0];
    intra_qs(5,:)=[-1,-1];
    intra_qs(6,:)=[0,-1];
    intra_qs(7,:)=[2,1];
    intra_qs(8,:)=[1,2];
    intra_qs(9,:)=[-1,1];
    intra_qs(10,:)=[-2,-1];
    intra_qs(11,:)=[-1,-2];
    intra_qs(12,:)=[1,-1];

    intra_qs_vec(1:num_intra_qs,1)=intra_qs(:,1)*hex_b1(1)+intra_qs(:,2)*hex_b2(1);
    intra_qs_vec(1:num_intra_qs,2)=intra_qs(:,1)*hex_b1(2)+intra_qs(:,2)*hex_b2(2);

    num_inter_qs=12;

    interall_given_qs(1,:)=(-2*hex_shift+hex_b2);
    interall_given_qs(2,:)=-hex_shift-(hex_shift+hex_b1);
    interall_given_qs(3,:)=(-hex_shift-hex_b1)-(hex_shift-hex_b2);

    interall_given_qs(4,:)=-2*hex_shift;
    interall_given_qs(5,:)=-2*hex_shift+2*hex_b2;
    interall_given_qs(6,:)=-2*hex_shift-2*hex_b1;

    interall_given_qs(7,:)=(-2*hex_shift+hex_b2)+hex_b1;
    interall_given_qs(8,:)=(-2*hex_shift+hex_b2)+hex_b1+hex_b2;

    interall_given_qs(9,:)=(-2*hex_shift+hex_b2)-hex_b1+hex_b2;
    interall_given_qs(10,:)=(-2*hex_shift+hex_b2)-2*hex_b1;

    interall_given_qs(11,:)=(-2*hex_shift+hex_b2)-(hex_b1+hex_b2)-hex_b2;
    interall_given_qs(12,:)=(-2*hex_shift+hex_b2)-(hex_b1+hex_b2)*2;

    L12_qvecs=interall_given_qs(:,1:2);


    hex_M=50;

    hex_table=zeros(2*hex_M+1);

    hex_index=0;
    hex_coor=0;
    hex_cut=8*HEX_BLEN;%3.21*HEX_BLEN;

    ind=1;
    for ind1=(-hex_M):hex_M
        for ind2=(-hex_M):hex_M
            vec=hex_b1*ind1+hex_b2*ind2+hex_shift;

            if sqrt(dot(vec,vec))<hex_cut
                %hex_index(ind,1:2)=[ind1,ind2];
                %hex_table(ind1+hex_M+1,ind2+hex_M+1)=ind;
                hex_coor(ind,1:2)=vec(1:2);

                ind=ind+1;
            end

        end
    end

    moire_k_vec1=hex_b1;
    moire_k_vec2=hex_b2;
    moire_k_vec1(3)=0;
    moire_k_vec2(3)=0;
    moire_k_vec3=[0,0,1];

    vvv=abs(dot(moire_k_vec1,cross(moire_k_vec2,moire_k_vec3)));
    moire_L_x1=2*pi*cross(moire_k_vec2,moire_k_vec3)/vvv;
    moire_L_x2=2*pi*cross(moire_k_vec3,moire_k_vec1)/vvv;
    moire_L_x3=2*pi*cross(moire_k_vec1,moire_k_vec2)/vvv;


    % checks out, 13nm
    moire_L_x2*2.46/10;

    num_hex=ind-1;

    %hex_delta=KTH*[0,-1];
    % 2 Layers
    tot_dim=num_hex*unit_dim*2;

    hex_all=zeros(num_hex,2);
    hex_all(1:num_hex,:)=hex_coor(:,:);


    hex_all_L1=hex_all;
    hex_all_L2=-hex_all;

    vecthres=1E-5;
    % create connection matrices
    connect_Mat_L1=zeros(num_hex,num_hex,num_intra_qs);
    connect_Mat_L2=zeros(num_hex,num_hex,num_intra_qs);

    tot_num_G12=num_inter_qs;
    connect_Mat_L12=zeros(num_hex,num_hex,tot_num_G12);


    for ind1=1:num_hex
        pvec1_L1=hex_all_L1(ind1,:);
        pvec1_L2=hex_all_L2(ind1,:);

        for ind2=1:num_hex
            pvec2_L1=hex_all_L1(ind2,:);
            pvec2_L2=hex_all_L2(ind2,:);



            for indty=1:num_inter_qs
                tmp_qvec=L12_qvecs(indty,1:2);
                dvec_L12=pvec2_L2-pvec1_L1-tmp_qvec;
                if sqrt(dot(dvec_L12,dvec_L12))<vecthres
                    connect_Mat_L12(ind2,ind1,indty)=1;
                end
            end


            for indty=1:num_intra_qs
                tmp_qvec=intra_qs_vec(indty,1:2);

                dvec_L1=pvec2_L1-pvec1_L1-tmp_qvec;
                dvec_L2=pvec2_L2-pvec1_L2-tmp_qvec;
                % dvec=pvec2-pvec1-(tmp_set(1)*sc_b1(1:2)+tmp_set(2)*sc_b2(1:2));

                if sqrt(dot(dvec_L1,dvec_L1))<vecthres
                    connect_Mat_L1(ind1,ind2,indty)=1;
                end
                if sqrt(dot(dvec_L2,dvec_L2))<vecthres
                    connect_Mat_L2(ind1,ind2,indty)=1;
                end
            end


        end
    end


    %

    % expfac=exp(-1*i*rot_theta);
    % 
    %  Layer1_ham=@(kknow) [0,Dirac_v*(-i*kknow(1)+kknow(2));Dirac_v*(i*kknow(1)+kknow(2)),0];
    %  Layer2_ham=@(kknow) [0,expfac*Dirac_v*(-i*kknow(1)+kknow(2));conj(expfac)*Dirac_v*(i*kknow(1)+kknow(2)),0];
    % 
    expfac1=exp(1*i*rot_theta/2);
     expfac2=exp(-1*i*rot_theta/2);

     Layer1_ham=@(kknow) [0,expfac1*Dirac_v*(-i*kknow(1)+kknow(2));conj(expfac1)*Dirac_v*(i*kknow(1)+kknow(2)),0];
     Layer2_ham=@(kknow) [0,expfac2*Dirac_v*(-i*kknow(1)+kknow(2));conj(expfac2)*Dirac_v*(i*kknow(1)+kknow(2)),0];



    all_index_L1=reshape(1:(num_hex*unit_dim),unit_dim,num_hex);
    all_index_L2=reshape((num_hex*unit_dim+1):(2*num_hex*unit_dim),unit_dim,num_hex);
    all_index=reshape(1:tot_dim,unit_dim,tot_dim/unit_dim);

    % indL1=1:(tot_dim/2);
    % indL2=(tot_dim/2+1):tot_dim;

    % K G M K

    % K
    kk1a=+(hex_b1-hex_b2)/3;
    kk1b=+(hex_b1+2*hex_b2)/3;
    kk1c=+(-2*hex_b1-hex_b2)/3;

    % Ks'
    kk2a=-(hex_b1-hex_b2)/3;
    kk2b=-(hex_b1+2*hex_b2)/3;
    kk2c=-(-2*hex_b1-hex_b2)/3;

    kk3=0*hex_b1;
    kk4=hex_b1/2;

    scan_klist1=[kk1a;kk3;kk4;kk1a];
    scan_klist2=-[kk1a;kk3;kk4;kk1a];
    scan_klist1(:,3)=0;
    scan_klist2(:,3)=0;


    knum=51;

    [ all_kpts1, scale_axis1] = generate_k_line( knum, scan_klist1 );
    [ all_kpts2, scale_axis2] = generate_k_line( knum, scan_klist2 );
    knum_tot=size(all_kpts1);
    knum_tot=knum_tot(1);

    Hmat_strain_L1=zeros(tot_dim);
    Hmat_strain_L2=zeros(tot_dim);

    Hmat_strain_L1_kplus=zeros(tot_dim);
    Hmat_strain_L2_kplus=zeros(tot_dim);
    Hmat_strain_L1_kminus=zeros(tot_dim);
    Hmat_strain_L2_kminus=zeros(tot_dim);

    for indq=1:num_intra_qs
        indq
        T_mat_L1=squeeze(All_Eff_intra_bot(:,:,indq));
        T_mat_L2=squeeze(All_Eff_intra_top(:,:,indq));
        Hmat_strain_L1(all_index_L1(:),all_index_L1(:))=Hmat_strain_L1(all_index_L1(:),all_index_L1(:))+kron(squeeze(connect_Mat_L1(:,:,indq)),T_mat_L1);
        Hmat_strain_L2(all_index_L2(:),all_index_L2(:))=Hmat_strain_L2(all_index_L2(:),all_index_L2(:))+kron(squeeze(connect_Mat_L2(:,:,indq)),T_mat_L2);

        T_mat_L1_kplus=squeeze(All_Eff_intra_bot_kplus(:,:,indq));
        T_mat_L2_kplus=squeeze(All_Eff_intra_top_kplus(:,:,indq));
        Hmat_strain_L1_kplus(all_index_L1(:),all_index_L1(:))=Hmat_strain_L1_kplus(all_index_L1(:),all_index_L1(:))+kron(squeeze(connect_Mat_L1(:,:,indq)),T_mat_L1_kplus);
        Hmat_strain_L2_kplus(all_index_L2(:),all_index_L2(:))=Hmat_strain_L2_kplus(all_index_L2(:),all_index_L2(:))+kron(squeeze(connect_Mat_L2(:,:,indq)),T_mat_L2_kplus);

        T_mat_L1_kminus=squeeze(All_Eff_intra_bot_kminus(:,:,indq));
        T_mat_L2_kminus=squeeze(All_Eff_intra_top_kminus(:,:,indq));
        Hmat_strain_L1_kminus(all_index_L1(:),all_index_L1(:))=Hmat_strain_L1_kminus(all_index_L1(:),all_index_L1(:))+kron(squeeze(connect_Mat_L1(:,:,indq)),T_mat_L1_kminus);
        Hmat_strain_L2_kminus(all_index_L2(:),all_index_L2(:))=Hmat_strain_L2_kminus(all_index_L2(:),all_index_L2(:))+kron(squeeze(connect_Mat_L2(:,:,indq)),T_mat_L2_kminus);



    end

    Hmat_strain_L1=(Hmat_strain_L1+Hmat_strain_L1')/2;
    Hmat_strain_L2=(Hmat_strain_L2+Hmat_strain_L2')/2;


    Hmat_inter=zeros(tot_dim);
    Hmat_inter_kplus0=zeros(tot_dim/2);
    Hmat_inter_kminus0=zeros(tot_dim/2);

    for indq=1:num_inter_qs
        indq
        T_tmp=squeeze(All_Eff_inter(:,:,indq));
        Hmat_inter(all_index_L2(:),all_index_L1(:))=Hmat_inter(all_index_L2(:),all_index_L1(:))+kron(squeeze(connect_Mat_L12(:,:,indq)),T_tmp);

        T_tmp_kplus=squeeze(All_Eff_inter_kplus(:,:,indq))*lattice_a;
        Hmat_inter_kplus0(:,:)=Hmat_inter_kplus0(:,:)+kron(squeeze(connect_Mat_L12(:,:,indq)),T_tmp_kplus);

        T_tmp_kminus=squeeze(All_Eff_inter_kminus(:,:,indq))*lattice_a;
        Hmat_inter_kminus0(:,:)=Hmat_inter_kminus0(:,:)+kron(squeeze(connect_Mat_L12(:,:,indq)),T_tmp_kminus);


    end
    Hmat_inter=Hmat_inter+Hmat_inter';


    %

    strain_fac=1.0;
    inter_fac=1.0;
    allbands1=zeros(tot_dim,knum_tot);
    allbands2=zeros(tot_dim,knum_tot);

    for indk=1:knum_tot
        indk
        know=all_kpts1(indk,1:2);
        shift_klist_L1=zeros(num_hex,2);
        shift_klist_L2=zeros(num_hex,2);
        shift_klist_L1(:,1)=hex_all_L1(:,1)+know(1);
        shift_klist_L1(:,2)=hex_all_L1(:,2)+know(2);
        shift_klist_L2(:,1)=hex_all_L2(:,1)+know(1);
        shift_klist_L2(:,2)=hex_all_L2(:,2)+know(2);

        kplus_L1=kron(diag(shift_klist_L1(:,1)+i*shift_klist_L1(:,2)),eye(2));
        kminus_L1=kron(diag(shift_klist_L1(:,1)-i*shift_klist_L1(:,2)),eye(2));
        kplus_L2=kron(diag(shift_klist_L2(:,1)+i*shift_klist_L2(:,2)),eye(2));
        kminus_L2=kron(diag(shift_klist_L2(:,1)-i*shift_klist_L2(:,2)),eye(2));

        Hmat_inter_kk=zeros(tot_dim);

        Hmat_inter_kk(all_index_L2(:),all_index_L1(:))=Hmat_inter_kk(all_index_L2(:),all_index_L1(:))+Hmat_inter_kplus0*kplus_L1+Hmat_inter_kminus0*kminus_L1;

        Hmat_inter_kk=Hmat_inter_kk+Hmat_inter_kk';


        Hmat = zeros(tot_dim);
        % construct the hamiltonian
        Hmat=Hmat+strain_fac*(Hmat_strain_L1+Hmat_strain_L2)+inter_fac*(Hmat_inter)+Hmat_inter_kk*1;
        %Hmat=(Hmat+Hmat')/2;
        % intralayer part
        % layer 1
        for indh=1:num_hex
            Hmat(all_index_L1(:,indh),all_index_L1(:,indh))=Hmat(all_index_L1(:,indh),all_index_L1(:,indh))+Layer1_ham(shift_klist_L1(indh,:));
            Hmat(all_index_L2(:,indh),all_index_L2(:,indh))=Hmat(all_index_L2(:,indh),all_index_L2(:,indh))+Layer2_ham(shift_klist_L2(indh,:));
        end
        [newbands,indband]=sort(real(eig(Hmat)),'ascend');
        allbands1(:,indk)=newbands;


        know=all_kpts2(indk,1:2);
        shift_klist_L1=zeros(num_hex,2);
        shift_klist_L2=zeros(num_hex,2);
        shift_klist_L1(:,1)=hex_all_L1(:,1)+know(1);
        shift_klist_L1(:,2)=hex_all_L1(:,2)+know(2);
        shift_klist_L2(:,1)=hex_all_L2(:,1)+know(1);
        shift_klist_L2(:,2)=hex_all_L2(:,2)+know(2);

        kplus_L1=kron(diag(shift_klist_L1(:,1)+i*shift_klist_L1(:,2)),eye(2));
        kminus_L1=kron(diag(shift_klist_L1(:,1)-i*shift_klist_L1(:,2)),eye(2));
        kplus_L2=kron(diag(shift_klist_L2(:,1)+i*shift_klist_L2(:,2)),eye(2));
        kminus_L2=kron(diag(shift_klist_L2(:,1)-i*shift_klist_L2(:,2)),eye(2));

        Hmat_inter_kk=zeros(tot_dim);

        Hmat_inter_kk(all_index_L2(:),all_index_L1(:))=Hmat_inter_kk(all_index_L2(:),all_index_L1(:))+Hmat_inter_kplus0*kplus_L1+Hmat_inter_kminus0*kminus_L1;

        Hmat_inter_kk=Hmat_inter_kk+Hmat_inter_kk';


        Hmat = zeros(tot_dim);
        % construct the hamiltonian
        Hmat=Hmat+strain_fac*(Hmat_strain_L1+Hmat_strain_L2)+inter_fac*(Hmat_inter)+Hmat_inter_kk*1;
        %Hmat=(Hmat+Hmat')/2;
        % intralayer part
        % layer 1
        for indh=1:num_hex
            Hmat(all_index_L1(:,indh),all_index_L1(:,indh))=Hmat(all_index_L1(:,indh),all_index_L1(:,indh))+Layer1_ham(shift_klist_L1(indh,:));
            Hmat(all_index_L2(:,indh),all_index_L2(:,indh))=Hmat(all_index_L2(:,indh),all_index_L2(:,indh))+Layer2_ham(shift_klist_L2(indh,:));
        end
        [newbands,indband]=sort(real(eig(Hmat)),'ascend');



        allbands2(:,indk)=newbands;


    end
    
    fprintf(['saving file: ' bands_filename_here '.mat \n']);
    save([bands_data_dir '/' bands_filename_here],'scale_axis1','all_kpts1','allbands1','scale_axis2','all_kpts2','allbands2');    

    
end

%%
%{


  benchmark_bands=load('Eff_kp_benchmark_30_29_xy.mat','allbands','scale_axis');
%  benchmark_bands=load('Eff_kp_benchmark_30_29_xy.mat','allbands','scale_axis');



num_bb=size(benchmark_bands.allbands,1)/2;
num_bbscan=100;

EEmid=mean(benchmark_bands.allbands(num_bb:(num_bb+1),1));





figure(1);


pltmax=0.3;

subplot(1,3,1)



centmp1=mean(allbands1((tot_dim/2):(tot_dim/2+1),1));
centmp2=mean(allbands2((tot_dim/2):(tot_dim/2+1),1));

plot(scale_axis1,allbands1'-centmp1,'r','LineWidth',2)
hold on;
 plot(scale_axis2,allbands2'-centmp2,'r','LineWidth',2)
hold on;

% plot(benchmark_bands.scale_axis,benchmark_bands.allbands((num_bb-num_bbscan):(num_bb+num_bbscan),:)-EEmid,'k','LineWidth',2)


line([scale_axis1(knum),scale_axis1(knum)],[-100,100],'Color','k');

line([scale_axis1(knum*2),scale_axis1(knum*2)],[-100,100],'Color','k');
set(gca,'XTick',[0,scale_axis1(knum),scale_axis1(2*knum),1]);
set(gca,'XTickLabel',{'K','\Gamma','M','K'});
title('Effective KP')
% plot(scale_axis2,allbands2','b','LineWidth',2)
box on;

hold off;

axis([-inf,inf,-pltmax,pltmax])
% ylabel('Energy (eV)')
set(gca,'FontSize',20)

subplot(1,3,2)


centmp1=mean(allbands1((tot_dim/2):(tot_dim/2+1),1));
centmp2=mean(allbands2((tot_dim/2):(tot_dim/2+1),1));

% plot(scale_axis1,allbands1'-centmp1,'r','LineWidth',2)
hold on;
 % plot(scale_axis2,allbands2'-centmp2,'r','LineWidth',2)
hold on;

 plot(benchmark_bands.scale_axis,benchmark_bands.allbands((num_bb-num_bbscan):(num_bb+num_bbscan),:)-EEmid,'k','LineWidth',2)


line([scale_axis1(knum),scale_axis1(knum)],[-100,100],'Color','k');
title('Full TBH')
line([scale_axis1(knum*2),scale_axis1(knum*2)],[-100,100],'Color','k');
set(gca,'XTick',[0,scale_axis1(knum),scale_axis1(2*knum),1]);
set(gca,'XTickLabel',{'K','\Gamma','M','K'});
box on;
% plot(scale_axis2,allbands2','b','LineWidth',2)


hold off;

axis([-inf,inf,-pltmax,pltmax])
% ylabel('Energy (eV)')
set(gca,'FontSize',20)

subplot(1,3,3)


centmp1=mean(allbands1((tot_dim/2):(tot_dim/2+1),1));
centmp2=mean(allbands2((tot_dim/2):(tot_dim/2+1),1));

plot(scale_axis1,allbands1'-centmp1,'r','LineWidth',2)
hold on;
 plot(scale_axis2,allbands2'-centmp2,'r','LineWidth',2)
hold on;

 plot(benchmark_bands.scale_axis,benchmark_bands.allbands((num_bb-num_bbscan):(num_bb+num_bbscan),:)-EEmid,'k','LineWidth',2)


line([scale_axis1(knum),scale_axis1(knum)],[-100,100],'Color','k');

line([scale_axis1(knum*2),scale_axis1(knum*2)],[-100,100],'Color','k');
set(gca,'XTick',[0,scale_axis1(knum),scale_axis1(2*knum),1]);
set(gca,'XTickLabel',{'K','\Gamma','M','K'});
title('Eff. KP vs. Full TBH')
% plot(scale_axis2,allbands2','b','LineWidth',2)


hold off;

axis([-inf,inf,-pltmax,pltmax])
% ylabel('Energy (eV)')
set(gca,'FontSize',20)
%}