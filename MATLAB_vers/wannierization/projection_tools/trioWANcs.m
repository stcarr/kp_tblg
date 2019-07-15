function [filename] = trioWANcs(tar_theta,tar_folder,numk,numxy,side)

    if side == 0 % top 3 bands
        side_text  = 'top';
        sym_sign   = -1; % p_z type
        start_band = 2; % 2 above nb/2
    elseif side == 1 % lower 3 bands
        side_text  = 'bot';
        sym_sign   = 1; % s type
        start_band = -3; % 3 below nb/2
    end
    
    filename = [tar_folder '/trioWANcs' side_text '_data_' strrep(num2str(tar_theta,'%.2f'),'.','p') '_' num2str(numk) 'k_' num2str(numxy) 'xy.mat'];

    % the twist angle converted into radian
    % this MacDonald effective theory has magic angle ~1.1 degree
    rot_theta=tar_theta*(pi/180);
    %rot_theta=1.5/180*pi;
    %rot_theta=0.63/180*pi;

    % cutoff in momentum space, which would affect accuracies for the solution
    hex_cut=4.51;


    lattice_a=1.42*sqrt(3);

    KD=4*pi/3/lattice_a;

    KTH=2*KD*sin(rot_theta/2);

    bilayer_w0=0.0;
    bilayer_w1=0.11;
    fermivf=1*2.7*sqrt(3)/2*lattice_a;

    % basic variables
    sigx=[0,1;1,0];
    sigy=[0,-i;i,0];
    sigz=[1,0;0,-1];



    ph1=exp(i*2*pi/3);
    ph2=exp(-i*2*pi/3);

    TT1_0=[bilayer_w0,bilayer_w1;bilayer_w1,bilayer_w0];
    TT2_0=[bilayer_w0,bilayer_w1*ph2;bilayer_w1*ph1,bilayer_w0];
    TT3_0=[bilayer_w0,bilayer_w1*ph1;bilayer_w1*ph2,bilayer_w0];

    % hamiltonian
    % no theta dependence yet
    gen_dirac_ham0=@(qvec,dummy_theta) fermivf*(qvec(1)*sigy+qvec(2)*sigx);

    gen_dirac_ham=@(qvec,angle_theta) fermivf*(qvec(1)*sigy+qvec(2)*sigx)*diag([exp(-i*angle_theta),exp(i*angle_theta)]);

    % define T matrix


    % create hex table

    hex_a1=[sqrt(3)/2,1/2];
    hex_a2=[-sqrt(3)/2,1/2];
    %hex_a3=[0,1];

    hex_vertshift=(hex_a1+2*hex_a2)/3;

    hex_M=40;
    hex_table=zeros(2*hex_M+1);
    hex_index=0;
    hex_coor=0;
    ind=1;
    for ind1=(-hex_M):hex_M
        for ind2=(-hex_M):hex_M
            vec=hex_a1*ind1+hex_a2*ind2+hex_vertshift;

            if sqrt(dot(vec,vec))<hex_cut
                hex_index(ind,1:2)=[ind1,ind2];
                hex_table(ind1+hex_M+1,ind2+hex_M+1)=ind;
                hex_coor(ind,1:2)=KTH*sqrt(3)*vec;

                ind=ind+1;
            end

        end
    end

    moire_k_vec1=KTH*sqrt(3)*[hex_a1,0];
    moire_k_vec2=KTH*sqrt(3)*[hex_a2,0];
    moire_k_vec3=[0,0,1];

    vvv=abs(dot(moire_k_vec1,cross(moire_k_vec2,moire_k_vec3)));
    moire_L_x1=2*pi*cross(moire_k_vec2,moire_k_vec3)/vvv;
    moire_L_x2=2*pi*cross(moire_k_vec3,moire_k_vec1)/vvv;
    moire_L_x3=2*pi*cross(moire_k_vec1,moire_k_vec2)/vvv;

    % checks out, 13nm
    %sqrt(dot(moire_L_x2,moire_L_x2))/10


    num_hex=ind-1;
    %hex_delta=KTH*[0,-1];
    tot_dim=2*num_hex*2;

    hex_all=zeros(2*num_hex,2);
    hex_all(1:num_hex,:)=hex_coor(:,:);
    hex_all((num_hex+1):2*num_hex,1:2)=-hex_coor(:,:);




    %scatter(hex_all(:,1),hex_all(:,2))
    %axis equal;




    group1=1:num_hex;
    group2=(num_hex+1):2*num_hex;

    vecthres=1E-7;
    % create connection matrices
    connect_Mat1=zeros(num_hex);
    connect_Mat2=zeros(num_hex);
    connect_Mat3=zeros(num_hex);

    for ind1=1:num_hex
        pvec1=hex_all(group1(ind1),:);

        for ind2=1:num_hex
            pvec2=hex_all(group2(ind2),:);

            dvec1=pvec2-pvec1-KTH*[1,0];
            dvec2=pvec2-pvec1-KTH*[-0.5,-0.5*sqrt(3)];
            dvec3=pvec2-pvec1-KTH*[-0.5,0.5*sqrt(3)];

            if sqrt(dot(dvec1,dvec1))<vecthres
                connect_Mat1(ind1,ind2)=1;
            end
            if sqrt(dot(dvec2,dvec2))<vecthres
                connect_Mat2(ind1,ind2)=1;
            end
            if sqrt(dot(dvec3,dvec3))<vecthres
                connect_Mat3(ind1,ind2)=1;
            end




        end
    end


    %
    % figure(1)
    % scatter(hex_all(:,1),hex_all(:,2))
    % axis equal;
    % hold on;
    % for ind1=1:num_hex
    %     pos1=hex_all(group1(ind1),:);
    %
    %     for ind2=1:num_hex
    %         pos2=hex_all(group2(ind2),:);
    %
    %         if connect_Mat1(ind2,ind1)==1
    %             line([pos1(1),pos2(1)],[pos1(2),pos2(2)],'Color','r');
    %         end
    %         if connect_Mat2(ind2,ind1)==1
    %             line([pos1(1),pos2(1)],[pos1(2),pos2(2)],'Color','g');
    %         end
    %         if connect_Mat3(ind2,ind1)==1
    %             line([pos1(1),pos2(1)],[pos1(2),pos2(2)],'Color','b');
    %         end
    %
    %     end
    % end
    % hold off;




    all_index=reshape(1:tot_dim,2,tot_dim/2);

    indL1=1:(tot_dim/2);
    indL2=(tot_dim/2+1):tot_dim;

    indL1A=indL1(1:2:(tot_dim/2));
    indL1B=indL1(2:2:(tot_dim/2));

    indL2A=indL2(1:2:(tot_dim/2));
    indL2B=indL2(2:2:(tot_dim/2));

    % build trial states

    LLmoire=sqrt(dot(moire_L_x1,moire_L_x1));
    %numxy=51;%401; passed to function
    allxy=2.0*LLmoire*linspace(-1,1,numxy);



    kagome_center=(moire_L_x1-moire_L_x2)/2;


    midpt=(numxy+1)/2;
    [mmX,mmY]=meshgrid(allxy,allxy);

    all_wfAL1=zeros(numxy,numxy,3);
    all_wfBL1=zeros(numxy,numxy,3);
    all_wfAL2=zeros(numxy,numxy,3);
    all_wfBL2=zeros(numxy,numxy,3);


    wfAL1=0*mmX;
    wfAL2=0*mmX;
    wfBL1=0*mmX;
    wfBL2=0*mmX;

    sigmax=35;
    sigmay=35;

    sigmax=28;
    sigmay=28;

    tmpwf=exp(-((mmX-kagome_center(1)).^2/2/(sigmax^2)+(mmY-kagome_center(2)).^2/2/(sigmay^2)));
    tmpwf=tmpwf/sqrt(sum(abs(tmpwf(:)).^2));

    shift_L1A=exp(1*i*2*pi/3);
    shift_L1B=exp(1*i*2*pi/3);
    shift_L2A=exp(-1*i*2*pi/3);
    shift_L2B=exp(-1*i*2*pi/3);

    wfAL1=tmpwf*shift_L1A;

    % mirror
    wfBL2=sym_sign*tmpwf*shift_L2B;

    % C2z T
    wfBL1=tmpwf*shift_L1B;

    % mirror from above
    wfAL2=sym_sign*tmpwf*shift_L2A;

    all_wfAL1(:,:,1)=wfAL1(:,:);
    all_wfBL1(:,:,1)=wfBL1(:,:);
    all_wfAL2(:,:,1)=wfAL2(:,:);
    all_wfBL2(:,:,1)=wfBL2(:,:);


    kagome_center=(moire_L_x2)/2;

    wfAL1=0*mmX;
    wfAL2=0*mmX;
    wfBL1=0*mmX;
    wfBL2=0*mmX;

    sigmax=35;
    sigmay=35;

    sigmax=28;
    sigmay=28;

     tmpwf=exp(-((mmX-kagome_center(1)).^2/2/(sigmax^2)+(mmY-kagome_center(2)).^2/2/(sigmay^2)));
     tmpwf=tmpwf/sqrt(sum(abs(tmpwf(:)).^2));


     shift_L1A=exp(1*i*2*pi/3)*exp(-1*i*2*pi/3);
     shift_L1B=exp(1*i*2*pi/3)*exp(1*i*2*pi/3);
     shift_L2A=exp(-1*i*2*pi/3)*exp(-1*i*2*pi/3);
     shift_L2B=exp(-1*i*2*pi/3)*exp(1*i*2*pi/3);

     wfAL1=tmpwf*shift_L1A;

     % mirror
     wfBL2=sym_sign*tmpwf*shift_L2B;

     % C2z T
     wfBL1=tmpwf*shift_L1B;

     % mirror from above 
     wfAL2=sym_sign*tmpwf*shift_L2A;

    all_wfAL1(:,:,2)=wfAL1(:,:);
    all_wfBL1(:,:,2)=wfBL1(:,:);
    all_wfAL2(:,:,2)=wfAL2(:,:);
    all_wfBL2(:,:,2)=wfBL2(:,:);


    kagome_center=(-moire_L_x1)/2;

    wfAL1=0*mmX;
    wfAL2=0*mmX;
    wfBL1=0*mmX;
    wfBL2=0*mmX;

    sigmax=35;
    sigmay=35;

    sigmax=28;
    sigmay=28;

     tmpwf=exp(-((mmX-kagome_center(1)).^2/2/(sigmax^2)+(mmY-kagome_center(2)).^2/2/(sigmay^2)));
     tmpwf=tmpwf/sqrt(sum(abs(tmpwf(:)).^2));

     shift_L1A=exp(1*i*2*pi/3)*exp(1*i*2*pi/3);
     shift_L1B=exp(1*i*2*pi/3)*exp(-1*i*2*pi/3);
     shift_L2A=exp(-1*i*2*pi/3)*exp(1*i*2*pi/3);
     shift_L2B=exp(-1*i*2*pi/3)*exp(-1*i*2*pi/3);

     wfAL1=tmpwf*shift_L1A;

     % mirror
     wfBL2=sym_sign*tmpwf*shift_L2B;

     % C2z T
     wfBL1=tmpwf*shift_L1B;

     % mirror from above 
     wfAL2=sym_sign*tmpwf*shift_L2A;


     all_wfAL1(:,:,3)=wfAL1(:,:);
    all_wfBL1(:,:,3)=wfBL1(:,:);
    all_wfAL2(:,:,3)=wfAL2(:,:);
    all_wfBL2(:,:,3)=wfBL2(:,:);



    for ind=1:3
        tmp_AL1=squeeze(all_wfAL1(:,:,ind));
        tmp_AL2=squeeze(all_wfAL2(:,:,ind));
        tmp_BL1=squeeze(all_wfBL1(:,:,ind));
        tmp_BL2=squeeze(all_wfBL2(:,:,ind));

        tmp_sum=sqrt(sum(abs(tmp_AL1(:)).^2)+sum(abs(tmp_BL1(:)).^2)+sum(abs(tmp_AL2(:)).^2)+sum(abs(tmp_BL2(:)).^2));

        all_wfAL1(:,:,ind)=tmp_AL1/tmp_sum;
        all_wfAL2(:,:,ind)=tmp_AL2/tmp_sum;
        all_wfBL1(:,:,ind)=tmp_BL1/tmp_sum;
        all_wfBL2(:,:,ind)=tmp_BL2/tmp_sum;

    end

    knum=numk; % =9, passed to function

    allk12=linspace(0,1,knum+1);
    allk12=allk12(1:knum);
    knum_tot=knum^2;

    [meshk1,meshk2]=meshgrid(allk12,allk12);
    all_kpts=zeros(knum_tot,3);
    all_kpts(:,1)=moire_k_vec1(1)*meshk1(:)+moire_k_vec2(1)*meshk2(:);
    all_kpts(:,2)=moire_k_vec1(2)*meshk1(:)+moire_k_vec2(2)*meshk2(:);

    all_new_basis=zeros(tot_dim,3,knum_tot);
    all_new_hmat=zeros(3,3,knum_tot);


    all_index=reshape(1:tot_dim,2,tot_dim/2);




    fff=tot_dim/2;

    % sel_bandind=(fff-3):(fff+1);
    sel_bandind=(fff+start_band):(fff+start_band+2);

    % save('Ten_trial_WF','mmX','mmY','all_wfAL1','all_wfAL2','all_wfBL1','all_wfBL2')

    %load('Ten_WAN/WAN_Total_data.mat');

    all_vecs=zeros(tot_dim,length(sel_bandind),knum_tot);

    all_ssdiag=zeros(3,knum_tot);

    allbands=zeros(tot_dim,knum_tot);

    allbands_eff=zeros(3,knum_tot);

    testfac=0;
    for indk=1:knum_tot
        indk

        %tic

        know=all_kpts(indk,1:2);

        shift_klist(:,1)=hex_all(:,1)+know(1);
        shift_klist(:,2)=hex_all(:,2)+know(2);

        Hmat = zeros(tot_dim);

        % construct the hamiltonian

        Hmat(indL2,indL1)=kron(connect_Mat1,TT1_0)+kron(connect_Mat2,TT2_0)+kron(connect_Mat3,TT3_0);


        % interlayer part

        Hmat=Hmat+Hmat';


        % layer 1
        for indh=1:num_hex
            Hmat(all_index(:,indh),all_index(:,indh))=gen_dirac_ham(shift_klist(indh,:),testfac*rot_theta/2);

        end

        % layer 2
        for indh=(num_hex+1):2*num_hex
            Hmat(all_index(:,indh),all_index(:,indh))=gen_dirac_ham(shift_klist(indh,:),-testfac*rot_theta/2);

        end


        [V,D]=eig(Hmat);


        [allbands(:,indk),tmpind]=sort(real(diag(D)),'ascend');
        V(:,:)=V(:,tmpind);

        all_vecs(:,:,indk)=V(:,sel_bandind);


        PPk_proj=V(:,sel_bandind)*(V(:,sel_bandind)');



        % trial WF
        kktmpL1=shift_klist(1:(tot_dim/4),:);
        kktmpL2=shift_klist((tot_dim/4+1):(tot_dim/2),:);

        trial_wf_bloch=zeros(tot_dim,3);

        for indw=1:3

            wfAL1=squeeze(all_wfAL1(:,:,indw));
            wfBL1=squeeze(all_wfBL1(:,:,indw));
            wfAL2=squeeze(all_wfAL2(:,:,indw));
            wfBL2=squeeze(all_wfBL2(:,:,indw));


            for indt=1:(tot_dim/4)
                kkkL1=kktmpL1(indt,:);
                kkkL2=kktmpL2(indt,:);

                tmp1=exp(-i*(kkkL1(1)*mmX+kkkL1(2)*mmY));
                tmp2=exp(-i*(kkkL2(1)*mmX+kkkL2(2)*mmY));

                tmpo1=sum(wfAL1(:).*tmp1(:));
                tmpo2=sum(wfBL1(:).*tmp1(:));
                tmpo3=sum(wfAL2(:).*tmp2(:));
                tmpo4=sum(wfBL2(:).*tmp2(:));

                trial_wf_bloch(indL1A(indt),indw)=tmpo1;
                trial_wf_bloch(indL1B(indt),indw)=tmpo2;
                trial_wf_bloch(indL2A(indt),indw)=tmpo3;
                trial_wf_bloch(indL2B(indt),indw)=tmpo4;


            end
            trial_wf_bloch(:,indw)=trial_wf_bloch(:,indw)/sqrt(trial_wf_bloch(:,indw)'*trial_wf_bloch(:,indw));

        end

        trial_wf_proj=PPk_proj*trial_wf_bloch;

        %diag(trial_wf_proj'*trial_wf_proj);

        %     for indw=1:10
        %         tmpwf=trial_wf_proj(:,indw);
        %         trial_wf_proj(:,indw)=trial_wf_proj(:,indw)/sqrt(abs(tmpwf'*tmpwf));
        %
        %     end
        Band_state=V(:,sel_bandind);

        Amat=Band_state'*trial_wf_proj;


        [Umat,Smat,Vmat] = svd(Amat);

        ssdiag=(diag(Smat'*Smat));
        all_ssdiag(:,indk)=sort(diag(Smat));
        %sort(diag(Smat))
        Sinv_root=Vmat*diag(1./sqrt(ssdiag))*Vmat';

        %Amat-Umat*Smat*Vmat';
        %Amat'*Amat-Vmat*Smat*Smat*Vmat';

        %Sinv_root*Vmat*Smat*Smat*Vmat';

        %Amat'*Amat*Sinv_root*Sinv_root';

        rotate_state_wf=trial_wf_proj*Sinv_root;


        all_new_basis(:,:,indk)=rotate_state_wf(:,:);

        %mesh(abs(ppp'*ppp))


        Heff=rotate_state_wf'*Hmat*rotate_state_wf;

        all_new_hmat(:,:,indk)=Heff;

        allbands_eff(:,indk)=sort(real(eig(Heff)));

        %toc

    end


    %{

    kindnow=10;

    Vstate=squeeze(all_new_basis(:,:,kindnow));

    Vstate'*Vstate-eye(3)


    %}
    
    %{
    figure(1)
    clf
    scatter3(all_kpts(:,1),all_kpts(:,2),squeeze(all_ssdiag(1,:)));
    hold on;


    scatter3(all_kpts(:,1),all_kpts(:,2),squeeze(all_ssdiag(3,:)));

    hold off;


    axis([-inf,inf,-inf,inf,0,inf])
    min(all_ssdiag(1,:))
    max(all_ssdiag(1,:))

    % 
    % 
    % figure(2);
    % 
    % scatter3(all_kpts(:,1),all_kpts(:,2),imag(all_new_hmat(1,4,:)));
    % 
    % axis([-inf,inf,-inf,inf,-inf,inf]);
    %}

    mmm=size(mmX,1);
    all_wfAL1=zeros(mmm,mmm,3);
    all_wfBL1=zeros(mmm,mmm,3);
    all_wfAL2=zeros(mmm,mmm,3);
    all_wfBL2=zeros(mmm,mmm,3);


    for wannier_ind=1:3
        wannier_ind
        %mmXp=mmX;
        %mmYp=mmY;


        Proj_wf_R=zeros(numxy,numxy,4);


        for indk=1:knum_tot
            %indk

            know=all_kpts(indk,1:2);
            shift_klist(:,1)=hex_all(:,1)+know(1);
            shift_klist(:,2)=hex_all(:,2)+know(2);

            kktmpL1=shift_klist(1:(tot_dim/4),:);
            kktmpL2=shift_klist((tot_dim/4+1):(tot_dim/2),:);

            Rtmp_wf=0*mmX;


            wfnow=squeeze(all_new_basis(:,wannier_ind,indk));


            for indt=1:(tot_dim/4)
                kkkL1=kktmpL1(indt,:);
                kkkL2=kktmpL2(indt,:);

                tmp1=exp(i*(kkkL1(1)*mmX+kkkL1(2)*mmY));
                tmp2=exp(i*(kkkL2(1)*mmX+kkkL2(2)*mmY));


                Proj_wf_R(:,:,1)= Proj_wf_R(:,:,1)+wfnow(indL1A(indt))*tmp1;
                Proj_wf_R(:,:,2)= Proj_wf_R(:,:,2)+wfnow(indL1B(indt))*tmp1;

                Proj_wf_R(:,:,3)= Proj_wf_R(:,:,3)+wfnow(indL2A(indt))*tmp2;
                Proj_wf_R(:,:,4)= Proj_wf_R(:,:,4)+wfnow(indL2B(indt))*tmp2;
            end




        end

        all_wfAL1(:,:,wannier_ind)=Proj_wf_R(:,:,1);
        all_wfBL1(:,:,wannier_ind)=Proj_wf_R(:,:,2);
        all_wfAL2(:,:,wannier_ind)=Proj_wf_R(:,:,3);
        all_wfBL2(:,:,wannier_ind)=Proj_wf_R(:,:,4);



    end


    for ind=1:3
        tmp_AL1=squeeze(all_wfAL1(:,:,ind));
        tmp_AL2=squeeze(all_wfAL2(:,:,ind));
        tmp_BL1=squeeze(all_wfBL1(:,:,ind));
        tmp_BL2=squeeze(all_wfBL2(:,:,ind));

        tmp_sum=sqrt(sum(abs(tmp_AL1(:)).^2)+sum(abs(tmp_BL1(:)).^2)+sum(abs(tmp_AL2(:)).^2)+sum(abs(tmp_BL2(:)).^2));

        all_wfAL1(:,:,ind)=tmp_AL1/tmp_sum;
        all_wfAL2(:,:,ind)=tmp_AL2/tmp_sum;
        all_wfBL1(:,:,ind)=tmp_BL1/tmp_sum;
        all_wfBL2(:,:,ind)=tmp_BL2/tmp_sum;

    end


    save(filename,'mmX','mmY','all_wfAL1','all_wfAL2','all_wfBL1','all_wfBL2','moire_L_x1','moire_L_x2');

end


% old figure plotting code from script version
%{

figure(3)
clf
hold on;
wannier_ind=3;
allcaption={'L1A','L1B','L2A','L2B'};

Proj_wf_R=zeros(numxy,numxy,4);

Proj_wf_R(:,:,1)=all_wfAL1(:,:,wannier_ind);
Proj_wf_R(:,:,2)=all_wfBL1(:,:,wannier_ind);
Proj_wf_R(:,:,3)=all_wfAL2(:,:,wannier_ind);
Proj_wf_R(:,:,4)=all_wfBL2(:,:,wannier_ind);

for ind=1:4
    subplot(2,2,ind)
    mesh(mmX,mmY,abs(squeeze(Proj_wf_R(:,:,ind))));
    title(allcaption{ind})
    set(gca,'FontSize',16)
    colorbar;
    axis equal
    view([0,0,1])
end

hold off;




%%


figure(4)
clf
hold on;

allcaption={'L1A','L1B','L2A','L2B'};
for ind=1:4
    subplot(2,2,ind)
    contourf(mmX,mmY,abs(squeeze(Proj_wf_R(:,:,ind))),50,'LineStyle','none');
    title(allcaption{ind})
    set(gca,'FontSize',16)
    colorbar;
    axis equal
    
    view([0,0,1])
end

hold off;
%}



