function [wf_filename, h_filename] = eightWANfull_fromWF(wf_in_filename,tar_theta,tar_folder,numk)
    % the twist angle converted into radian
    % this MacDonald effective theory has magic angle ~1.1 degree
    
    wf_filename = [tar_folder '/eightWANnew_wf_' strrep(num2str(tar_theta,'%.2f'),'.','p') '_' num2str(numk) 'k.mat'];
    h_filename = [tar_folder '/eightWANnew_H_' strrep(num2str(tar_theta,'%.2f'),'.','p') '_' num2str(numk) 'k.mat'];


    %tar_theta = 1.05;
    rot_theta=tar_theta*pi/180;

    lattice_a=1.42*sqrt(3);

    KD=4*pi/3/lattice_a;

    KTH=2*KD*sin(rot_theta/2);

    % create hex table

    hex_a1=[sqrt(3)/2,1/2];
    hex_a2=[-sqrt(3)/2,1/2];

    hex_vertshift=(hex_a1+2*hex_a2)/3;

    moire_k_vec1=KTH*sqrt(3)*[hex_a1,0];
    moire_k_vec2=KTH*sqrt(3)*[hex_a2,0];
    moire_k_vec3=[0,0,1];

    vvv=abs(dot(moire_k_vec1,cross(moire_k_vec2,moire_k_vec3)));
    moire_L_x1=2*pi*cross(moire_k_vec2,moire_k_vec3)/vvv;
    moire_L_x2=2*pi*cross(moire_k_vec3,moire_k_vec1)/vvv;
    moire_L_x3=2*pi*cross(moire_k_vec1,moire_k_vec2)/vvv;
    
    %{
    knum=numk; % previously knum = 15, now passed to function

    allk12=linspace(0,1,knum+1);
    allk12=allk12(1:knum);
    knum_tot=knum^2;
    
    [meshk1,meshk2]=meshgrid(allk12,allk12);
    all_kpts=zeros(knum_tot,3);
    all_kpts(:,1)=moire_k_vec1(1)*meshk1(:)+moire_k_vec2(1)*meshk2(:);
    all_kpts(:,2)=moire_k_vec1(2)*meshk1(:)+moire_k_vec2(2)*meshk2(:);
    %}
    
    allk12 = ((0:(numk-1))/numk)-1/2;
    [meshk1, meshk2] = meshgrid(allk12,allk12);
    all_kpts(:,1)=moire_k_vec1(1)*meshk1(:)+moire_k_vec2(1)*meshk2(:);
    all_kpts(:,2)=moire_k_vec1(2)*meshk1(:)+moire_k_vec2(2)*meshk2(:);
    knum_tot = size(all_kpts,1);
    
    % performs a single k-point calculation
    %{
    knum = 1;
    knum_tot = 1;
    kpt_tar = all_kpts(7,:);
    all_kpts = zeros(1,3);
    all_kpts = kpt_tar;
    %}
    
    [sweep_H,hex_all] = tblg_kp_calc_ext_getH('tar_theta',tar_theta,'k_list',all_kpts);
    tot_dim = 2*size(hex_all,1);


    all_index=reshape(1:tot_dim,2,tot_dim/2);

    indL1=1:(tot_dim/2);
    indL2=(tot_dim/2+1):tot_dim;

    indL1A=indL1(1:2:(tot_dim/2));
    indL1B=indL1(2:2:(tot_dim/2));

    indL2A=indL2(1:2:(tot_dim/2));
    indL2B=indL2(2:2:(tot_dim/2));

    all_new_basis=zeros(tot_dim,8,knum_tot);
    all_new_hmat=zeros(8,8,knum_tot);


    all_index=reshape(1:tot_dim,2,tot_dim/2);


    allbands=zeros(tot_dim,knum_tot);

    allbands_eff_trio=zeros(3,knum_tot);
    allbands_eff_ext=zeros(8,knum_tot);

    fff=tot_dim/2;
    %sel_bandind=(fff-4):(fff+5);
    % sel_bandind=(fff-6):(fff+7);
    
    % 8 seems to work well for theta ~0.85
    % might need a higher count for theta > 1.3
    sel_band_count = 20;%20;
    sel_bandind=(fff-(sel_band_count-1)):(fff+sel_band_count);


    %proj_band_ind=(fff-8):(fff+9);
    %proj_band_ind=(fff-13):(fff+14);

    %~ 0.05 for 0.85 seems OK
    % ~0.08 for 1.00
    % ~0.14 for 1.30
    EEcut_point = 0.05/(0.85/tar_theta)^2.5;

    EEcut_low_top  = -1;
    EEcut_high_bot =  1;
    
    EEcut_high_top =  EEcut_point;%0.08
    
    EEcut_low_bot  = -EEcut_point;%-0.12;
    
    EEcut_sigma = .1*sqrt(EEcut_point/0.05);
    
    %save('WAN_Total_fiveband_trio_data','mmX','mmY','all_wfAL1','all_wfAL2','all_wfBL1','all_wfBL2','moire_L_x1','moire_L_x2');
    %{
    load(trioWANcs_top_filename);

    mmX_trio=mmX;
    mmY_trio=mmY;

    all_wfAL1_trio_top=all_wfAL1;
    all_wfBL1_trio_top=all_wfBL1;
    all_wfAL2_trio_top=all_wfAL2;
    all_wfBL2_trio_top=all_wfBL2;
    
    load(trioWANcs_bot_filename);

    all_wfAL1_trio_bot=all_wfAL1;
    all_wfBL1_trio_bot=all_wfBL1;
    all_wfAL2_trio_bot=all_wfAL2;
    all_wfBL2_trio_bot=all_wfBL2;
    %}
    % save('Ten_trial_WF','mmX','mmY','all_wfAL1','all_wfAL2','all_wfBL1','all_wfBL2')

    %load('WAN_Total_fiveband_data_51xy');
    %{
    load(eightWANcs_filename);
    mmX_eight=mmX;
    mmY_eight=mmY;

    all_wfAL1_eight=all_wfAL1;
    all_wfBL1_eight=all_wfBL1;
    all_wfAL2_eight=all_wfAL2;
    all_wfBL2_eight=all_wfBL2;
    %}
    
    load(wf_in_filename);
    all_wfAL1_start=all_wfAL1;
    all_wfBL1_start=all_wfBL1;
    all_wfAL2_start=all_wfAL2;
    all_wfBL2_start=all_wfBL2;
    
    % redefine mmX and mmY for the current theta
    % (in case we are loading from a nearby result).
    %allxy = 3*linspace(-1,1,3*numxy);
    %allxy = 3*linspace(-1,1,6*numxy+1);

    %[sc_grid1,sc_grid2]=meshgrid(allxy,allxy);
    
    mmX = sc_grid1*moire_L_x1(1) + sc_grid2*moire_L_x2(1);
    mmY = sc_grid1*moire_L_x1(2) + sc_grid2*moire_L_x2(2);
    mmm=size(mmX,1);


    %all_vecs=zeros(tot_dim,length(sel_bandind),knum_tot);

    all_ssdiag_trio_top=zeros(3,knum_tot);
    all_ssdiag_trio_bot=zeros(3,knum_tot);

    all_ssdiag_eight=zeros(8,knum_tot);

    testfac=0;

    fffmid=tot_dim/2;
    all_fffind_top=(1:(fffmid+1));
    all_fffind_bot=(fffmid:tot_dim);
    %all_fffind_top=(1:(fffmid-1));
    %all_fffind_bot=(fffmid+2:tot_dim);
    allbands_eff_ext8=zeros(8,knum_tot);


    ext_state_wf=zeros(tot_dim,8);
    all_ext_state_wf8=zeros(tot_dim,8,knum_tot);
    all_ext_state_wf=zeros(tot_dim,8,knum_tot);

    for indk=1:knum_tot
        fprintf("%d / %d k-points \n",indk,knum_tot);

        %tic

        know=all_kpts(indk,1:2);

        shift_klist(:,1)=hex_all(:,1)+know(1);
        shift_klist(:,2)=hex_all(:,2)+know(2);

        Hmat = sweep_H{indk};

        [V,D]=eig(Hmat);


        [allbands(:,indk),tmpind]=sort(real(diag(D)),'ascend');
        V(:,:)=V(:,tmpind);


        all_vecs(:,:,indk)=V(:,sel_bandind);


        PPk_proj=V(:,sel_bandind)*(V(:,sel_bandind)');

        % trial WF
        kktmpL1=shift_klist(1:(tot_dim/4),:);
        kktmpL2=shift_klist((tot_dim/4+1):(tot_dim/2),:);

        trial_wf_bloch=zeros(tot_dim,8);

        for indw=1:8

            wfAL1=squeeze(all_wfAL1_start(:,:,indw));
            wfBL1=squeeze(all_wfBL1_start(:,:,indw));
            wfAL2=squeeze(all_wfAL2_start(:,:,indw));
            wfBL2=squeeze(all_wfBL2_start(:,:,indw));


            for indt=1:(tot_dim/4)
                kkkL1=kktmpL1(indt,:);
                kkkL2=kktmpL2(indt,:);

                tmp1=exp(-1j*(kkkL1(1)*mmX+kkkL1(2)*mmY));
                tmp2=exp(-1j*(kkkL2(1)*mmX+kkkL2(2)*mmY));

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

        allbands_eff_ext8(:,indk)=sort(real(eig(Heff)));

        %toc

    end
    
    % build wannier states

    %LLmoire=sqrt(dot(moire_L_x1,moire_L_x1));
    %numxy=51;%501;
    %allxy=2.0*LLmoire*linspace(-1,1,numxy);

    %midpt=(numxy+1)/2;
    %[mmX,mmY]=meshgrid(allxy,allxy);

    %all_new_basis = all_ext_state_wf8;

    LLmoire=sqrt(dot(moire_L_x1,moire_L_x1));
    %numxy=51;%401; passed to function
    %allxy=2.0*LLmoire*linspace(-1,1,numxy);
    %{
    allxy = 3*linspace(-1,1,3*numxy);
    

    [sc_grid1,sc_grid2]=meshgrid(allxy,allxy);
    
    mmX = sc_grid1*moire_L_x1(1) + sc_grid2*moire_L_x2(1);
    mmY = sc_grid1*moire_L_x1(2) + sc_grid2*moire_L_x2(2);
    mmm=size(mmX,1);
    %}
    
    all_wfAL1=zeros(mmm,mmm,8);
    all_wfBL1=zeros(mmm,mmm,8);
    all_wfAL2=zeros(mmm,mmm,8);
    all_wfBL2=zeros(mmm,mmm,8);


    for wannier_ind=1:8
        fprintf('starting %d/%d wannier constructions \n',wannier_ind,8);
        %mmXp=mmX;
        %mmYp=mmY;


        Proj_wf_R=zeros(mmm,mmm,4);


        for indk=1:knum_tot
            %indk

            know=all_kpts(indk,1:2);
            shift_klist(:,1)=hex_all(:,1)+know(1);
            shift_klist(:,2)=hex_all(:,2)+know(2);

            kktmpL1=shift_klist(1:(tot_dim/4),:);
            kktmpL2=shift_klist((tot_dim/4+1):(tot_dim/2),:);

            wfnow=squeeze(all_new_basis(:,wannier_ind,indk));

            for indt=1:(tot_dim/4)
                kkkL1=kktmpL1(indt,:);
                kkkL2=kktmpL2(indt,:);

                tmp1=exp(1j*(kkkL1(1)*mmX+kkkL1(2)*mmY));
                tmp2=exp(1j*(kkkL2(1)*mmX+kkkL2(2)*mmY));
                % there are 8 psi_nk(r)_L1A's (one for each wannier)
                %
                % psi_nk(r)_L1A = sum(wfnow(:).*exp(1j*(internal_kpts(:,1:2)*r))
                % I think wfnow may be what is needed to sum over for M_{mn}
                % to give accurate overlap matricies.
                
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

    for ind=1:8
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


    %{

    figure(3)
    clf
    hold on;
    wannier_ind=1;
    allcaption={'L1A','L1B','L2A','L2B'};

    Proj_wf_R=zeros(mmm,mmm,4);

    Proj_wf_R(:,:,1)=all_wfAL1(:,:,wannier_ind);
    Proj_wf_R(:,:,2)=all_wfBL1(:,:,wannier_ind);
    Proj_wf_R(:,:,3)=all_wfAL2(:,:,wannier_ind);
    Proj_wf_R(:,:,4)=all_wfBL3(:,:,wannier_ind);

    z_max = max(abs(Proj_wf_R(:)));


    for ind=1:4
        subplot(2,2,ind)
        hold on
        surf(mmX,mmY,abs(squeeze(Proj_wf_R(:,:,ind))),'EdgeColor','none');
        title(allcaption{ind})
        set(gca,'FontSize',16)
        caxis([0 z_max])
        colorbar;
        axis([mmX(1) mmX(end) mmX(1) mmX(end)])
        axis equal
        z_max_h = max(max(abs(squeeze(Proj_wf_R(:,:,ind)))));
        plot3([0 moire_L_x1(1)],[0 moire_L_x1(2)],z_max_h+[0 0],'-w')
        plot3([0 moire_L_x2(1)],[0 moire_L_x2(2)],z_max_h+[0 0],'-w')

        view([0,0,1])
    end

    hold off;
    %}


    save(wf_filename,'mmX','mmY','sc_grid1','sc_grid2','all_wfAL1','all_wfAL2','all_wfBL1','all_wfBL2','moire_L_x1','moire_L_x2');
    %save(h_filename,'all_new_hmat','all_ext_state_wf8','all_kpts','moire_L_x1','moire_L_x2','moire_k_vec1','moire_k_vec2','knum','knum_tot')
    save(h_filename,'all_new_hmat','all_new_basis','all_kpts','hex_all','moire_L_x1','moire_L_x2','moire_k_vec1','moire_k_vec2','numk','knum_tot')

    %{
    figure()
    plot(all_ssdiag','LineWidth',2);
    ylabel('SVD eigenvalue (diag(S) from A=USV)')
    xlabel('k')
    set(gca,'FontSize',20)
    %}
    %flipud(all_ssdiag_trio_combined)
    %flipud(all_ssdiag_eight)



end
