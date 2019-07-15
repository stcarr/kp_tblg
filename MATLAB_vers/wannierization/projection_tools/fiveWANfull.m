function [wf_filename, h_filename] = fiveWANfull(trioWANcs_filename,fiveWANcs_filename,tar_theta,tar_folder,numk,numxy, side)
    % the twist angle converted into radian
    % this MacDonald effective theory has magic angle ~1.1 degree
    
    if side == 0 % top 3 bands
        side_text  = 'top';
        %sym_sign   = -1; % p_z type
        %start_band = 2; % 2 above nb/2
    elseif side == 1 % lower 3 bands
        side_text  = 'bot';
        %sym_sign   = 1; % s type
        %start_band = -3; % 3 below nb/2
    end
    
    wf_filename = [tar_folder '/fiveWANfull' side_text '_wf_' strrep(num2str(tar_theta,'%.2f'),'.','p') '_' num2str(numk) 'k.mat'];
    h_filename = [tar_folder '/fiveWANfull' side_text '_H_' strrep(num2str(tar_theta,'%.2f'),'.','p') '_' num2str(numk) 'k.mat'];


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

    knum=numk; % previously knum = 15, now passed to function

    allk12=linspace(0,1,knum+1);
    allk12=allk12(1:knum);
    knum_tot=knum^2;

    [meshk1,meshk2]=meshgrid(allk12,allk12);
    all_kpts=zeros(knum_tot,3);
    all_kpts(:,1)=moire_k_vec1(1)*meshk1(:)+moire_k_vec2(1)*meshk2(:);
    all_kpts(:,2)=moire_k_vec1(2)*meshk1(:)+moire_k_vec2(2)*meshk2(:);

    [sweep_H,hex_all] = tblg_kp_calc_ext_getH('tar_theta',tar_theta,'k_list',all_kpts);
    tot_dim = 2*size(hex_all,1);


    all_index=reshape(1:tot_dim,2,tot_dim/2);

    indL1=1:(tot_dim/2);
    indL2=(tot_dim/2+1):tot_dim;

    indL1A=indL1(1:2:(tot_dim/2));
    indL1B=indL1(2:2:(tot_dim/2));

    indL2A=indL2(1:2:(tot_dim/2));
    indL2B=indL2(2:2:(tot_dim/2));

    all_new_basis=zeros(tot_dim,5,knum_tot);
    all_new_hmat=zeros(5,5,knum_tot);


    all_index=reshape(1:tot_dim,2,tot_dim/2);


    allbands=zeros(tot_dim,knum_tot);

    allbands_eff_trio=zeros(3,knum_tot);
    allbands_eff_ext=zeros(5,knum_tot);

    fff=tot_dim/2;
    %sel_bandind=(fff-4):(fff+5);
    % sel_bandind=(fff-6):(fff+7);



    %proj_band_ind=(fff-8):(fff+9);
    %proj_band_ind=(fff-13):(fff+14);

    if (side == 0) % top 3 bands
        EEcut_low = -1;
        EEcut_high = 0.12;
    elseif (side == 1) % bottom 3 bands
        EEcut_low = -0.12;
        EEcut_high = 1;        
    end


    %save('WAN_Total_fiveband_trio_data','mmX','mmY','all_wfAL1','all_wfAL2','all_wfBL1','all_wfBL2','moire_L_x1','moire_L_x2');
    load(trioWANcs_filename);
    mmX_trio=mmX;
    mmY_trio=mmY;

    all_wfAL1_trio=all_wfAL1;
    all_wfBL1_trio=all_wfBL1;
    all_wfAL2_trio=all_wfAL2;
    all_wfBL2_trio=all_wfBL2;


    % save('Ten_trial_WF','mmX','mmY','all_wfAL1','all_wfAL2','all_wfBL1','all_wfBL2')

    %load('WAN_Total_fiveband_data_51xy');
    load(fiveWANcs_filename);
    mmX_five=mmX;
    mmY_five=mmY;

    all_wfAL1_five=all_wfAL1;
    all_wfBL1_five=all_wfBL1;
    all_wfAL2_five=all_wfAL2;
    all_wfBL2_five=all_wfBL2;



    %all_vecs=zeros(tot_dim,length(sel_bandind),knum_tot);

    all_ssdiag_trio=zeros(3,knum_tot);
    all_ssdiag_five=zeros(5,knum_tot);

    testfac=0;

    fffmid=tot_dim/2;
    if (side == 0) % top 3 bands
        all_fffind=(1:(fffmid+1));
    elseif (side == 1) % bottom 3 bands
        all_fffind=(fffmid:tot_dim);        
    end

    allbands_eff_ext5=zeros(5,knum_tot);


    ext_state_wf=zeros(tot_dim,5);
    all_ext_state_wf5=zeros(tot_dim,5,knum_tot);
    all_ext_state_wf=zeros(tot_dim,5,knum_tot);

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


        %all_vecs(:,:,indk)=V(:,sel_bandind);


        %PPk_proj=V(:,sel_bandind)*(V(:,sel_bandind)');

        PPk_bande=allbands(:,indk);

        %PPk_proj=zeros(tot_dim);
        band_weight=zeros(tot_dim,1);

        tmp_ind=(PPk_bande>=(EEcut_low))&(PPk_bande<=(EEcut_high));
        band_weight(tmp_ind)=1;

        if (side == 0)
            tmp_ind=(PPk_bande<(EEcut_low));
            %tmpe=(-EEcut1-PPk_bande(tmp_ind));
            band_weight(tmp_ind)=0;

            tmp_ind=(PPk_bande>(EEcut_high));
            tmpe=(PPk_bande(tmp_ind)-EEcut_high);
            band_weight(tmp_ind)=exp(-4*(tmpe/0.2).^2);
        elseif (side == 1)
            tmp_ind=(PPk_bande<(EEcut_low));
            tmpe=(PPk_bande(tmp_ind)-EEcut_low);
            band_weight(tmp_ind)=exp(-4*(tmpe/0.2).^2);
            
            tmp_ind=(PPk_bande>(EEcut_high));
            band_weight(tmp_ind)=0;
        end

        band_weight(all_fffind)=0;

        PPk_proj=V*(diag(band_weight))*V';

        % Trial WF part is slow
        % trial WF

        kktmpL1=shift_klist(1:(tot_dim/4),:);
        kktmpL2=shift_klist((tot_dim/4+1):(tot_dim/2),:);

        % 3 band


        trial_wf_bloch_trio=zeros(tot_dim,3);

        for indw=1:3

            wfAL1=squeeze(all_wfAL1_trio(:,:,indw));
            wfBL1=squeeze(all_wfBL1_trio(:,:,indw));
            wfAL2=squeeze(all_wfAL2_trio(:,:,indw));
            wfBL2=squeeze(all_wfBL2_trio(:,:,indw));


            for indt=1:(tot_dim/4)
                kkkL1=kktmpL1(indt,:);
                kkkL2=kktmpL2(indt,:);

                tmp1=exp(-1i*(kkkL1(1)*mmX_trio+kkkL1(2)*mmY_trio));
                tmp2=exp(-1i*(kkkL2(1)*mmX_trio+kkkL2(2)*mmY_trio));

                tmpo1=sum(wfAL1(:).*tmp1(:));
                tmpo2=sum(wfBL1(:).*tmp1(:));
                tmpo3=sum(wfAL2(:).*tmp2(:));
                tmpo4=sum(wfBL2(:).*tmp2(:));

                trial_wf_bloch_trio(indL1A(indt),indw)=tmpo1;
                trial_wf_bloch_trio(indL1B(indt),indw)=tmpo2;
                trial_wf_bloch_trio(indL2A(indt),indw)=tmpo3;
                trial_wf_bloch_trio(indL2B(indt),indw)=tmpo4;


            end
            trial_wf_bloch_trio(:,indw)=trial_wf_bloch_trio(:,indw)/sqrt(trial_wf_bloch_trio(:,indw)'*trial_wf_bloch_trio(:,indw));

        end

        trial_wf_proj_trio=PPk_proj*trial_wf_bloch_trio;
        % end trial WF part

        %diag(trial_wf_proj'*trial_wf_proj);

    %     for indw=1:10
    %         tmpwf=trial_wf_proj(:,indw);
    %         trial_wf_proj(:,indw)=trial_wf_proj(:,indw)/sqrt(abs(tmpwf'*tmpwf));
    %         
    %     end
        Band_state=V(:,:);

        Amat=Band_state'*trial_wf_proj_trio;


        [Umat,Smat,Vmat] = svd(Amat);

        ssdiag=(diag(Smat'*Smat));
        all_ssdiag_trio(:,indk)=sort(diag(Smat));
        sort(diag(Smat));
        Sinv_root=Vmat*diag(1./sqrt(ssdiag))*Vmat';

        %Amat-Umat*Smat*Vmat';
        %Amat'*Amat-Vmat*Smat*Smat*Vmat';

        %Sinv_root*Vmat*Smat*Smat*Vmat';

        %Amat'*Amat*Sinv_root*Sinv_root';

        rotate_state_wf=trial_wf_proj_trio*Sinv_root;

        %mesh(abs(ppp'*ppp))


        Heff=rotate_state_wf'*Hmat*rotate_state_wf;


        allbands_eff_trio(:,indk)=sort(real(eig(Heff)));

        ext_state_wf(:,1:3)=rotate_state_wf;

        ext_state_wf(:,4:5)=V(:,(tot_dim/2):(tot_dim/2+1));

        all_ext_state_wf(:,:,indk)=ext_state_wf(:,:);


        Heff_ext=ext_state_wf'*Hmat*ext_state_wf;

        allbands_eff_ext(:,indk)=sort(real(eig(Heff_ext)));


        % five band

        trial_wf_bloch_five=zeros(tot_dim,5);

        for indw=1:5

            wfAL1=squeeze(all_wfAL1_five(:,:,indw));
            wfBL1=squeeze(all_wfBL1_five(:,:,indw));
            wfAL2=squeeze(all_wfAL2_five(:,:,indw));
            wfBL2=squeeze(all_wfBL2_five(:,:,indw));


            for indt=1:(tot_dim/4)
                kkkL1=kktmpL1(indt,:);
                kkkL2=kktmpL2(indt,:);

                tmp1=exp(-1i*(kkkL1(1)*mmX_five+kkkL1(2)*mmY_five));
                tmp2=exp(-1i*(kkkL2(1)*mmX_five+kkkL2(2)*mmY_five));

                tmpo1=sum(wfAL1(:).*tmp1(:));
                tmpo2=sum(wfBL1(:).*tmp1(:));
                tmpo3=sum(wfAL2(:).*tmp2(:));
                tmpo4=sum(wfBL2(:).*tmp2(:));

                trial_wf_bloch_five(indL1A(indt),indw)=tmpo1;
                trial_wf_bloch_five(indL1B(indt),indw)=tmpo2;
                trial_wf_bloch_five(indL2A(indt),indw)=tmpo3;
                trial_wf_bloch_five(indL2B(indt),indw)=tmpo4;


            end
            trial_wf_bloch_five(:,indw)=trial_wf_bloch_five(:,indw)/sqrt(trial_wf_bloch_five(:,indw)'*trial_wf_bloch_five(:,indw));

        end

        trial_wf_bloch_five_proj=ext_state_wf*ext_state_wf'*trial_wf_bloch_five;

        Amat=ext_state_wf'*trial_wf_bloch_five_proj;


        [Umat,Smat,Vmat] = svd(Amat);

        ssdiag=(diag(Smat'*Smat));
        all_ssdiag_five(:,indk)=sort(diag(Smat));
        sort(diag(Smat));
        Sinv_root=Vmat*diag(1./sqrt(ssdiag))*Vmat';

        rotate_state_wf_ext=trial_wf_bloch_five_proj*Sinv_root;


        all_ext_state_wf5(:,:,indk)=rotate_state_wf_ext(:,:);


        Heff_ext5=rotate_state_wf_ext'*Hmat*rotate_state_wf_ext;

        allbands_eff_ext5(:,indk)=sort(real(eig(Heff_ext5)));

        all_new_hmat(:,:,indk)=Heff_ext5(:,:);

        %toc

    end
    
    % build wannier states

    %LLmoire=sqrt(dot(moire_L_x1,moire_L_x1));
    %numxy=51;%501;
    %allxy=2.0*LLmoire*linspace(-1,1,numxy);

    %midpt=(numxy+1)/2;
    %[mmX,mmY]=meshgrid(allxy,allxy);

    all_new_basis = all_ext_state_wf5;

    LLmoire=sqrt(dot(moire_L_x1,moire_L_x1));
    %numxy=51;%401; passed to function
    %allxy=2.0*LLmoire*linspace(-1,1,numxy);
    %[mmX,mmY]=meshgrid(allxy,allxy);
    
    %allxy = 3*linspace(-1,1,3*numxy);
     allxy = 3*linspace(-1,1,6*numxy+1);

    [sc_grid1,sc_grid2]=meshgrid(allxy,allxy);
    
    mmX = sc_grid1*moire_L_x1(1) + sc_grid2*moire_L_x2(1);
    mmY = sc_grid1*moire_L_x1(2) + sc_grid2*moire_L_x2(2);
    mmm=size(mmX,1);

    mmm=size(mmX,1);
    all_wfAL1=zeros(mmm,mmm,5);
    all_wfBL1=zeros(mmm,mmm,5);
    all_wfAL2=zeros(mmm,mmm,5);
    all_wfBL2=zeros(mmm,mmm,5);


    for wannier_ind=1:5
        fprintf('done with %d/%d wannier constructions \n',wannier_ind,5);
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


    for ind=1:5
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
    Proj_wf_R(:,:,4)=all_wfBL2(:,:,wannier_ind);

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
    %save(h_filename,'all_new_hmat','all_ext_state_wf5','all_kpts','moire_L_x1','moire_L_x2','moire_k_vec1','moire_k_vec2','knum','knum_tot')
    save(h_filename,'all_new_hmat','all_kpts','moire_L_x1','moire_L_x2','moire_k_vec1','moire_k_vec2','knum','knum_tot')
    
end


    % 
    % figure(2)
    % subplot(1,2,2)
    % 
    % plot(scale_axis,all_ssdiag_trio(:,:)','LineWidth',2)
    % axis([-inf,inf,0,inf])
    % xlabel('k')
    % ylabel('SVD eigenvalue (diag(S) from A=USV)')
    % set(gca,'FontSize',20)
    % 
    % 
    % subplot(1,2,1)
    % 
    % plot(scale_axis,allbands','r','LineWidth',3);
    % hold on;
    % 
    % %plot(scale_axis,allbands(sel_bandind,:),'b','LineWidth',2);
    % 
    % plot(scale_axis,allbands_eff_trio','k','LineWidth',2);
    % hold off;
    % axis([-inf,inf,-0.35,0.35])
    % 
    % 
    % xlabel('k')
    % ylabel('Energy (eV)')
    % set(gca,'FontSize',20)
    % 






