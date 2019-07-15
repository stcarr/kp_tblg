function [wf_filename] = fiveWAN_WF_modifyU(h_in_filename,wf_in_filename,modify_U_filename,tar_theta,tar_folder,side,numk)
    % the twist angle converted into radian
    % this MacDonald effective theory has magic angle ~1.1 degree
    
    if (side == 0)
        side_str = 'top';
    elseif (side == 1)
        side_str = 'bot';
    end
    num_orbs = 5;
    
    wf_filename = [tar_folder '/fiveWAN' side_str '_modifyU_wf_' strrep(num2str(tar_theta,'%.2f'),'.','p') '_' num2str(numk) 'k.mat'];
    %h_filename = [tar_folder '/eightWANnew_H_' strrep(num2str(tar_theta,'%.2f'),'.','p') '_' num2str(numk) 'k.mat'];
    
    hex_all = 0;
    all_new_basis = 0;
    U_store = 0;
    
    load(h_in_filename,'hex_all','all_new_basis','all_kpts','moire_L_x1','moire_L_x2');
    load(wf_in_filename,'sc_grid1','sc_grid2');
    load(modify_U_filename,'U_store');
    
    % Add constant phase to all WF after REWAN, gives clearer sym. properties
    % Should be same phase as in real_space_sym_5band.m "Rephase WFs" section
    if (side == 0)
        WF_phase = [exp(1j*pi/4),exp(1j*3*pi/4),1,1j,1j];
    elseif (side == 1)
        WF_phase = [exp(1j*pi/4),exp(1j*3*pi/4),1,-1j,1j];        
    end    
    %{
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
    %}
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
    
    %[~,hex_all] = tblg_kp_calc_ext_getH('tar_theta',tar_theta,'k_list',all_kpts);
    tot_dim = 2*size(hex_all,1);

    %all_index=reshape(1:tot_dim,2,tot_dim/2);

    indL1=1:(tot_dim/2);
    indL2=(tot_dim/2+1):tot_dim;

    indL1A=indL1(1:2:(tot_dim/2));
    indL1B=indL1(2:2:(tot_dim/2));

    indL2A=indL2(1:2:(tot_dim/2));
    indL2B=indL2(2:2:(tot_dim/2));
    
    % redefine mmX and mmY for the current theta
    % (in case we are loading from a nearby result).
    %allxy = 3*linspace(-1,1,3*numxy);
    %allxy = 3*linspace(-1,1,6*numxy+1);

    %[sc_grid1,sc_grid2]=meshgrid(allxy,allxy);
    
    knum_tot = size(all_kpts,1);
    
    mmX = sc_grid1*moire_L_x1(1) + sc_grid2*moire_L_x2(1);
    mmY = sc_grid1*moire_L_x1(2) + sc_grid2*moire_L_x2(2);
    mmm=size(mmX,1);
    
    % build wannier states
    modify_U_basis = zeros(size(all_new_basis));
    for indk = 1:knum_tot
        modify_U_basis(:,:,indk) = all_new_basis(:,:,indk)*U_store(:,:,indk);
    end
    
    all_wfAL1=zeros(mmm,mmm,num_orbs);
    all_wfBL1=zeros(mmm,mmm,num_orbs);
    all_wfAL2=zeros(mmm,mmm,num_orbs);
    all_wfBL2=zeros(mmm,mmm,num_orbs);


    for wannier_ind=1:num_orbs
        fprintf('starting %d/%d wannier constructions \n',wannier_ind,num_orbs);
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

            wfnow=squeeze(modify_U_basis(:,wannier_ind,indk));

            for indt=1:(tot_dim/4)
                kkkL1=kktmpL1(indt,:);
                kkkL2=kktmpL2(indt,:);

                tmp1=exp(1j*(kkkL1(1)*mmX+kkkL1(2)*mmY));
                tmp2=exp(1j*(kkkL2(1)*mmX+kkkL2(2)*mmY));
                
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


    % normalize
    for ind=1:num_orbs
        tmp_AL1=squeeze(all_wfAL1(:,:,ind));
        tmp_AL2=squeeze(all_wfAL2(:,:,ind));
        tmp_BL1=squeeze(all_wfBL1(:,:,ind));
        tmp_BL2=squeeze(all_wfBL2(:,:,ind));

        tmp_sum=sqrt(sum(abs(tmp_AL1(:)).^2)+sum(abs(tmp_BL1(:)).^2)+sum(abs(tmp_AL2(:)).^2)+sum(abs(tmp_BL2(:)).^2));

        ph = WF_phase(ind);
        all_wfAL1(:,:,ind)=ph*tmp_AL1/tmp_sum;
        all_wfAL2(:,:,ind)=ph*tmp_AL2/tmp_sum;
        all_wfBL1(:,:,ind)=ph*tmp_BL1/tmp_sum;
        all_wfBL2(:,:,ind)=ph*tmp_BL2/tmp_sum;

    end

    save(wf_filename,'mmX','mmY','sc_grid1','sc_grid2','all_wfAL1','all_wfAL2','all_wfBL1','all_wfBL2','moire_L_x1','moire_L_x2');
    %save(h_filename,'all_new_hmat','all_ext_state_wf8','all_kpts','moire_L_x1','moire_L_x2','moire_k_vec1','moire_k_vec2','knum','knum_tot')
    %save(h_filename,'all_new_hmat','all_new_basis','all_kpts','moire_L_x1','moire_L_x2','moire_k_vec1','moire_k_vec2','knum','knum_tot')

end
