function [perturb_filename] = eightWAN_perturb_calc(h_in_filename,tar_theta,tar_folder,numk,disp_str, sub_str_sym, sub_str_asym)
    % the twist angle converted into radian

    %wf_filename = [tar_folder '/eightWANnew_wf_' strrep(num2str(tar_theta,'%.2f'),'.','p') '_' num2str(numk) 'k.mat'];
    %h_filename = [tar_folder '/eightWANnew_H_' strrep(num2str(tar_theta,'%.2f'),'.','p') '_' num2str(numk) 'k.mat'];
    perturb_filename = [tar_folder '/eightWAN_perturb_' strrep(num2str(tar_theta,'%.2f'),'.','p') '_' num2str(numk) 'k.mat'];

    load(h_in_filename);
    
    %tar_theta = 1.05;
    %{
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

    allk12=linspace(0,1,numk+1);
    allk12=allk12(1:numk);
    knum_tot=numk^2;

    [meshk1,meshk2]=meshgrid(allk12,allk12);
    all_kpts=zeros(knum_tot,3);
    all_kpts(:,1)=moire_k_vec1(1)*meshk1(:)+moire_k_vec2(1)*meshk2(:);
    all_kpts(:,2)=moire_k_vec1(2)*meshk1(:)+moire_k_vec2(2)*meshk2(:);
    %}

    % performs a single k-point calculation
    %{
    numk = 1;
    knum_tot = 1;
    kpt_tar = all_kpts(7,:);
    all_kpts = zeros(1,3);
    all_kpts = kpt_tar;
    %}
    
    [sweep_H_base,hex_all] = tblg_kp_calc_ext_getH('tar_theta',tar_theta,'k_list',all_kpts);
    [sweep_H_perturb,hex_all] = tblg_kp_calc_ext_getH('tar_theta',tar_theta,'k_list',all_kpts, ...
                                'displacement_strength',disp_str, ...
                                'sublattice_strength_sym',sub_str_sym, ...
                                'sublattice_strength_asym',sub_str_asym);
                            
    %tot_dim = 2*size(hex_all,1);
    
    h_eff_perturb = zeros(8,8,knum_tot);
    all_new_hmat = zeros(8,8,knum_tot);
    
    for indk=1:knum_tot
        %fprintf("%d / %d k-points \n",indk,knum_tot);
        %know=all_kpts(indk,1:2);
        
        rotate_state_wf = all_new_basis(:,:,indk);
        Hmat_diff = sweep_H_perturb{indk} - sweep_H_base{indk};
        
        h_eff_perturb(:,:,indk) = rotate_state_wf'*Hmat_diff*rotate_state_wf;
        
        all_new_hmat(:,:,indk) = rotate_state_wf'*sweep_H_perturb{indk}*rotate_state_wf;

    end
    

    %save(wf_filename,'mmX','mmY','sc_grid1','sc_grid2','all_wfAL1','all_wfAL2','all_wfBL1','all_wfBL2','moire_L_x1','moire_L_x2');
    %save(h_filename,'all_new_hmat','all_ext_state_wf8','all_kpts','moire_L_x1','moire_L_x2','moire_k_vec1','moire_k_vec2','knum','knum_tot')
    %save(h_filename,'all_new_hmat','all_new_basis','all_kpts','moire_L_x1','moire_L_x2','moire_k_vec1','moire_k_vec2','knum','knum_tot')
    save(perturb_filename,'all_new_hmat','h_eff_perturb','all_kpts','moire_L_x1','moire_L_x2','moire_k_vec1','moire_k_vec2','numk','knum_tot');


end






