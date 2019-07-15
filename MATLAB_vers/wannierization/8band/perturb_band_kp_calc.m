function [kp_bands, scale_axis] = perturb_band_kp_calc(tar_theta,numk,disp_str, sub_str_sym, sub_str_asym)
    % the twist angle converted into radian

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

    kk1=0*moire_k_vec1;
    kk2=-moire_k_vec1/2;
    kk3=-(moire_k_vec1-moire_k_vec2)/3;


    %kscan_list=[kk1;kk2;kk3;kk1];
    kscan_list=[kk3;kk1;kk2;kk3];
    scan_klist(:,3)=0;


    [ all_kpts, scale_axis] = generate_k_line( knum, kscan_list );
    knum_tot = size(all_kpts,1);
    
    %{
    allk12=linspace(0,1,knum+1);
    allk12=allk12(1:knum);
    knum_tot=knum^2;

    [meshk1,meshk2]=meshgrid(allk12,allk12);
    all_kpts=zeros(knum_tot,3);
    all_kpts(:,1)=moire_k_vec1(1)*meshk1(:)+moire_k_vec2(1)*meshk2(:);
    all_kpts(:,2)=moire_k_vec1(2)*meshk1(:)+moire_k_vec2(2)*meshk2(:);
    %}
    
    % performs a single k-point calculation
    %{
    knum = 1;
    knum_tot = 1;
    kpt_tar = all_kpts(7,:);
    all_kpts = zeros(1,3);
    all_kpts = kpt_tar;
    %}
    
    %[sweep_H_base,hex_all] = tblg_kp_calc_ext_getH('tar_theta',tar_theta,'k_list',all_kpts);
    [sweep_H_perturb,hex_all] = tblg_kp_calc_ext_getH('tar_theta',tar_theta,'k_list',all_kpts, ...
                                'displacement_strength',disp_str, ...
                                'sublattice_strength_sym',sub_str_sym, ...
                                'sublattice_strength_asym',sub_str_asym);
                            
    tot_dim = 2*size(hex_all,1);
    kp_bands = zeros(tot_dim, knum_tot);
    
    for indk=1:knum_tot
        %fprintf("%d / %d k-points \n",indk,knum_tot);
        %know=all_kpts(indk,1:2);
        
        kp_bands(:,indk) = sort(real(eig(sweep_H_perturb{indk})),'ascend');
    end
    

end






