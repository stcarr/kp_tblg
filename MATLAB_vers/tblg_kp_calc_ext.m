function [sweep_vals, scaleaxis, sweep_kpts] = tblg_kp_calc_ext(varargin)

    % Computes eigenvalues for selected k-points according to the Name/Value
    % pair list supplied in varargin.
    % Here we let number of bands = nb, and the number of k points = kp.
    %
    %   sweep_vals: Cell containing the nb x nk eigenvalues for each theta
    %   scaleaxis:  Cell containing the proper scaling for high-sym linecut
    %   sweep_kpts: Cell containing the kpts sampled for each theta

    % Some useful k-point info:
    %   our moire BZ gamma point is at [0 0]
    %   reciprocal vectors are given by hex_b1 and hex_b2 later in code
    %   kk1a/b/c give the K points, kk4 gives M point (depend on hex_b's)
    %   search for kk1a or kk4 in code to find their definition
    
    % option structure copied from Robert Cain on stackoverflow
    opts = struct(  ... % kp-model assumptions
                        'relax_type','full_relax', ... %type of relaxation
                        ... % ^ 'full_relax', 'flat_relax', or 'no_relax' [not implemented yet]
                        'use_bmd',0,  ... % Use BMD 2011 model [not im[plemented yet]
                    ... % k-sampling options
                        'vf_only',0,  ... % only computes fermi velocity
                        'gam_only',0, ... % only compute at gamma point
                        'full_bz',0, ...  % does a full BZ sampling
                    ... % calculation settings
                        'theta_list',[1.00], ... % list of desired thetas
                        'knum',41, ... % how many k-point sampling points
                        'ax_m',0.1, ... % Energy max in plotting tool
                        'plot_on',0, ... % set to 1 to turn on band plots
                    ... % controls kp terms in construction of H
                        ...%'inter_shells',3, ... % keep up to 3rd NN inter couplings
                        'inter_qdep_shells',1, ... % keep NN inter k+- terms
                        'inter_fac',1.0, ... % controls strength of inter couplings
                        'strain_fac',1.0 ... % controls strength of in-plane strain gauge field
                    );

    
    % read the acceptable names
    optionNames = fieldnames(opts);

    % count arguments
    nArgs = length(varargin);
    if round(nArgs/2)~=nArgs/2
       error('tblg_kp_calc takes only propertyName/propertyValue pairs')
    end

    for pair = reshape(varargin,2,[]) %# pair is {propName;propValue}
       inpName = lower(pair{1}); %# make case insensitive

       if any(strcmp(inpName,optionNames))
          % overwrite options. If you want you can test for the right class here
          % Also, if you find out that there is an option you keep getting wrong,
          % you can use "if strcmp(inpName,'problemOption'),testMore,end"-statements
          opts.(inpName) = pair{2};
       else
          error('%s is not a recognized parameter name',inpName)
       end
    end
    
    tar_theta_list = opts.theta_list;

    for tar_theta_idx = 1:length(tar_theta_list)

        clearvars -except tar_theta_idx tar_theta_list opts ...
                    sweep_vals scaleaxis sweep_kpts
        
        if (strcmp(opts.relax_type, 'full_relax'))
            %load('dft_full_relax_data_02-04-2019.mat');
            filename = '../data/full_relax_kp_01-06-2019.dat';
            [thetas, inter_kp, intra_bot_kp, intra_top_kp, inter_shells, intra_shells] = parse_datafile(filename);
        elseif (strcmp(opts.relax_type, 'flat_relax'))
            %load('dft_flat_relax_kp_data_11-02-2018.mat')
            fprintf('FLAT RELAX NOT IMPLEMENTED YET \n');
        else % load no_relax as default otherwise
            %load('dft_no_relax_kp_data_11-02-2018.mat');
            fprintf('NO RELAX NOT IMPLEMENTED YET \n');
        end


        % set a cut-off condition for the k-p model
        % that scales automatically with theta
        hex_cut_fac0 = 5;
        theta_cut_0 = 1.5;

        % next two terms control the weight of strain and interlayer coupling
        %strain_fac = 1.00;   % in-plane strain correction scaling
        %inter_fac = 1.00; % inter-layer coupling correction scaling

        % terms which control number of k-p coupling terms to include
        %{
        if opts.inter_shells == 1
            max_inter_q = 3;
        elseif opts.inter_shells == 2
            max_inter_q = 6;
        else
            max_inter_q = 12; % [#: includes] 3: nearest, 6: next nearest, 12: 3rd nearest neighbors
        end
        %}

        if opts.inter_qdep_shells == 0
            max_inter_q_plusminus = 0;
        else
            max_inter_q_plusminus = 3; % only do k+- correction for 3 nearest interlayer q's
        end

        max_intra_q_plusminus = 0; % no k+- correction for intralayer

        
        % get the angle and then interpolate the kp-model parameters
        tar_theta = tar_theta_list(tar_theta_idx);

        thetas_deg = thetas*180/pi;

        hex_cut_fac = hex_cut_fac0*sqrt(theta_cut_0/tar_theta);
        hex_cut_fac = max(hex_cut_fac,hex_cut_fac0);

        All_Eff_inter_ext = squeeze(interp1(thetas_deg,inter_kp,tar_theta));
        All_Eff_intra_bot_ext = squeeze(interp1(thetas_deg,intra_bot_kp,tar_theta));
        All_Eff_intra_top_ext = squeeze(interp1(thetas_deg,intra_top_kp,tar_theta));
        rot_theta = tar_theta*pi/180;

        
        % clean up NAN errors during fitting (small terms -> singular)
        All_Eff_inter=squeeze(All_Eff_inter_ext(3,:,:,:));
        All_Eff_intra_bot=squeeze(All_Eff_intra_bot_ext(3,:,:,:))*1;
        All_Eff_intra_top=squeeze(All_Eff_intra_top_ext(3,:,:,:))*1;
        All_Eff_inter_kplus=squeeze(All_Eff_inter_ext(1,:,:,:));
        All_Eff_inter_kminus=squeeze(All_Eff_inter_ext(2,:,:,:));
        
        All_Eff_inter(isnan(All_Eff_inter)) = 0;
        All_Eff_intra_bot(isnan(All_Eff_intra_bot)) = 0;
        All_Eff_intra_top(isnan(All_Eff_intra_top)) = 0;
        All_Eff_inter_kplus(isnan(All_Eff_inter_kplus)) = 0;
        All_Eff_inter_kminus(isnan(All_Eff_inter_kminus)) = 0;
        
        All_Eff_intra_shell_indices = intra_shells;
        All_Eff_inter_shell_indices = inter_shells;
        
        % Symm enforcer, to make sure the numerically calculated terms
        % do not break any assumed symmetery of our k-p Hamiltonian
        [All_Eff_inter,All_Eff_inter_kplus,All_Eff_inter_kminus] = TwBLG_KP_sym_enforce_inter(All_Eff_inter,All_Eff_inter_kplus,All_Eff_inter_kminus);

        % we turn off k-dependent terms for intralayer coupling for now
        kdep_fac2=0.0;

        All_Eff_intra_bot_kplus=squeeze(All_Eff_intra_bot_ext(1,:,:,:))*kdep_fac2;
        All_Eff_intra_top_kplus=squeeze(All_Eff_intra_top_ext(1,:,:,:))*kdep_fac2;

        All_Eff_intra_bot_kminus=squeeze(All_Eff_intra_bot_ext(2,:,:,:))*kdep_fac2;
        All_Eff_intra_top_kminus=squeeze(All_Eff_intra_top_ext(2,:,:,:))*kdep_fac2;

        % this uses the Bistritzer/MacDonald (BMD) model from 2011 PNAS paper
        if (opts.use_bmd == 1)
            bmd_w=0.11;
            All_Eff_inter=zeros(2,2,12);

            All_Eff_inter(:,:,1)=bmd_w*[1,1;1,1];
            All_Eff_inter(:,:,2)=bmd_w*[1,exp(-i*2*pi/3);exp(i*2*pi/3),1];
            All_Eff_inter(:,:,3)=bmd_w*[1,exp(i*2*pi/3);exp(-i*2*pi/3),1];

            All_Eff_intra_bot=zeros(2,2,12);
            All_Eff_intra_top=zeros(2,2,12);
            opts.strain_fac = 0;

            All_Eff_inter_kplus=zeros(2,2,12);
            All_Eff_intra_bot_kplus=zeros(2,2,12);
            All_Eff_intra_top_kplus=zeros(2,2,12);

            All_Eff_inter_kminus=zeros(2,2,12);
            All_Eff_intra_bot_kminus=zeros(2,2,12);
            All_Eff_intra_top_kminus=zeros(2,2,12);    
        end


        %[All_Eff_inter,All_Eff_intra_bot,All_Eff_intra_top] = Fit_Model_parms(fname_parm,rot_theta);

        lattice_a=1.42*sqrt(3);

        unit_dim=2;


        % sets the linear and quadratic term of the Dirac cones
        Dirac_v1=5.2268;
        Dirac_v2=2.2450*0; % no quadratic for now
        % Dirac_v=5.2;
        Dirac_diag=-1.1142*1;

        % some important length scales for the Brill. Zone
        KD=4*pi/3/lattice_a;
        KTH=2*KD*sin(rot_theta/2);

        HEX_BLEN=KTH*sqrt(3);

        % basic pauli matrices
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



        % create hex table of possible intralayer interactions
        %num_intra_qs=12;
        
        
        % keep all inter terms for now
        max_inter_q = size(All_Eff_inter_shell_indices,2);%size(inter_kp,5);

        max_intra_q = size(All_Eff_intra_shell_indices,2);

        if tar_theta > 1.2
            max_intra_q = 3;
            max_inter_q = 3;
        elseif tar_theta > 0.6
            max_intra_q = 5;
            max_inter_q = 5;
        elseif tar_theta > 0.3
            max_intra_q = 12;
            max_inter_q = 12;            
        end
        
        %max_intra_q = 3;
        %max_inter_q = 3;
        
        %num_intra_qs = max_intra_q;%size(All_Eff_intra_bot,3);
        
        intra_qs_idx = 1;
        for shell_idx = 1:max_intra_q
            shell_here = All_Eff_intra_shell_indices{shell_idx};
            n_momenta = size(shell_here,1);
            if n_momenta > 0
                intra_qs(intra_qs_idx:intra_qs_idx+n_momenta-1,:) = shell_here(:,1:2);
            end
            intra_qs_idx = intra_qs_idx + n_momenta;
            
        end
        num_intra_qs = intra_qs_idx-1;
        
        % create hex table of possible interlayer interactions
        %num_inter_qs = max_inter_q;%size(All_Eff_inter,3);
        inter_qs_idx = 1;
        for shell_idx = 1:max_inter_q
            shell_here = All_Eff_inter_shell_indices{shell_idx};
            n_momenta = size(shell_here,1);
            if n_momenta > 0
                inter_qs(inter_qs_idx:inter_qs_idx+n_momenta-1,:) = shell_here;
            end
            inter_qs_idx = inter_qs_idx + n_momenta;
            
        end
        num_inter_qs = inter_qs_idx - 1;
        
        %{
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
        %}

        miller_b1 = hex_b1;
        miller_b2 = hex_b1 + hex_b2;
        intra_qs_vec(1:num_intra_qs,1)=intra_qs(:,1)*miller_b1(1)+intra_qs(:,2)*miller_b2(1);
        intra_qs_vec(1:num_intra_qs,2)=intra_qs(:,1)*miller_b1(2)+intra_qs(:,2)*miller_b2(2);
        %norm(intra_qs_vec(end,:))
        %{
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
        %}
        
          % interq miller construct
        a_rot = 120;

        r_mat = [cosd(a_rot) sind(a_rot); -sind(a_rot) cosd(a_rot)];

        a1 = (-2*hex_shift+hex_b2)'; % first INTERQ
        a2 = r_mat*a1;
        a3 = r_mat*a2;

        interall_given_qs=inter_qs(:,1)*a1'+inter_qs(:,2)*a2' + inter_qs(:,3)*a3';
        %norm(interall_given_qs(end,:))
        L12_qvecs=interall_given_qs(:,1:2);

        % loop over a large grid of hex points to find all valid couplings
        hex_M=100;

        hex_table=zeros(2*hex_M+1);

        hex_index=0;
        hex_coor=0;
        hex_cut=hex_cut_fac*HEX_BLEN;%3.21*HEX_BLEN;

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


        % This gives length in nm, should be ~14 nm for magic-angle
        %moire_L_x2*2.46/10;

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


                dvec_L12 = pvec2_L2 - pvec1_L1 - L12_qvecs(:,1:2);
                dvec_L12_match = ( sqrt(dvec_L12(:,1).^2 + dvec_L12(:,2).^2) < vecthres );
                [match_val, match_idx] = max(dvec_L12_match);
                if (match_val == 1)
                    connect_Mat_L12(ind2,ind1,match_idx) = 1;
                end
                
                %{
                for indty=1:num_inter_qs
                    tmp_qvec=L12_qvecs(indty,1:2);
                    dvec_L12=pvec2_L2-pvec1_L1-tmp_qvec;
                    if sqrt(dot(dvec_L12,dvec_L12))<vecthres
                        connect_Mat_L12(ind2,ind1,indty)=1;
                    end
                end
                %}

                dvec_L1 = pvec2_L1 - pvec1_L1 - intra_qs_vec(:,1:2);
                dvec_L2 = pvec2_L2 - pvec1_L2 - intra_qs_vec(:,1:2);

                dvec_L1_match = ( sqrt(dvec_L1(:,1).^2 + dvec_L1(:,2).^2) < vecthres );
                dvec_L2_match = ( sqrt(dvec_L2(:,1).^2 + dvec_L2(:,2).^2) < vecthres );

                [match_val_L1, match_idx_L1] = max(dvec_L1_match);
                [match_val_L2, match_idx_L2] = max(dvec_L2_match);

                if (match_val_L1 == 1)
                    connect_Mat_L1(ind2,ind1,match_idx_L1) = 1;
                end
                if (match_val_L2 == 1)
                    connect_Mat_L2(ind2,ind1,match_idx_L2) = 1;
                end                
                %{
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
                %}


            end
        end


        %

        % expfac=exp(-1*i*rot_theta);
        % 
        %  Layer1_ham=@(kknow) [0,Dirac_v*(-i*kknow(1)+kknow(2));Dirac_v*(i*kknow(1)+kknow(2)),0];
        %  Layer2_ham=@(kknow) [0,expfac*Dirac_v*(-i*kknow(1)+kknow(2));conj(expfac)*Dirac_v*(i*kknow(1)+kknow(2)),0];
        % 

        % this term is not needed, can break the symmetry at K-point
        expfac1=exp(0*i*rot_theta/2);
        expfac2=exp(-0*i*rot_theta/2);

        Layer1_ham=@(kknow) [Dirac_diag*dot(kknow,kknow),expfac1*Dirac_v1*(-i*kknow(1)+kknow(2))+(expfac2^2)*Dirac_v2*(kknow(1)-i*kknow(2))^2;conj(expfac1)*Dirac_v1*(i*kknow(1)+kknow(2))+(expfac1^2)*Dirac_v2*(kknow(1)+i*kknow(2))^2,Dirac_diag*dot(kknow,kknow)];
        Layer2_ham=@(kknow) [Dirac_diag*dot(kknow,kknow),expfac2*Dirac_v1*(-i*kknow(1)+kknow(2))+(expfac2^2)*Dirac_v2*(kknow(1)-i*kknow(2))^2;conj(expfac2)*Dirac_v1*(i*kknow(1)+kknow(2))+(expfac1^2)*Dirac_v2*(kknow(1)+i*kknow(2))^2,Dirac_diag*dot(kknow,kknow)];

        all_index_L1=reshape(1:(num_hex*unit_dim),unit_dim,num_hex);
        all_index_L2=reshape((num_hex*unit_dim+1):(2*num_hex*unit_dim),unit_dim,num_hex);
        all_index=reshape(1:tot_dim,unit_dim,tot_dim/unit_dim);

        % indL1=1:(tot_dim/2);
        % indL2=(tot_dim/2+1):tot_dim;

        
        % Default K point sampling, goes along:
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


        knum=opts.knum;

        [ all_kpts1, scale_axis1] = generate_k_line( knum, scan_klist1 );
        [ all_kpts2, scale_axis2] = generate_k_line( knum, scan_klist2 );
        % combine the two
        all_kpts1 = [all_kpts1; all_kpts2];

        % special k-point sampling methods
        if (opts.vf_only == 1) % only try to compute fermi velocity
            all_kpts1 = all_kpts1([1,2,size(all_kpts1,1)-1],:);
        elseif (opts.gam_only == 1) % only compute at SC gamma point
            all_kpts1 = [0 0 0];
        elseif (opts.full_bz == 1) % samples over entire Brill. zone
            
            % use these for periodicity in interpolated DoS calc
            bz_base_1 = hex_b1;
            bz_base_2 = hex_b2;

            % use these for irred. BZ
            %bz_base_1 = -kk1b;
            %bz_base_2 = -kk1c;

            bz_n = opts.knum;
            dk = 1/bz_n;
            grid_mesh = 0:dk:1;
            [x_grid,y_grid] = meshgrid(grid_mesh,grid_mesh);
            kx_grid = bz_base_1(1)*x_grid + bz_base_2(1)*y_grid;
            ky_grid = bz_base_1(2)*x_grid + bz_base_2(2)*y_grid;
            all_kpts1 = [kx_grid(:) ky_grid(:)];
            all_kpts1(:,3) = 0;
        end

        knum_tot=size(all_kpts1);
        knum_tot=knum_tot(1);

        Hmat_strain_L1=zeros(tot_dim);
        Hmat_strain_L2=zeros(tot_dim);

        Hmat_strain_L1_kplus=zeros(tot_dim);
        Hmat_strain_L2_kplus=zeros(tot_dim);
        Hmat_strain_L1_kminus=zeros(tot_dim);
        Hmat_strain_L2_kminus=zeros(tot_dim);

        for indq=1:num_intra_qs
            %indq
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
            %indq
            T_tmp=squeeze(All_Eff_inter(:,:,indq));
            Hmat_inter(all_index_L2(:),all_index_L1(:))=Hmat_inter(all_index_L2(:),all_index_L1(:))+kron(squeeze(connect_Mat_L12(:,:,indq)),T_tmp);

            if (indq <= max_inter_q_plusminus) % only keep k+- for nearest inter terms
            T_tmp_kplus=squeeze(All_Eff_inter_kplus(:,:,indq))*lattice_a;
            Hmat_inter_kplus0(:,:)=Hmat_inter_kplus0(:,:)+kron(squeeze(connect_Mat_L12(:,:,indq)),T_tmp_kplus);

            T_tmp_kminus=squeeze(All_Eff_inter_kminus(:,:,indq))*lattice_a;
            Hmat_inter_kminus0(:,:)=Hmat_inter_kminus0(:,:)+kron(squeeze(connect_Mat_L12(:,:,indq)),T_tmp_kminus);
            end

        end
        Hmat_inter=Hmat_inter+Hmat_inter';

        Hmat_inter_qplus0=zeros(tot_dim);
        Hmat_inter_qminus0=zeros(tot_dim);

        for indq=1:max_inter_q_plusminus 

            given_q=interall_given_qs(indq,:);
            given_qplus=given_q(1)+i*given_q(2);
            given_qminus=given_q(1)-i*given_q(2);

            T_tmp_kplus=squeeze(All_Eff_inter_kplus(:,:,indq))*lattice_a;
            Hmat_inter_qplus0(all_index_L2(:),all_index_L1(:))=Hmat_inter_qplus0(all_index_L2(:),all_index_L1(:))+kron(squeeze(connect_Mat_L12(:,:,indq)),T_tmp_kplus)*given_qplus/2;


            T_tmp_kminus=squeeze(All_Eff_inter_kminus(:,:,indq))*lattice_a;
            Hmat_inter_qminus0(all_index_L2(:),all_index_L1(:))=Hmat_inter_qminus0(all_index_L2(:),all_index_L1(:))+kron(squeeze(connect_Mat_L12(:,:,indq)),T_tmp_kminus)*given_qminus/2;

        end

        Hmat_inter_qdep=Hmat_inter+(Hmat_inter_qplus0+Hmat_inter_qminus0+Hmat_inter_qplus0'+Hmat_inter_qminus0')*1;

        %

        allbands1=zeros(tot_dim,knum_tot);
        allbands2=zeros(tot_dim,knum_tot);

        for indk=1:knum_tot
            % just to see how long it is taking, can be turned off
            fprintf("%d/%d k sampling \n",indk,knum_tot);
            
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

            Hmat_inter_kk(all_index_L2(:),all_index_L1(:))=Hmat_inter_kk(all_index_L2(:),all_index_L1(:))+(Hmat_inter_kplus0*kplus_L1+Hmat_inter_kminus0*kminus_L1)*1;

            Hmat_inter_kk=Hmat_inter_kk+Hmat_inter_kk';


            Hmat = zeros(tot_dim);
            % construct the hamiltonian
            Hmat=Hmat+opts.strain_fac*(Hmat_strain_L1+Hmat_strain_L2)+opts.inter_fac*(Hmat_inter_qdep + Hmat_inter_kk);
            %Hmat=(Hmat+Hmat')/2;
            % intralayer part
            % layer 1
            for indh=1:num_hex
                Hmat(all_index_L1(:,indh),all_index_L1(:,indh))=Hmat(all_index_L1(:,indh),all_index_L1(:,indh))+Layer1_ham(shift_klist_L1(indh,:));
                Hmat(all_index_L2(:,indh),all_index_L2(:,indh))=Hmat(all_index_L2(:,indh),all_index_L2(:,indh))+Layer2_ham(shift_klist_L2(indh,:));
            end
            [newbands,indband]=sort(real(eig(Hmat)),'ascend');
            allbands1(:,indk)=newbands;

            
            % don't use all_kpts2 (from old vers, TR partner of all_kpts1)
            %{
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
            Hmat=Hmat+opts.strain_fac*(Hmat_strain_L1+Hmat_strain_L2)+opts.inter_fac*(Hmat_inter_qdep)+Hmat_inter_kk*1;
            %Hmat=(Hmat+Hmat')/2;
            % intralayer part
            % layer 1
            for indh=1:num_hex
                Hmat(all_index_L1(:,indh),all_index_L1(:,indh))=Hmat(all_index_L1(:,indh),all_index_L1(:,indh))+Layer1_ham(shift_klist_L1(indh,:));
                Hmat(all_index_L2(:,indh),all_index_L2(:,indh))=Hmat(all_index_L2(:,indh),all_index_L2(:,indh))+Layer2_ham(shift_klist_L2(indh,:));
            end
            [newbands,indband]=sort(real(eig(Hmat)),'ascend');



            allbands2(:,indk)=newbands;
            %}

        end

        % if we did high-symm k-sampling, make a plot of bandstructre
        if (opts.plot_on == 1 && opts.vf_only ~= 1 && opts.gam_only ~= 1 && opts.full_bz ~= 1)
            figure;
            hold on
            size(allbands1)
            plot(scale_axis1,allbands1(:,1:knum_tot/2),'r')
            plot(scale_axis1,allbands1(:,knum_tot/2+1:end),'k')
            axis([0 1 -opts.ax_m opts.ax_m])
            xlabel('k')
            ylabel('Energy (eV)')
            set(gca,'xtick',[])
            title(['relaxed \theta = ' num2str(tar_theta,'%.2f') '^\circ'])

            %band_frames(tar_theta_idx) = getframe(gcf);
        end

        if (opts.plot_on == 1 && opts.full_bz == 1)
            band_mesh_h = zeros(bz_n,bz_n);
            band_mesh_e = zeros(bz_n,bz_n);

            kx_mesh = zeros(bz_n,bz_n);
            ky_mesh = zeros(bz_n,bz_n);

            nb = size(allbands1,1);
            bz_n = sqrt(size(allbands1,2));
            idx = 1;
            for y = 1:bz_n
                for x = 1:bz_n
                    band_mesh_h(x,y) = allbands1(nb/2,idx);
                    band_mesh_e(x,y) = allbands1(nb/2+1,idx);
                    kx_mesh(x,y) = all_kpts1(idx,1);
                    ky_mesh(x,y) = all_kpts1(idx,2);
                    idx = idx+1;
                end
            end

            clf
            hold on
            surf(kx_mesh,ky_mesh,band_mesh_h,'EdgeColor','none')
            surf(kx_mesh,ky_mesh,band_mesh_e,'EdgeColor','none')

            axis equal        
        end

        sweep_vals{tar_theta_idx} = allbands1;
        scaleaxis{tar_theta_idx} = scale_axis1;
        sweep_kpts{tar_theta_idx} = all_kpts1;
        fprintf("%d / %d done with theta sweep \n",tar_theta_idx,length(tar_theta_list));

    end
    
end

%

%filename = 'extended_kp_full_bz_0p505.mat';
%save(filename,'allbands1','all_kpts1','scale_axis1','tar_theta');