function [h_filename, postproc_filename] = rewannier_aapm_5band(filename,tar_folder,tar_theta,numk,k_radius,side,filetype)

    %clearvars;
    % Re-Wannierization of the effective TB

    if (side == 0)
        side_str = 'top';
        tar_bands = 1:2; %active bands
        dis_bands = 1:4;
    elseif (side == 1)
       side_str = 'bot'; 
       tar_bands = 4:5; %active bands
       dis_bands = 2:5;
    end
    
    h_filename = [tar_folder '/hmat_5band_' side_str '_' strrep(num2str(tar_theta,'%.2f'),'.','p') '_rewan.dat'];
    postproc_filename = [tar_folder '/hmat_5band_' side_str '_' strrep(num2str(tar_theta,'%.2f'),'.','p') '_rewan_postproc.mat'];
    
    num_orbs = 5;
    %tar_theta = 1.1;

    Nk = 101;
    NkGrid = numk;%300;

    %k_radius = 0.25; % radius around Gamma used in disentangling (ratio to reciprocal lattice constant)

    Err = 1E-12; % error tolerance
    %% Load bond data
    %tar_folder = 'hmats_8band';
    %tar_filename = [tar_folder '/hmat_8band_' strrep(num2str(tar_theta,'%.2f'),'.','p') '_RS.dat']; 
    % input is the real-space symmetrized model

    
    if (filetype == 0)

        fid = fopen(filename);

        temp = fgetl(fid);  % skip: Moire lattice vectors
        temp = str2num(fgetl(fid));  
        moire_L_x1 = [temp(1) temp(2) 0];
        temp = str2num(fgetl(fid));  
        moire_L_x2 = [temp(1) temp(2) 0];

        temp = fgetl(fid);  % skip: Moire reciprocal vectors
        temp = str2num(fgetl(fid));  
        moire_k_vec1 = [temp(1) temp(2) 0];
        temp = str2num(fgetl(fid));  
        moire_k_vec2 = [temp(1) temp(2) 0];

        temp = fgetl(fid);  % skip: Orbital locations
        for o = 1:num_orbs
            temp = str2num(fgetl(fid));  
            all_wan_xyz(o,:) = [temp(1) temp(2)];
        end

        temp = fgetl(fid);  % skip: Hamiltonian
        temp = fgetl(fid);  % skip: R_x R_y ... (etc)
        idx = 1;
        while (1)
            temp = fgetl(fid); % load the Hamiltonian data lines!
            if (temp == -1)
                break
            end
            H(idx,:) = str2num(temp);  
            idx = idx+1;
        end

        fclose(fid);
        H_in = H;
    elseif (filetype == 1)
        
        H_in = load(filename);
        
    end


    %% Decompose into harmonics
    % e.g. build the Bloch theory for the given TBH

    bonds = zeros(size(H_in,1),5);
    bonds(:,1:5) = H_in(:,1:5);
    bonds(:,5) = bonds(:,5) + 1j*H_in(:,6);

    [Har_list,~,Har_ind] = unique(bonds(:,1:2),'rows');

    H_Har = zeros(num_orbs,num_orbs,max(Har_ind));
    for brun = 1:length(Har_ind)
        bond_now = bonds(brun,3:5);

        H_Har(bond_now(1),bond_now(2),Har_ind(brun)) = ...
            H_Har(bond_now(1),bond_now(2),Har_ind(brun)) + bond_now(3);
        % provision for repeated bonds although it shouldn't happen
    end

    Hermitianize = @(h) (h+h')/2;
    Hk = @(g) Hermitianize(sum(...
                bsxfun(@times,H_Har,reshape(exp(-1j*2*pi*Har_list*g.'),1,1,[])),...
                    3)); 


    % g = [g1,g2] when the physical momentum is g1 b1 + g2 b2, where b1,b2 are
    % reciprocal lattice vectors

    %% E plot
    kk1 = [0,0,0];
    kk2 = [-1/2,0,0];
    kk3 = [-1/3,1/3,0];


    %kscan_list=[kk1;kk2;kk3;kk1];
    kscan_list=[kk1;kk2;kk3;kk1];



    [ all_kpts, scale_axis] = generate_k_line( Nk, kscan_list);

    all_kpts = all_kpts(:,1:2);



    bands = zeros(num_orbs,size(all_kpts,1));
    for krun = 1:size(all_kpts,1)

        Hk_now = Hk(all_kpts(krun,:));
        bands(:,krun) = sort(eig(Hk_now));
    end


    %% re-Wannierize

    kDis = @(k) norm(k*[sqrt(3),1;-sqrt(3),1]/2); 
    % length of a k point in units of reciprocal lattice constant


    % Symmetric sampling around Gamma
    [k1, k2] = meshgrid((0:(NkGrid-1))/NkGrid-1/2,(0:(NkGrid-1))/NkGrid-1/2);
    
    % For computing perturbation terms
    %allk12=linspace(0,1,NkGrid+1);
    %allk12=allk12(1:NkGrid);
    %[k1,k2]=meshgrid(allk12,allk12);
    
    % kX = k1*b1(1) + k2*b2(1);
    % kY = k1*b1(2) + k2*b2(2);
    % 
    % kGridSampPts = [kX(:)*kScale+kCen(1), kY(:)*kScale+kCen(2)];

    kGridSampPts = [k1(:),k2(:)];
    NkGridPts = size(kGridSampPts,1);

    SkStore = zeros(NkGridPts,1);


    PDOS = zeros(NkGridPts*num_orbs,5,2);
    % PDOS[x,j,m]:
    %       x = [k_index*num_orbs + orbital_index]
    %       j: 1 = Energy eigenvalue, 2:5 = PDOS_sum value (see below)
    %       m: 1 = original model, 2 = new model (AA_pm at flat bands)

    PDOS_sum = zeros(3,num_orbs);
    PDOS_sum(1,1:2)=[1,1];  % AA_pm
    PDOS_sum(2,3)=1;        % AA_z
    PDOS_sum(3,4:5)=[1,1];  % AB/BA

    Hk_store = zeros(num_orbs,num_orbs,NkGridPts);


    print_div = 20;

    for krun = 1:NkGridPts


        if (mod(krun,NkGridPts/print_div) == 0)
           pc = (100/print_div)*floor(krun/(NkGridPts/print_div));
           fprintf("%d%% done with ReWan \n",floor(pc)); 
        end

        k_now = kGridSampPts(krun,:);
        Hk_now = Hk(k_now);
        [Psi_orig,E0] = H_eig(Hk_now,Err);

        % PDOS-original
        PDOS((krun-1)*num_orbs + (1:num_orbs),1,1) = E0;
        PDOS((krun-1)*num_orbs + (1:num_orbs),2:4,1) = (PDOS_sum*power(abs(Psi_orig),2)).'/NkGridPts;

        % Re-Wannier

        wk = exp(-(kDis(k_now)/k_radius)^2);
        w_mat = -ones(1,num_orbs);
        w_mat(dis_bands) = wk; % disentanglement bands are weighted by the
                               % current distance from Gamma, e.g. only have
                               % noticable weight near Gamma.

        w_mat(tar_bands) = 1;  % but we want full weighting for the flat bands
        w_mat(w_mat==-1) = [];

        PsiNow = Psi_orig(:,dis_bands)*diag(w_mat);

        [U,Sig,~] = svd(PsiNow(1:2,:)'); % projection to AA_\pm
        PsiNow = PsiNow*U(:,1:2);
        Sig = diag(Sig);

        [PsiNew,~,~] = svd(PsiNow); % re-orthonormalize

        Psi_tar = PsiNew(:,1:2);
        Psi_perp = PsiNew(:,3:num_orbs);



        Sk0 = Psi_tar(1:2,:)';
        Sk1 = Psi_perp(3:num_orbs,:)';


        [U0,Sig0,V0] = svd(Sk0);
        [U1,Sig1,V1] = svd(Sk1);

        Sig0 = diag(Sig0);
        Sig1 = diag(Sig1);

        SkStore(krun) = min([Sig;Sig0;Sig1]);

        U_Wan = PsiNew*blkdiag(U0*V0',U1*V1');

        % PDOS-new
        PDOS((krun-1)*num_orbs + (1:num_orbs),1,2) = E0;
        PDOS((krun-1)*num_orbs + (1:num_orbs),2:4,2) =...
                            (PDOS_sum*power(abs(U_Wan'*Psi_orig),2)).'/NkGridPts;


        Hk_store(:,:,krun) = U_Wan'*Hk_now*U_Wan;
        U_store(:,:,krun) = U_Wan;

    end

    %% PDOS
    Nbins = 3000;

    PDOS_plot0 = zeros(Nbins,2,3);
    PDOS_plot1 = zeros(Nbins,2,3);
    for channel = 1:3
        PDOS_plot0(:,:,channel) = DOS_Fn(PDOS(:,[1,channel+1],1),Nbins);
        PDOS_plot1(:,:,channel) = DOS_Fn(PDOS(:,[1,channel+1],2),Nbins); 
    end

    tri_weights0 = DOS_Fn(PDOS(:,[1,2],1),500);
    tri_weights1 = DOS_Fn(PDOS(:,[1,2],2),500);

    % weight containment
    %EMin = mean([min(bands(tar_bands(1),:)),max(bands(tar_bands(1)-1,:))]);
    %EMax = mean([min(bands(tar_bands(2)+1,:)),max(bands(tar_bands(2),:))]);
    EMin = min(bands(tar_bands(1),:));
    EMax = max(bands(tar_bands(2),:));

    tri_weights0_trunc = tri_weights0(tri_weights0(:,1)>EMin,:);
    tri_weights0_trunc = tri_weights0_trunc(tri_weights0_trunc(:,1)<EMax,:);

    tri_weights1_trunc = tri_weights1(tri_weights1(:,1)>EMin,:);
    tri_weights1_trunc = tri_weights1_trunc(tri_weights1_trunc(:,1)<EMax,:);


    %% Convert to real-space             
    Hr = @(r_vec) sum(...
                bsxfun(@times,Hk_store,reshape(exp(1j*2*pi*kGridSampPts*r_vec.'),1,1,[])),...
                    3)/NkGridPts;

    rMax = 5;            
    [r1, r2] = meshgrid(-rMax:rMax,-rMax:rMax);   
    Har_list_new = [r1(:),r2(:)];
    NrPts = size(Har_list_new,1);

    H_Har_new = zeros(num_orbs,num_orbs,NrPts);
    for rrun = 1:NrPts
        H_Har_new(:,:,rrun) = Hr(Har_list_new(rrun,:));
    end

    % range analysis
    a1 = [1/2,sqrt(3)/2];
    a2 = [-1/2,sqrt(3)/2];

    basis_vecs = [ [0,0]; [0,0]; [0,0];...
                    (a1+a2)/3; -(a1+a2)/3;...
                    (a1-a2)/2; a2/2; -a1/2];

    [t_list0,dis_list0] = t_sort([a1;a2], basis_vecs, H_Har,Har_list,1E-8);
    [~,dis_list1] = t_sort([a1;a2], basis_vecs, H_Har_new,Har_list_new,1E-8);

    % truncate to the same max bond distance as original
    H_Har_new(dis_list1(:)>max(t_list0(:,1)))=0;
    [t_list1,dis_list1] = t_sort([a1;a2], basis_vecs, H_Har_new,Har_list_new,1E-8);
    
    
    H_save = zeros(num_orbs*num_orbs*NrPts,6);
    idx = 1;
    for rrun = 1:NrPts
        for o1 = 1:num_orbs
            for o2 = 1:num_orbs
                H_save(idx,1) = Har_list_new(rrun,1);
                H_save(idx,2) = Har_list_new(rrun,2);
                H_save(idx,3) = o1;
                H_save(idx,4) = o2;
                H_save(idx,5) = real(H_Har_new(o1,o2,rrun));
                H_save(idx,6) = imag(H_Har_new(o1,o2,rrun)); 
                idx = idx+1;
            end
        end
    end
    
    %% New Ek
    Hk_Wan = @(g) Hermitianize(sum(...
                bsxfun(@times,H_Har_new,reshape(exp(-1j*2*pi*Har_list_new*g.'),1,1,[])),...
                    3)); 

    bands_Wan = zeros(num_orbs,size(all_kpts,1));
    for krun = 1:size(all_kpts,1)
        Hk_now = Hk_Wan(all_kpts(krun,:));
        bands_Wan(:,krun) = sort(eig(Hk_now));
    end

    %% Save
    
    save(h_filename,'H_save','-ascii');
    save(postproc_filename, ...
        'U_store', 'kGridSampPts', ...
        'SkStore','scale_axis','bands','bands_Wan','t_list0','t_list1', ...
        'PDOS_plot0','PDOS_plot1','EMin','EMax', ...
        'tri_weights0','tri_weights1','tri_weights0_trunc','tri_weights1_trunc');



    %% Plots
    %{
    fprintf('Worst singular value = %g \n',min(SkStore));

    figure;
    subplot(1,2,1);
    hold on
        plot(scale_axis,bands','--k');
        plot(scale_axis,bands_Wan','-r');
    hold off
    xticks([]);
    title('Before (black) vs After (red)');
    ylabel('E_k (eV)');
    xlabel('k')
    set(gca,'fontsize',20);


    subplot(1,2,2);
    hold on
        scatter(t_list0(:,1),log(t_list0(:,2)),'ko');
        scatter(t_list1(:,1),log(t_list1(:,2)),'rs','filled');
    hold off
    legend('Before','After');
    ylim([-12,0]);
    title('Hopping strength vs distance');
    ylabel('log(|t_{ij}|)');
    xlabel('moire lattice costant')
    set(gca,'fontsize',20);

    set(gcf, 'PaperPosition', [0 0 16 9]);
    set(gcf,'PaperSize',[16,9]);
    %saveas(gcf,'Before_after','pdf');


    %%
    figure;
    plot_style = {'--xr','--sk','--ob','--dg'};


    subplot(3,1,1);
    hold on
        for channel = 4:-1:1
            plot(PDOS_plot0(:,1,channel),PDOS_plot0(:,2,channel),plot_style{channel});
        end
    hold off
    xlim([EMin, EMax]);
    ylim([0, max(max(PDOS_plot0(:,2,:)))]);
    title('Before');
    set(gca,'fontsize',20);

    subplot(3,1,2);
    hold on
        for channel = 4:-1:1
            plot(PDOS_plot1(:,1,channel),PDOS_plot1(:,2,channel),plot_style{channel});
        end
    hold off
    xlim([EMin, EMax]);
    ylim([0, max(max(PDOS_plot1(:,2,:)))]);
    title('After');
    set(gca,'fontsize',20);

    subplot(3,1,3);
    hold on
        for channel = 4:-1:1
            plot(PDOS_plot1(:,1,channel),PDOS_plot1(:,2,channel),plot_style{channel});
        end
    hold off
    xlim([EMin, EMax]);
    ylim([0,0.003]);
    title('After (zoomed)');
    xlabel('E (eV)')

    set(gca,'fontsize',20);

    set(gcf, 'PaperPosition', [0 0 16 9]);
    set(gcf,'PaperSize',[16,9]);
    %saveas(gcf,'PDOS','pdf');
    %%

    figure;
    hold on
        plot(tri_weights0(:,1),tri_weights0(:,2),'--xk');
        plot(tri_weights1(:,1),tri_weights1(:,2),'--sr');
    hold off
    legend('Before','After');
    ylim([0,0.01]);
    xlim([tri_weights0(1,1),tri_weights0(end,1)]);
    xlabel('E (eV)')
    title(sprintf('Weight contained %.2f \\rightarrow %.2f \n',...
        sum(tri_weights0_trunc(:,2))/2,sum(tri_weights1_trunc(:,2))/2));

    set(gca,'fontsize',20);
    set(gcf, 'PaperPosition', [0 0 16 9]);
    set(gcf,'PaperSize',[16,9]);
    %saveas(gcf,'PDOS_AApm','pdf');
    %}
end

%% Functions
function [PsiNow,ENow] = H_eig(H_Now,Err)
    % assuming H_Now is manifestly Hermitian
    [PsiNow,ENow] = eig(H_Now);
    ENow = diag(ENow);
    
    if norm(ENow-sort(ENow))>Err || norm(PsiNow*PsiNow'-eye(length(ENow)))>Err
        fprintf('issue with diagonalization. \n');
        return
        % could replace by E reordering and Psi orthogonalization
    end

end


function [DOS_out] = DOS_Fn(E_weights,Nbins)
    % E_weights a Nx2 matrix with the first column be energies and 2nd
    % the weights
    
    [~,EOrd] = sort(E_weights(:,1));
    E_weights = E_weights(EOrd,:);

    dE = (E_weights(end,1)-E_weights(1,1))/(Nbins-1);

    ERef = E_weights(1,1) + dE/2;
    DOS_out = zeros(Nbins,2);

    cnt = 1;
    DOS_out(cnt,1) = ERef-dE/2;
    for ERun = 1:size(E_weights,1)
        while E_weights(ERun,1) > ERef
            cnt = cnt+1;
            ERef = ERef + dE;
            DOS_out(cnt,1) = ERef-dE/2;
        end

        DOS_out(cnt,2) = DOS_out(cnt,2) + E_weights(ERun,2);
    end   
end


function [t_list,dis_list] = t_sort(a_vecs, basis, Hr,r_list,thres)
    % Each row of a_vecs is a lattice vec
    % Each row of basis is the position of the corresponding orbital in the
    % home cell
    % size(r_list,1) = size(Hr,3) should be true
    
    t_list = zeros(length(Hr(:)),2);
    num_orbs = size(Hr,1);
    dis_list = zeros(size(Hr));
    
    cnt = 0;
    for rrun = 1:size(r_list,1)
        r0 = r_list(rrun,:)*a_vecs;
        for to = 1:num_orbs
            for fr = 1:num_orbs
                cnt = cnt+1;
                t_list(cnt,1) = norm(r0+basis(to,:)-basis(fr,:));
                t_list(cnt,2) = abs(Hr(to,fr,rrun));
                
                dis_list(to,fr,rrun) = t_list(cnt,1);
            end
        end 
    end
    
    [~,re_ord] = sort(t_list(:,1));
    t_list = t_list(re_ord,:);
    
    t_list(t_list(:,2)<thres,:) = [];
end