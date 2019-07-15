function [bands, scale_axis, allbands_nosym] = Generate_8band_TBH_sym_ver4(filename)

    %clear all;

    %tar_folder = 'test_runs_03-22-2019';

    %load('eightband_25x25k_1p00_deg_Hk.mat')
    %load([tar_folder '/eightWANfull_H_1p10_15k.mat']);
    load(filename);
    
    % WAN 1,2,3: at 0
    % WAN 4 at (moire_L_x1+moire_L_x2)/3;
    % WAN 5 at -(moire_L_x1+moire_L_x2)/3;

    %moire_L_x1
    %moire_L_x2

    all_wan_xyz=zeros(5,3);
    all_wan_xyz(4,:)=(moire_L_x1+moire_L_x2)/3;
    all_wan_xyz(5,:)=-(moire_L_x1+moire_L_x2)/3;

    all_wan_xyz(6,:)=(moire_L_x1-moire_L_x2)/2;
    all_wan_xyz(7,:)=(moire_L_x2)/2;
    all_wan_xyz(8,:)=(-moire_L_x1)/2;


    %all_wan_xyz



    %{
    figure(1);
    clf

    scatter3(all_kpts(:,1),all_kpts(:,2),imag(all_new_hmat(1,4,:)));

    axis([-inf,inf,-inf,inf,-inf,inf]);
    %}

    %

    hop_cut=norm(moire_L_x1)*9.9;
    %hop_cut

    MMM=20;

    clear hex_R;
    indh=1;
    for ind1=(-MMM):MMM
        for ind2=(-MMM):MMM
            vectmp=ind1*moire_L_x1+ind2*moire_L_x2;


            if norm(vectmp)<=hop_cut
                hex_R(indh,:)=vectmp(:);
                indh=indh+1;

            end


        end
    end
    num_hex=indh-1;


    %scatter(hex_R(:,1),hex_R(:,2))

    BZ_area=cross(moire_k_vec1,moire_k_vec2);
    BZ_area=abs(BZ_area(3));

    dBZ_area=BZ_area/numk/numk;


    all_R_mat=zeros(8,8,num_hex);



    for indr=1:num_hex
        R_hop=hex_R(indr,1:2);
        for indk=1:knum_tot
            know=all_kpts(indk,1:2);

            Heff_now=squeeze(all_new_hmat(1:8,1:8,indk));

            all_R_mat(:,:,indr)=all_R_mat(:,:,indr)+Heff_now*exp(1j*dot(know,R_hop))/(numk^2);

        end

        tmp=squeeze(all_R_mat(:,:,indr));
        maxRR(indr)=sum(abs(tmp(:)));

    end


    % scatter3(hex_R(:,1),hex_R(:,2),abs(all_R_mat(5,4,:)));

    % symmetric cutoff
    TBH_Rcut=norm(moire_L_x1)*6;

    Rcut_mask=zeros(8,8,num_hex);

    for indr=1:num_hex
        bigR=hex_R(indr,1:2);
        for indp1=1:8
            pos1=all_wan_xyz(indp1,1:2);
            for indp2=1:8
                pos2=all_wan_xyz(indp2,1:2);


                r_connect=bigR+pos2-pos1;
                if norm(r_connect)<TBH_Rcut
                    Rcut_mask(indp2,indp1,indr)=1;
                end

            end

        end
    end

    all_R_mat_cut=zeros(8,8,num_hex);
    for indr=1:num_hex
        tmph=squeeze(all_R_mat(:,:,indr));
        tmpm=squeeze(Rcut_mask(:,:,indr));

        all_R_mat_cut(:,:,indr)=tmph.*tmpm;


    end



    % know0=[2,3];
    % know1=[-0.5*know0(1)-0.5*sqrt(3)*know0(2),0.5*sqrt(3)*know0(1)-0.5*know0(2)];
    % know2=[-0.5*know0(1)+0.5*sqrt(3)*know0(2),-0.5*sqrt(3)*know0(1)-0.5*know0(2)];
    %     
    %     
    % norm(know0)
    % norm(know1)
    % norm(know2)
    % 
    % dot(know0,know1)/norm(know0)/norm(know0)
    % dot(know1,know2)/norm(know0)/norm(know0)
    % dot(know2,know0)/norm(know0)/norm(know0)



    % generate the band strcuture


    kk1=0*moire_k_vec1;
    kk2=-moire_k_vec1/2;
    kk3=-(moire_k_vec1-moire_k_vec2)/3;


    %kscan_list=[kk1;kk2;kk3;kk1];
    kscan_list=[kk3;kk1;kk2;kk3];

    knum=31;

    [ all_kpts, scale_axis] = generate_k_line( knum, kscan_list);

    knum_tot=size(all_kpts);
    knum_tot=knum_tot(1);
    allbands_sym=zeros(8,knum_tot);
    allbands_nosym=zeros(8,knum_tot);

    %knum_tot=1;
    %all_kpts=norm(moire_k_vec1)*[0.3,0.17,0];

    phiphase=exp(1j*2*pi/3);


    % old symmetrization code (in momentum space), not used anymore.
    % likely inconsintent with the WF Rephase in real_space_sym.
    Mmat=zeros(8);
    Mmat(1:2,1:2)=[0,i;-i,0];
    Mmat(3,3)=1;
    Mmat(4,4)=-1;
    Mmat(5,5)=-1;
    Mmat(6:8,6:8)=[1,0,0;0,0,1;0,1,0];

    Rmat=zeros(8);
    Rmat(1:5,1:5)=diag([phiphase',phiphase,1,1,1]);
    Rmat(6:8,6:8)=[0,1,0;0,0,1;1,0,0];

    C2Tmat=zeros(8);
    C2Tmat(1:2,1:2)=[0,1;1,0];
    C2Tmat(3,3)=1;
    C2Tmat(4:5,4:5)=[0,-1;-1,0];
    C2Tmat(6:8,6:8)=eye(3);

    for indk=1:knum_tot
        know0=all_kpts(indk,1:2);


        sym_know=zeros(6,2);
        sym_know(1,:)=know0;
        % counterclock
        sym_know(2,:)=[-0.5*know0(1)-0.5*sqrt(3)*know0(2),0.5*sqrt(3)*know0(1)-0.5*know0(2)];
        sym_know(3,:)=[-0.5*know0(1)+0.5*sqrt(3)*know0(2),-0.5*sqrt(3)*know0(1)-0.5*know0(2)];

        sym_know(4:6,1)=-sym_know(1:3,1);
        sym_know(4:6,2)=sym_know(1:3,2);

        allHmat=zeros(8,8,6);
        allHmat_sym=zeros(8,8,6);

        % generate the set of Hamiltonians with k vectors related by symmetry
        for indr=1:num_hex

            for indks=1:6
                know=sym_know(indks,1:2);
                rvec_now=hex_R(indr,1:2);

                allHmat(:,:,indks)=allHmat(:,:,indks)+squeeze(all_R_mat_cut(:,:,indr))*exp(-1j*dot(know,rvec_now));

            end
        end


        Hmat_nosym=squeeze(allHmat(:,:,1));

        %all_gaugemat=zeros(5,5,6);

        % perform the necessary gauge transformation and 
        % enforce the C2T symmetry

        for indks=1:6
            know=sym_know(indks,1:2);
            %all_gaugemat(:,:,indks)
            %tmpgauge=diag([1,1,1,exp(i*dot(know,all_wan_xyz(4,1:2))),exp(i*dot(know,all_wan_xyz(5,1:2)))]);
            tmpgauge0=all_wan_xyz(:,1:2)*know';
            tmpgauge=diag(exp(i*tmpgauge0));

            tmpH=squeeze(allHmat(:,:,indks));
            tmpH2=tmpgauge'*tmpH*tmpgauge;
            C2TtmpH=(tmpH2+C2Tmat*conj(tmpH2)*C2Tmat)/2;
            allHmat(:,:,indks)=C2TtmpH;


        end

        % ensure the Mirror and R3 rotation symmetries
        allHmat_sym(:,:,1)=squeeze(allHmat(:,:,1));
        allHmat_sym(:,:,2)=Rmat*squeeze(allHmat(:,:,2))*Rmat';
        allHmat_sym(:,:,3)=Rmat'*squeeze(allHmat(:,:,3))*Rmat;
        allHmat_sym(:,:,4)=Mmat'*squeeze(allHmat(:,:,4))*Mmat;
        allHmat_sym(:,:,5)=Rmat*Mmat'*squeeze(allHmat(:,:,5))*Mmat*Rmat';
        allHmat_sym(:,:,6)=Rmat'*Mmat'*squeeze(allHmat(:,:,6))*Mmat*Rmat;


        % symmetric Hamiltonian is then derived by the average
        Hmatsym=mean(allHmat_sym,3);
        allbands_sym(:,indk)=sort(real(eig(Hmatsym)),'ascend');
        allbands_nosym(:,indk)=sort(real(eig(Hmat_nosym)),'ascend');



    end
    
    bands = allbands_sym;


    %{

    pmat1=squeeze(allHmat_sym(:,:,4));
    pmat2=squeeze(allHmat_sym(:,:,6));
    mesh(abs(pmat1-pmat2))


    %}

    %{
    figure(1);
    clf
    subplot(1,2,1);
    hold on;
    plot(scale_axis,allbands_nosym','b','LineWidth',2);
    plot(scale_axis,allbands_sym','r','LineWidth',1);

    hold off;


    subplot(1,2,2);
    hold on;
    plot(scale_axis,allbands_nosym','b','LineWidth',2);
    plot(scale_axis,allbands_sym','r','LineWidth',1);

    hold off;
    axis([-inf,inf,-0.01,0.01])
    %}

    %{

    % check how gapless it is at K 
    tt=allbands_sym(:,1);
    tt(5)-tt(4)

    tt=allbands_nosym(:,1);
    tt(5)-tt(4)
    %}
    
end

