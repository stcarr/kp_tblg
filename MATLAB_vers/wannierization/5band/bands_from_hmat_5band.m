
clear all;

% plot the symmetry enforced (non-symm.) bands in black (red)
% just pick folder location, hamiltonian size, and twist angle...

num_orbs = 5;
tar_theta = 1.1;
side = 1;
if (side == 0)
    side_text = 'top';
elseif (side == 1)
    side_text = 'bot';
end

tar_folder = 'run_folder_5band_final_07-12-2019';
tar_filename = [tar_folder '/hmat_5band_' side_text '_' strrep(num2str(tar_theta,'%.2f'),'.','p') '.dat'];

% everything after this is automated!

fid = fopen(tar_filename);

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

% generate the band structure

kk1=0*moire_k_vec1;
kk2=-moire_k_vec1/2;
kk3=-(moire_k_vec1-moire_k_vec2)/3;


%kscan_list=[kk1;kk2;kk3;kk1];
kscan_list=[kk3;kk1;kk2;kk3];

knum=31;

[ all_kpts, scale_axis] = generate_k_line( knum, kscan_list);

knum_tot=size(all_kpts);
knum_tot=knum_tot(1);
allbands_sym=zeros(5,knum_tot);
allbands_nosym=zeros(5,knum_tot);

%knum_tot=1;
%all_kpts=norm(moire_k_vec1)*[0.3,0.17,0];

phiphase=exp(i*2*pi/3);


Mmat=zeros(5);
Mmat(1:2,1:2)=[0,i;-i,0];
if (side == 0)
    Mmat(3,3)=1;
    Mmat(4,4)=-1;
    Mmat(5,5)=-1;
elseif (side == 1)
    Mmat(3,3)=-1;
    Mmat(4,4)=1;
    Mmat(5,5)=1;
end

Rmat=diag([phiphase',phiphase,1,1,1]);


C2Tmat=zeros(5);
C2Tmat(1:2,1:2)=[0,1;1,0];

if (side == 0)
    C2Tmat(3,3)=1;
    C2Tmat(4:5,4:5)=[0,-1;-1,0];
elseif( side == 1)
    C2Tmat(3,3)= 1;
    C2Tmat(4:5,4:5)=[0,1;1,0];
end

for indk=1:knum_tot
    know0=all_kpts(indk,1:2);


    sym_know=zeros(6,2);
    sym_know(1,:)=know0;
    % counterclock
    sym_know(2,:)=[-0.5*know0(1)-0.5*sqrt(3)*know0(2),0.5*sqrt(3)*know0(1)-0.5*know0(2)];
    sym_know(3,:)=[-0.5*know0(1)+0.5*sqrt(3)*know0(2),-0.5*sqrt(3)*know0(1)-0.5*know0(2)];

    sym_know(4:6,1)=-sym_know(1:3,1);
    sym_know(4:6,2)=sym_know(1:3,2);

    allHmat=zeros(5,5,6);
    allHmat_sym=zeros(5,5,6);

    % generate the set of Hamiltonians with k vectors related by symmetry
    for indr=1:size(H,1)

        for indks=1:6
            know = sym_know(indks,1:2);
            rvec_now = moire_L_x1(1:2)*H(indr,1) + moire_L_x2(1:2)*H(indr,2);
            m = H(indr,3);
            n = H(indr,4);
            t = H(indr,5) + 1j*H(indr,6);
 
            allHmat(m,n,indks)=allHmat(m,n,indks)+t*exp(-1j*dot(know,rvec_now));

        end
    end

    Hmat_nosym = allHmat(:,:,1);


    %all_gaugemat=zeros(5,5,6);

    % perform the necessary gauge transformation and 
    % enforce the C2T symmetry

    for indks=1:6
        know=sym_know(indks,1:2);
        %all_gaugemat(:,:,indks)
        tmpgauge=diag([1,1,1,exp(i*dot(know,all_wan_xyz(4,1:2))),exp(i*dot(know,all_wan_xyz(5,1:2)))]);
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

    allbands_nosym(:,indk)=sort(real(eig(squeeze(allHmat(:,:,1)))),'ascend');
    %allbands_nosym(:,indk)=sort(real(eig(Hmat_nosym)),'ascend');

end

bands = allbands_sym;

clf
subplot(1,2,1)
hold on
plot(scale_axis,bands','k')
plot(scale_axis,allbands_nosym','r')
