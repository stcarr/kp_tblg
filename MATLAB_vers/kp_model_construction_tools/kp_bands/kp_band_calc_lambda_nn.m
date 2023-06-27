clearvars -except global_var_file;

run(global_var_file);


% the twist angle converted into radian

% this MacDonald effective theory has magic angle ~1.1 degree

rot_theta = 1.1213/180*pi;

%rot_theta=1.1/180*pi;

%rot_theta=1.5/180*pi;

%rot_theta=0.63/180*pi;

 

% cutoff in momentum space, which would affect accuracies for the solution

hex_cut=4.51;

 

L1_mass_gap=0.0;

L2_mass_gap=0.0;

 

lattice_a=1.42*sqrt(3);

 

KD=4*pi/3/lattice_a;

 

KTH=2*KD*sin(rot_theta/2);

 

bilayer_w0=0.10;
bilayer_w1=0.11;

% dft-inter at 28x27
bilayer_w0 = 0.075;
bilayer_w1 = 0.105;

% koshino-inter at 28x27
bilayer_w0 = 0.0761;
bilayer_w1 = 0.1151;

fermivf=1.03*1*2.7*sqrt(3)/2*lattice_a;

 

% basic variables

sigx=[0,1;1,0];

sigy=[0,-i;i,0];

sigz=[1,0;0,-1];

 

Lambda_kplus= 0*0.12*i;

Lambda_kminus=conj(Lambda_kplus);

 

 

ph1=exp(i*2*pi/3);

ph2=exp(-i*2*pi/3);

 

TT1_0=[bilayer_w0,bilayer_w1;bilayer_w1,bilayer_w0];

TT2_0=[bilayer_w0,bilayer_w1*ph2;bilayer_w1*ph1,bilayer_w0];

TT3_0=[bilayer_w0,bilayer_w1*ph1;bilayer_w1*ph2,bilayer_w0];

 

% TT1_0=[bilayer_w0,bilayer_w0;bilayer_w0,bilayer_w0];

% TT2_0=[bilayer_w0,bilayer_w0*ph2;bilayer_w0*ph1,bilayer_w0];

% TT3_0=[bilayer_w0,bilayer_w0*ph1;bilayer_w0*ph2,bilayer_w0];

 

TT1_kplus=[1,1;1,1];

TT2_kplus=[ph1,1;ph2,ph1];

TT3_kplus=[ph2,1;ph1,ph2];

 

TT1_kminus=[1,1;1,1];

TT2_kminus=[ph2,ph1;1,ph2];

TT3_kminus=[ph1,ph2;1,ph1];

 

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

sqrt(dot(moire_L_x2,moire_L_x2))/10

 

 

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

 

%%

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

 

 

%%

 

all_index=reshape(1:tot_dim,2,tot_dim/2);

 

indL1=1:(tot_dim/2);

indL2=(tot_dim/2+1):tot_dim;

 

indL1A=indL1(1:2:(tot_dim/2));

indL1B=indL1(2:2:(tot_dim/2));

 

indL2A=indL2(1:2:(tot_dim/2));

indL2B=indL2(2:2:(tot_dim/2));

 

 

scan_klist=KTH*[[0.5,-0.5*sqrt(3),0];[0,0,0];[0,-0.5*sqrt(3),0];[-0.5,-0.5*sqrt(3),0]];%;[0.5/2,0.5/2*sqrt(3),0];[0,0,0]];

 

scan_klist2=KTH*[[0.5,0.5*sqrt(3),0];[0,0,0];[0,0.5*sqrt(3),0];[-0.5,0.5*sqrt(3),0]];%;[0.5/2,0.5/2*sqrt(3),0];[0,0,0]];

 

 

 

knum=51;

[ all_kpts, scale_axis] = generate_k_line( knum, scan_klist );

[ all_kpts2, scale_axis2] = generate_k_line( knum, scan_klist2 );

knum_tot=size(all_kpts);

knum_tot=knum_tot(1);

 

all_index=reshape(1:tot_dim,2,tot_dim/2);

 

 

allbands1=zeros(tot_dim,knum_tot);

allbands2=zeros(tot_dim,knum_tot);

 

testfac=0;

 

TT_inter_mat=kron(connect_Mat1,TT1_0)+kron(connect_Mat2,TT2_0)+kron(connect_Mat3,TT3_0);

TT_inter_kplus=kron(connect_Mat1,TT1_kplus)+kron(connect_Mat2,TT2_kplus)+kron(connect_Mat3,TT3_kplus);

TT_inter_kminus=kron(connect_Mat1,TT1_kminus)+kron(connect_Mat2,TT2_kminus)+kron(connect_Mat3,TT3_kminus);

 

 

for indk=1:knum_tot

 

    know=all_kpts(indk,1:2);

    

    shift_klist(:,1)=hex_all(:,1)+know(1);

    shift_klist(:,2)=hex_all(:,2)+know(2);

    

    Hmat = zeros(tot_dim);

    

    % construct the hamiltonian

    

    kplus_layer1=kron(diag(shift_klist(1:num_hex,1)+i*shift_klist(1:num_hex,2)),eye(2));

    kminus_layer1=kron(diag(shift_klist(1:num_hex,1)-i*shift_klist(1:num_hex,2)),eye(2));

    Hmat(indL2,indL1)=TT_inter_mat+Lambda_kplus*TT_inter_kplus*kplus_layer1+Lambda_kminus*TT_inter_kminus*kminus_layer1;

    

    

    % interlayer part

    

    Hmat=Hmat+Hmat';

    

    

    % layer 1

    for indh=1:num_hex

        Hmat(all_index(:,indh),all_index(:,indh))=gen_dirac_ham(shift_klist(indh,:),testfac*rot_theta/2)+L1_mass_gap*sigz;

        

    end

    

    % layer 2

    for indh=(num_hex+1):2*num_hex

        Hmat(all_index(:,indh),all_index(:,indh))=gen_dirac_ham(shift_klist(indh,:),-testfac*rot_theta/2)+L2_mass_gap*sigz;

        

    end

    

    

    allbands1(:,indk)=sort(real(eig(Hmat)),'ascend');

    

    

    

    know=all_kpts2(indk,1:2);

    

    shift_klist(:,1)=hex_all(:,1)+know(1);

    shift_klist(:,2)=hex_all(:,2)+know(2);

    

    Hmat = zeros(tot_dim);

    

    % construct the hamiltonian

    kplus_layer1=kron(diag(shift_klist(1:num_hex,1)+i*shift_klist(1:num_hex,2)),eye(2));

    kminus_layer1=kron(diag(shift_klist(1:num_hex,1)-i*shift_klist(1:num_hex,2)),eye(2));

    Hmat(indL2,indL1)=TT_inter_mat+Lambda_kplus*TT_inter_kplus*kplus_layer1+Lambda_kminus*TT_inter_kminus*kminus_layer1;

    

    

    % interlayer part

    Hmat=Hmat+Hmat';

    % layer 1

    for indh=1:num_hex

        Hmat(all_index(:,indh),all_index(:,indh))=gen_dirac_ham(shift_klist(indh,:),testfac*rot_theta/2)+L1_mass_gap*sigz;

        

    end

    

    % layer 2

    for indh=(num_hex+1):2*num_hex

        Hmat(all_index(:,indh),all_index(:,indh))=gen_dirac_ham(shift_klist(indh,:),-testfac*rot_theta/2)+L2_mass_gap*sigz;

        

    end

    allbands2(:,indk)=sort(real(eig(Hmat)),'ascend');

end

 

 

%%

 

figure(1);

 

sizetmp=size(allbands1);

tmps=sizetmp(2)/2;

 

ccmean1=mean(allbands1(tot_dim/2:(tot_dim/2+1),1));

ccmean2=mean(allbands2(tot_dim/2:(tot_dim/2+1),1));

 

 

 

plot(scale_axis,(allbands1'-ccmean1),'r','LineWidth',3);

 

hold on;

plot(scale_axis,(allbands2'-ccmean2),'b','LineWidth',2);

 

axis([-inf,inf,-0.15,0.15]);

%axis([-inf,inf,-.25,.25]);

ffmmm=0.03;

axis([-inf,inf,-ffmmm,ffmmm]);

 

set(gca,'FontSize',24);

 

ylabel('Energy (eV)')

box on;