clear all;

% set the "twist angle" with the supercell pair [M,N]:
sc_m=29;%[7:29, 30:2:38, 40:5:80];
sc_n=sc_m-1;

% turns on/off relaxation
fft_relax_onoff_factor = 1;

%  ---------------- !! NOTE !! ------------------------
% both layers have been "shifted" downwards a bit in order to make the
% rotation center is at a honeycomb center (this causes D6 symm, instead of
% D3 symmetry if rotated at an "AA" stacking dimer). See the variables
% "shift_bot_vec" and "shift_top_vec".


% setup basic environment for supercell
% graphene_setup_supercell.m

% input: sc_m, sc_n

% expect sc_m >=2    sc_n>=1


lattice_a=1.42*sqrt(3);

layer_d=[0,0,3.35];

% lattice_a*
a1=[sqrt(3)/2,-1/2]';
a2=[sqrt(3)/2,+1/2]';


sc_int_bottom_a1=[sc_n,sc_m];
sc_int_bottom_a2=[-sc_m,sc_n+sc_m];

sc_int_top_a1=[sc_m,sc_n];
sc_int_top_a2=[-sc_n,sc_n+sc_m];


num_pc=(sc_m^2+sc_m*sc_n+sc_n^2);
total_num=4*num_pc;

rot_theta=acos((sc_m^2+sc_n^2+4*sc_m*sc_n)/2/(sc_m^2+sc_m*sc_n+sc_n^2));
% show in degrees
fprintf("constructing %.2f degree TBH (%d atoms) \n", rot_theta/pi*180, total_num);

rot_mat=[cos(rot_theta),-sin(rot_theta);sin(rot_theta),cos(rot_theta)];

ra1=rot_mat*a1;
ra2=rot_mat*a2;

ra_mat=zeros(2);
ra_mat(:,1)=ra1;
ra_mat(:,2)=ra2;


% for both bottom and top unit
sc_t1=sc_n*a1+sc_m*a2;
sc_t2=-sc_m*a1+(sc_m+sc_n)*a2;


sc_ft1=[sc_t1(1),sc_t1(2),0];
sc_ft2=[sc_t2(1),sc_t2(2),0];
sc_ft3=[0,0,1];

sc_v=abs(dot(sc_ft1,cross(sc_ft2,sc_ft3)));
sc_b1=2*pi*cross(sc_ft2,sc_ft3)/sc_v;
sc_b2=2*pi*cross(sc_ft3,sc_ft1)/sc_v;
sc_b3=2*pi*cross(sc_ft1,sc_ft2)/sc_v;

sc_vec1=sc_b1(1:2);
sc_vec2=sc_b2(1:2)+sc_b1(1:2);
sc_gamma=sc_vec1*0;
sc_kpoint=(sc_vec1+sc_vec2)/3;
sc_mpoint=sc_vec1*0.5;

% primitive cell coordinate

pc_vec_a1=[a1',0];
pc_vec_a2=[a2',0];
pc_vec_a3=[0,0,1];

pc_vec_ra1=[ra1',0];
pc_vec_ra2=[ra2',0];
pc_vec_ra3=[0,0,1];


pc_vec_v=abs(dot(pc_vec_a1,cross(pc_vec_a2,pc_vec_a3)));
pc_vec_b1=2*pi*cross(pc_vec_a2,pc_vec_a3)/pc_vec_v;
pc_vec_b2=2*pi*cross(pc_vec_a3,pc_vec_a1)/pc_vec_v;
pc_vec_b3=2*pi*cross(pc_vec_a1,pc_vec_a2)/pc_vec_v;

pc_vec_rv=abs(dot(pc_vec_ra1,cross(pc_vec_ra2,pc_vec_ra3)));
pc_vec_rb1=2*pi*cross(pc_vec_ra2,pc_vec_ra3)/pc_vec_rv;
pc_vec_rb2=2*pi*cross(pc_vec_ra3,pc_vec_ra1)/pc_vec_rv;
pc_vec_rb3=2*pi*cross(pc_vec_ra1,pc_vec_ra2)/pc_vec_rv;

pc_bot_Mpoint=pc_vec_b1/2;
pc_bot_Kpoint=(2*pc_vec_b1+pc_vec_b2)/3;


pc_top_Mpoint=pc_vec_rb1/2;
pc_top_Kpoint=(2*pc_vec_rb1+pc_vec_rb2)/3;



% shift the layer
shift_top_vec=(-1/3)*(2*ra2-ra1);
shift_bot_vec=(-1/3)*(2*a2-a1);

% graphene_setup_supercell_2.m

% this is the file to establish atomic sites and the link between them
% both intra layer and inter layer

% input: sc_m, sc_n

code_check=false;


% expect sc_m >=2    sc_n>=1

% graphene_setup_supercell_env


% rescale now?

% inter_cutoff=3.99;

inter_cutoff=2.2975;


% bottom layer
inplane_vecs_bottom=zeros(3,2);
deltav_b=(a1+a2)/3;
inplane_vecs_bottom(1,:)=deltav_b;
inplane_vecs_bottom(2,:)=deltav_b-a1;
inplane_vecs_bottom(3,:)=deltav_b-a2;

% top layer
deltav_t=(ra1+ra2)/3;
inplane_vecs_top(1,:)=deltav_t;
inplane_vecs_top(2,:)=deltav_t-ra1;
inplane_vecs_top(3,:)=deltav_t-ra2;



ind=2;
all_pos_layer1=[0,0];
all_pos_layer2=[0,0];

% bottom layer atom positions

% boundary:     -M  to  N
% boundary:     0   to  N+2M

% actual boundary:     -M+1  to  N-1     (N+M-1)
% actual boundary:     0     to  N+2M-1  (N+2M)

% large boundary:     -M  to  N     (N+M+1)
% large boundary:     -1     to  N+2M  (N+2M+2)

offset_inds_bottom=[-sc_m+1,0];

index_mat_b=zeros(sc_m+sc_n+1,sc_n+2*sc_m+2);
inverse_index_mat_b=0;
index_mat_b(sc_m+1,2)=1; % first [0,0] Bravvis point
inverse_index_mat_b(1,1:2)=[0,0];

% go over all possible positions
for ind1=(-sc_m+1):(sc_n-1)
    for ind2=1:(sc_n+2*sc_m-1)
        
        % the point at (0,0) was already added, so skip
        specialpt = (ind1==0 && ind2==0);
        
        % keep vector positions in primitive supercell
        [ newvec ] = moveto_primitive_method2( [sc_n,sc_m],[-sc_m,sc_n+sc_m],[ind1,ind2]' );
        
        if newvec(1)==ind1 && newvec(2)==ind2  && ~specialpt
            tmpvec=a1*ind1+a2*ind2;
            inverse_index_mat_b(ind,1:2)=[ind1,ind2];
            all_pos_layer1(ind,1:2)=tmpvec(:);
            
            index_mat_b(ind1-(-sc_m+1)+2,ind2+2)=ind;
            
            ind=ind+1;
        end
        
    end
end

%sum(index_mat_b(:)>0)

num_pos_b=ind-1;


% top layer atom positions, same as above


% actual boundary:     -N+1  to  M-1     (N+M-1)
% actual boundary:     0     to  2N+M-1  (2N+M)

% large boundary:     -N  to  M     (N+M+1)
% large boundary:     -1     to  2N+M  (2N+M+2)

offset_inds_top=[-sc_n+1,0];

index_mat_t=zeros(sc_n+sc_m+1,2*sc_n+sc_m+2);
index_mat_t(sc_n+1,2)=1; % first [0,0] Bravvis point
inverse_index_mat_t=0;
ind=2;
inverse_index_mat_t(1,1:2)=[0,0];

for ind1=(-sc_n+1):(sc_m-1)
    for ind2=1:(2*sc_n+sc_m-1)
        
        specialpt = (ind1==0 && ind2==0);
        
        
        [ newvec ] = moveto_primitive_method2( [sc_m,sc_n],[-sc_n,sc_n+sc_m],[ind1,ind2]' );
        
        if newvec(1)==ind1 && newvec(2)==ind2  && ~specialpt
            tmpvec=ra1*ind1+ra2*ind2;
            
            inverse_index_mat_t(ind,1:2)=[ind1,ind2];
            
            all_pos_layer2(ind,1:2)=tmpvec(:);
            
            index_mat_t(ind1-(-sc_n+1)+2,ind2+2)=ind;
            
            ind=ind+1;
        end
        
    end
end

num_pos_t=ind-1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% all lattice points
tot_num=4*num_pos_t;
index_bot_atoms=1:(2*num_pos_b);
bot_A_atoms=1:2:(2*num_pos_b);
bot_B_atoms=2:2:(2*num_pos_b);
index_top_atoms=(2*num_pos_b+1):(4*num_pos_b);
top_A_atoms=(2*num_pos_b+1):2:(4*num_pos_b);
top_B_atoms=(2*num_pos_b+2):2:(4*num_pos_b);

sc_all_points=zeros(tot_num,2);

sc_all_points_shift=zeros(tot_num,2);

sc_all_points(bot_A_atoms,:)=all_pos_layer1(:,:);
sc_all_points(bot_B_atoms,:)=all_pos_layer1(:,:);
sc_all_points(bot_B_atoms,1)=sc_all_points(bot_B_atoms,1)+deltav_b(1);
sc_all_points(bot_B_atoms,2)=sc_all_points(bot_B_atoms,2)+deltav_b(2);

sc_all_points(top_A_atoms,:)=all_pos_layer2(:,:);
sc_all_points(top_B_atoms,:)=all_pos_layer2(:,:);
sc_all_points(top_B_atoms,1)=sc_all_points(top_B_atoms,1)+deltav_t(1);
sc_all_points(top_B_atoms,2)=sc_all_points(top_B_atoms,2)+deltav_t(2);

sc_all_points_shift=sc_all_points;

sc_all_points_shift(bot_A_atoms,1)=sc_all_points_shift(bot_A_atoms,1)+shift_bot_vec(1);
sc_all_points_shift(bot_A_atoms,2)=sc_all_points_shift(bot_A_atoms,2)+shift_bot_vec(2);
sc_all_points_shift(bot_B_atoms,1)=sc_all_points_shift(bot_B_atoms,1)+shift_bot_vec(1);
sc_all_points_shift(bot_B_atoms,2)=sc_all_points_shift(bot_B_atoms,2)+shift_bot_vec(2);

sc_all_points_shift(top_A_atoms,1)=sc_all_points_shift(top_A_atoms,1)+shift_top_vec(1);
sc_all_points_shift(top_A_atoms,2)=sc_all_points_shift(top_A_atoms,2)+shift_top_vec(2);
sc_all_points_shift(top_B_atoms,1)=sc_all_points_shift(top_B_atoms,1)+shift_top_vec(1);
sc_all_points_shift(top_B_atoms,2)=sc_all_points_shift(top_B_atoms,2)+shift_top_vec(2);




corner_pos=zeros(4,2);
corner_pos(2,:)=sc_t1;
corner_pos(3,:)=sc_t2;
corner_pos(4,:)=sc_t1+sc_t2;

% sc_all_points
% units of lattice constant
 
% these are in AA

pos_all_points=sc_all_points;


pos_all_points(:,3)=0;
pos_all_points=pos_all_points*lattice_a;

pos_all_points((tot_num/2+1):tot_num,3)=3.35;

pos_a1=sc_t1'*lattice_a;
pos_a2=sc_t2'*lattice_a;



if code_check
    
    for ind=1:7
        
        
        cc=actual_index_mat_b(:,:,ind);
        cc=cc(:);
        cc=sort(cc,'descend');
        %cc=cc(1:19);
        kk(ind,1)=sum(cc(:));
        % kk(ind,3)=prod(cc(:)+0.1);
        
        
        cc=actual_index_mat_t(:,:,ind);
        cc=cc(:);
        cc=sort(cc,'descend');
        %cc=cc(1:19);
        kk(ind,2)=sum(cc(:));
        % kk(ind,4)=prod(cc(:)+0.1);
        
    end
end


% graphene_setup_supercell_relax_ver2.m

% in Shiang's code, normally I use the length scale lattice_a = 1 !!!

% crystal relax:
newx=sc_b1/sqrt(dot(sc_b1,sc_b1));
newy=sc_t2/sqrt(dot(sc_t2,sc_t2));

% update atomic positions

sc_all_points_shift_relax=sc_all_points_shift;

% bottom layer
postmp1= Relax_model_Stephen_ver3(rot_theta,sc_b1(1:2)',sc_b2(1:2)',sc_t1,sc_t2,false,sc_all_points_shift(1:(tot_num/2),:));
% top layer
postmp2= Relax_model_Stephen_ver3(rot_theta,sc_b1(1:2)',sc_b2(1:2)',sc_t1,sc_t2,true,sc_all_points_shift((tot_num/2+1):(tot_num),:));

% bottom layer
Koshino_relax_shift(1:(tot_num/2),1:2)=fft_relax_onoff_factor*postmp1(:,1:2)/lattice_a;
% top layer
Koshino_relax_shift((tot_num/2+1):(tot_num),:)=fft_relax_onoff_factor*postmp2(:,1:2)/lattice_a;

sc_all_points_shift_relax=sc_all_points_shift_relax+Koshino_relax_shift;

height_relax=zeros(tot_num,1);

AB_Bernal_d=3.35;
% bottom layer (A)
height_relax(1:(tot_num/2))=fft_relax_onoff_factor*postmp1(:,3);
% top layer (A)
height_relax((tot_num/2+1):(tot_num))=AB_Bernal_d+fft_relax_onoff_factor*postmp2(:,3);

% plot the results

clf
hold on
scatter3(sc_all_points_shift_relax(1:tot_num/2,1),sc_all_points_shift_relax(1:tot_num/2,2),height_relax(1:tot_num/2),'.');
scatter3(sc_all_points_shift_relax(tot_num/2+1:end,1),sc_all_points_shift_relax(tot_num/2+1:end,2),height_relax(tot_num/2+1:end),'.');

plot([0 sc_t1(1)],[0 sc_t1(2)],'--k')
plot([0 sc_t2(1)],[0 sc_t2(2)],'--k')
plot(sc_t2(1)+[0 sc_t1(1)],sc_t2(2)+[0 sc_t1(2)],'--k')
plot(sc_t1(1)+[0 sc_t2(1)],sc_t1(2)+[0 sc_t2(2)],'--k')

view(2)
axis equal
