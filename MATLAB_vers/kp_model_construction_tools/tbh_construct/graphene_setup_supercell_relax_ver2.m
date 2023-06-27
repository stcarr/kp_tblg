
% input variables:
%
% in Shiang's code, normally I use the length scale lattice_a = 1 !!!



% crystal relax:

% 2.65

% numGG=3;
% Koshino_fft_Gs=zeros(numGG,2);
% 
% Koshino_fft_Gs(1,:)=sc_b1(1:2);
% Koshino_fft_Gs(2,:)=2*sc_b1(1:2);
% Koshino_fft_Gs(3,:)=2*sc_b1(1:2)+sc_b2(1:2);
% Koshino_fft_Gs(4,:)=3*sc_b1(1:2);
% Koshino_fft_Gs(5,:)=3*sc_b1(1:2)+sc_b2(1:2);
% Koshino_fft_Gs(6,:)=3*sc_b1(1:2)+2*sc_b2(1:2);
% Koshino_fft_Gs(7,:)=4*sc_b1(1:2)+0*sc_b2(1:2);
% Koshino_fft_Gs(8,:)=4*sc_b1(1:2)+1*sc_b2(1:2);
% Koshino_fft_Gs(9,:)=4*sc_b1(1:2)+2*sc_b2(1:2);
% Koshino_fft_Gs(10,:)=4*sc_b1(1:2)+3*sc_b2(1:2);


newx=sc_b1/sqrt(dot(sc_b1,sc_b1));
newy=sc_t2/sqrt(dot(sc_t2,sc_t2));


% A?
Koshino_r0=0.184*lattice_a;
% Koshino_pi=@(dd) -2.7*exp(-(dd*lattice_a-lattice_a/sqrt(3))/Koshino_r0);
% Koshino_sigma=@(dd) 0.48*exp(-(dd*lattice_a-layer_d(3))/Koshino_r0);

% correction to try to better match fermi velocity (turned off for now)
dft_based_scaling = 1.0;%1.0388;

Koshino_pi=@(dd) -2.7*dft_based_scaling*exp(-(dd-lattice_a/sqrt(3))/Koshino_r0);
Koshino_sigma=@(dd) 0.48*dft_based_scaling*exp(-(dd-layer_d(3))/Koshino_r0);


% Koshino_intra_t=@(vecvec) Koshino_pi(sqrt(vecvec(:,1).^2+vecvec(:,2).^2));
% Koshino_inter_t=@(vecvec) 0.0;
% Koshino_intra_t=@(vecrr) Koshino_pi(vecrr);


% define relaxation functional form

% first update atomic positions

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

top_index=[top_A_atoms,top_B_atoms];
bot_index=[bot_A_atoms,bot_B_atoms];

% scatter3(sc_all_points_shift(:,1),sc_all_points_shift(:,2),height_relax)

min(height_relax(top_index))-max(height_relax(bot_index));
max(height_relax(top_index))-min(height_relax(bot_index));

% max(height_relax(top_index))-min(height_relax(top_index))
% max(height_relax(bot_index))-min(height_relax(bot_index))







% inplaneGAB_table_bottom,inplaneGAB_vecs_bottom, 
% inplaneGAB_table_top,inplaneGAB_vecs_top,
% inplaneGAABB_table_bottom,inplaneGAABB_vecs_bottom,
% inplaneGAABB_table_top,inplaneGAABB_vecs_top,
%inplane_hops,inter_table,inter_vec,inter_hop);

sizetmp=size(inplaneGAB_table_bottom);

inplaneGAB_bottom_relax=zeros(sizetmp(1),6,sizetmp(3));
inplaneGAB_top_relax=zeros(sizetmp(1),6,sizetmp(3));


for indnn=1:sizetmp(3)
    del_tmp=inplaneGAB_vecs_bottom(indnn,:);
    pos_to_ind=squeeze(inplaneGAB_table_bottom(:,1,indnn));
    pos_from_ind=squeeze(inplaneGAB_table_bottom(:,2,indnn));
    
    %pos_to_tmp=sc_all_points_shift(pos_to_ind,1:2);
    %pos_from_tmp=sc_all_points_shift(pos_from_ind,1:2);
    
    pos_to_tmp_relax=Koshino_relax_shift(pos_to_ind,1:2);
    pos_from_tmp_relax=Koshino_relax_shift(pos_from_ind,1:2);
    
    pos_to_z_tmp=height_relax(pos_to_ind);
    pos_from_z_tmp=height_relax(pos_from_ind);
    
    % del_vec_relax=zeros(sizetmp(1),2);
    del_vec_relax=pos_to_tmp_relax-pos_from_tmp_relax;
    del_vec_relax(:,1)=del_vec_relax(:,1)+del_tmp(1);
    del_vec_relax(:,2)=del_vec_relax(:,2)+del_tmp(2);
    
    del_vec_relax_r=sqrt(del_vec_relax(:,1).^2+del_vec_relax(:,2).^2);
    
    %size(del_vec_relax_r)
    %hop_t_tmp=Koshino_intra_t(del_vec_relax_r);
    
    
    inplaneGAB_bottom_relax(:,1:2,indnn)=inplaneGAB_table_bottom(:,1:2,indnn);
    %inplaneGAB_bottom_relax(:,3,indnn)=hop_t_tmp(:);
    inplaneGAB_bottom_relax(:,4:5,indnn)=del_vec_relax(:,1:2)*lattice_a;
    inplaneGAB_bottom_relax(:,6,indnn)=pos_to_z_tmp-pos_from_z_tmp;
end



for indnn=1:sizetmp(3)
    del_tmp=inplaneGAB_vecs_top(indnn,:);
    pos_to_ind=squeeze(inplaneGAB_table_top(:,1,indnn));
    pos_from_ind=squeeze(inplaneGAB_table_top(:,2,indnn));
    
    %pos_to_tmp=sc_all_points_shift(pos_to_ind,1:2);
    %pos_from_tmp=sc_all_points_shift(pos_from_ind,1:2);
    
    pos_to_tmp_relax=Koshino_relax_shift(pos_to_ind,1:2);
    pos_from_tmp_relax=Koshino_relax_shift(pos_from_ind,1:2);
    
    pos_to_z_tmp=height_relax(pos_to_ind);
    pos_from_z_tmp=height_relax(pos_from_ind);
    
    % del_vec_relax=zeros(sizetmp(1),2);
    del_vec_relax=pos_to_tmp_relax-pos_from_tmp_relax;
    del_vec_relax(:,1)=del_vec_relax(:,1)+del_tmp(1);
    del_vec_relax(:,2)=del_vec_relax(:,2)+del_tmp(2);
    
    del_vec_relax_r=sqrt(del_vec_relax(:,1).^2+del_vec_relax(:,2).^2);
    
    %hop_t_tmp=Koshino_intra_t(del_vec_relax_r);
    
    
    inplaneGAB_top_relax(:,1:2,indnn)=inplaneGAB_table_top(:,1:2,indnn);
    %inplaneGAB_top_relax(:,3,indnn)=hop_t_tmp(:);
    inplaneGAB_top_relax(:,4:5,indnn)=del_vec_relax(:,1:2)*lattice_a;
    inplaneGAB_top_relax(:,6,indnn)=pos_to_z_tmp-pos_from_z_tmp;
end


NN_vec_table=zeros(tot_num,3,3);
NN_vec_table_index=ones(tot_num,1);
NN_pair_table=zeros(tot_num,1);

for ind11=1:sizetmp(1)
    tmpdd=inplaneGAB_bottom_relax(ind11,1:2,1);
    tmpvec=inplaneGAB_bottom_relax(ind11,4:6,1);
    tmpind1=tmpdd(1);
    tmpind2=tmpdd(2);
    
    tmpcnt1=NN_vec_table_index(tmpind1);
    tmpcnt2=NN_vec_table_index(tmpind2);
    NN_vec_table_index(tmpind1)=NN_vec_table_index(tmpind1)+1;
    NN_vec_table_index(tmpind2)=NN_vec_table_index(tmpind2)+1;
    
    NN_vec_table(tmpind1,tmpcnt1,:)=-tmpvec;
    NN_vec_table(tmpind2,tmpcnt2,:)=tmpvec;
    
    % Get nearest neighbor index (only need to do for first NN direction)
    NN_pair_table(tmpind1,1) = tmpind2;
    NN_pair_table(tmpind2,1) = tmpind1;
    
    tmpdd=inplaneGAB_bottom_relax(ind11,1:2,2);
    tmpvec=inplaneGAB_bottom_relax(ind11,4:6,2);
    tmpind1=tmpdd(1);
    tmpind2=tmpdd(2);
    
    tmpcnt1=NN_vec_table_index(tmpind1);
    tmpcnt2=NN_vec_table_index(tmpind2);
    NN_vec_table_index(tmpind1)=NN_vec_table_index(tmpind1)+1;
    NN_vec_table_index(tmpind2)=NN_vec_table_index(tmpind2)+1;
    
    NN_vec_table(tmpind1,tmpcnt1,:)=-tmpvec;
    NN_vec_table(tmpind2,tmpcnt2,:)=tmpvec;
    
    % Get nearest neighbor index (only need to do for first NN direction)
    NN_pair_table(tmpind1,2) = tmpind2;
    NN_pair_table(tmpind2,2) = tmpind1;
    
    tmpdd=inplaneGAB_bottom_relax(ind11,1:2,3);
    tmpvec=inplaneGAB_bottom_relax(ind11,4:6,3);
    tmpind1=tmpdd(1);
    tmpind2=tmpdd(2);
    
    tmpcnt1=NN_vec_table_index(tmpind1);
    tmpcnt2=NN_vec_table_index(tmpind2);
    NN_vec_table_index(tmpind1)=NN_vec_table_index(tmpind1)+1;
    NN_vec_table_index(tmpind2)=NN_vec_table_index(tmpind2)+1;
    
    NN_vec_table(tmpind1,tmpcnt1,:)=-tmpvec;
    NN_vec_table(tmpind2,tmpcnt2,:)=tmpvec;
    
    % Get nearest neighbor index (only need to do for first NN direction)
    NN_pair_table(tmpind1,3) = tmpind2;
    NN_pair_table(tmpind2,3) = tmpind1;
    
    
    
    % top
    tmpdd=inplaneGAB_top_relax(ind11,1:2,1);
    tmpvec=inplaneGAB_top_relax(ind11,4:6,1);
    tmpind1=tmpdd(1);
    tmpind2=tmpdd(2);
    
    tmpcnt1=NN_vec_table_index(tmpind1);
    tmpcnt2=NN_vec_table_index(tmpind2);
    NN_vec_table_index(tmpind1)=NN_vec_table_index(tmpind1)+1;
    NN_vec_table_index(tmpind2)=NN_vec_table_index(tmpind2)+1;
    
    NN_vec_table(tmpind1,tmpcnt1,:)=-tmpvec;
    NN_vec_table(tmpind2,tmpcnt2,:)=tmpvec;
    
    % Get nearest neighbor index (need all 3 NN directions)
    NN_pair_table(tmpind1,1) = tmpind2;
    NN_pair_table(tmpind2,1) = tmpind1;
    
    tmpdd=inplaneGAB_top_relax(ind11,1:2,2);
    tmpvec=inplaneGAB_top_relax(ind11,4:6,2);
    tmpind1=tmpdd(1);
    tmpind2=tmpdd(2);
    
    tmpcnt1=NN_vec_table_index(tmpind1);
    tmpcnt2=NN_vec_table_index(tmpind2);
    NN_vec_table_index(tmpind1)=NN_vec_table_index(tmpind1)+1;
    NN_vec_table_index(tmpind2)=NN_vec_table_index(tmpind2)+1;
    
    NN_vec_table(tmpind1,tmpcnt1,:)=-tmpvec;
    NN_vec_table(tmpind2,tmpcnt2,:)=tmpvec;
    
    % Get nearest neighbor index (need all 3 NN directions)
    NN_pair_table(tmpind1,2) = tmpind2;
    NN_pair_table(tmpind2,2) = tmpind1;
    
    tmpdd=inplaneGAB_top_relax(ind11,1:2,3);
    tmpvec=inplaneGAB_top_relax(ind11,4:6,3);
    tmpind1=tmpdd(1);
    tmpind2=tmpdd(2);
    
    tmpcnt1=NN_vec_table_index(tmpind1);
    tmpcnt2=NN_vec_table_index(tmpind2);
    NN_vec_table_index(tmpind1)=NN_vec_table_index(tmpind1)+1;
    NN_vec_table_index(tmpind2)=NN_vec_table_index(tmpind2)+1;
    
    NN_vec_table(tmpind1,tmpcnt1,:)=-tmpvec;
    NN_vec_table(tmpind2,tmpcnt2,:)=tmpvec;
    
    % Get nearest neighbor index (need all 3 NN directions)
    NN_pair_table(tmpind1,3) = tmpind2;
    NN_pair_table(tmpind2,3) = tmpind1;
    
end

% find normal vector
% convert to A unit

% NN_vec_table(:,:,1:2)=NN_vec_table(:,:,1:2)*lattice_a;

NN_norm_table=zeros(tot_num,3);

for indn=1:tot_num
    vec1=squeeze(NN_vec_table(indn,1,:));
    vec2=squeeze(NN_vec_table(indn,2,:));
    vec3=squeeze(NN_vec_table(indn,3,:));
    
    vecnorm1=cross(vec1,vec2);
    vecnorm2=cross(vec2,vec3);
    vecnorm3=cross(vec3,vec1);
    
    vecnorm1=vecnorm1*sign(vecnorm1(3));
    vecnorm2=vecnorm2*sign(vecnorm2(3));
    vecnorm3=vecnorm3*sign(vecnorm3(3));
    
    vecnorm1=vecnorm1/sqrt(dot(vecnorm1,vecnorm1));
    vecnorm2=vecnorm2/sqrt(dot(vecnorm2,vecnorm2));
    vecnorm3=vecnorm3/sqrt(dot(vecnorm3,vecnorm3));
    
    avg_norm_vec=(vecnorm1+vecnorm2+vecnorm3)/3;
    avg_norm_vec=avg_norm_vec/sqrt(dot(avg_norm_vec,avg_norm_vec));
   
    NN_norm_table(indn,:)=avg_norm_vec(:);
    
end

tmp_set_1=1:(tot_num/2);
tmp_set_2=(tot_num/2+1):(tot_num);
sel_set=tmp_set_2;
% 
% quiver(sc_all_points_shift_relax(sel_set,1),sc_all_points_shift_relax(sel_set,2),NN_norm_table(sel_set,1),NN_norm_table(sel_set,2))
% axis equal;

% inplaneGAB_top_relax and inplaneGAB_bottom_relax

% inplaneGAB_top/bottom_relax:
% (number of atoms) x 6 x (number of hopping dirs)

for indnn=1:sizetmp(3)
    ind_to=squeeze(inplaneGAB_top_relax(:,1,indnn));
    ind_from=squeeze(inplaneGAB_top_relax(:,2,indnn));
    
    hop_vec=squeeze(inplaneGAB_top_relax(:,4:6,indnn));
    hop_vec_r=sqrt(sum(hop_vec.^2,2));
    hop_vec_u=hop_vec;
    hop_vec_u(:,1)=hop_vec_u(:,1)./hop_vec_r;
    hop_vec_u(:,2)=hop_vec_u(:,2)./hop_vec_r;
    hop_vec_u(:,3)=hop_vec_u(:,3)./hop_vec_r;
    
    norm_vec_from=NN_norm_table(ind_from,:);
    norm_vec_to=NN_norm_table(ind_to,:);
    
    tmp_sigma=Koshino_sigma(hop_vec_r);
    tmp_pi=Koshino_pi(hop_vec_r);
    
    clear tmpfa2;
    clear tmpfa3;
    
    tmpfa1=sum(hop_vec_u.*norm_vec_from,2).*sum(hop_vec_u.*norm_vec_to,2);
    
    tmpfa2(:,1)=norm_vec_from(:,2).*hop_vec_u(:,3)-norm_vec_from(:,3).*hop_vec_u(:,2);
    tmpfa2(:,2)=norm_vec_from(:,3).*hop_vec_u(:,1)-norm_vec_from(:,1).*hop_vec_u(:,3);
    tmpfa2(:,3)=norm_vec_from(:,1).*hop_vec_u(:,2)-norm_vec_from(:,2).*hop_vec_u(:,1);
    
    tmpfa3(:,1)=norm_vec_to(:,2).*hop_vec_u(:,3)-norm_vec_to(:,3).*hop_vec_u(:,2);
    tmpfa3(:,2)=norm_vec_to(:,3).*hop_vec_u(:,1)-norm_vec_to(:,1).*hop_vec_u(:,3);
    tmpfa3(:,3)=norm_vec_to(:,1).*hop_vec_u(:,2)-norm_vec_to(:,2).*hop_vec_u(:,1);
    
    tmpfa4=sum(tmpfa2.*tmpfa3,2);
    
    hop_tmp1=tmp_sigma.*tmpfa1;
    hop_tmp2=tmp_pi.*tmpfa4;
    
    inplaneGAB_top_relax(:,3,indnn)=hop_tmp1+hop_tmp2;
    
    
    
    
    ind_to=squeeze(inplaneGAB_bottom_relax(:,1,indnn));
    ind_from=squeeze(inplaneGAB_bottom_relax(:,2,indnn));
    
    hop_vec=squeeze(inplaneGAB_bottom_relax(:,4:6,indnn));
    hop_vec_r=sqrt(sum(hop_vec.^2,2));
    hop_vec_u=hop_vec;
    hop_vec_u(:,1)=hop_vec_u(:,1)./hop_vec_r;
    hop_vec_u(:,2)=hop_vec_u(:,2)./hop_vec_r;
    hop_vec_u(:,3)=hop_vec_u(:,3)./hop_vec_r;
    
    norm_vec_from=NN_norm_table(ind_from,:);
    norm_vec_to=NN_norm_table(ind_to,:);
    
    tmp_sigma=Koshino_sigma(hop_vec_r);
    tmp_pi=Koshino_pi(hop_vec_r);
    
    clear tmpfa2;
    clear tmpfa3;
    
    tmpfa1=sum(hop_vec_u.*norm_vec_from,2).*sum(hop_vec_u.*norm_vec_to,2);
    
    tmpfa2(:,1)=norm_vec_from(:,2).*hop_vec_u(:,3)-norm_vec_from(:,3).*hop_vec_u(:,2);
    tmpfa2(:,2)=norm_vec_from(:,3).*hop_vec_u(:,1)-norm_vec_from(:,1).*hop_vec_u(:,3);
    tmpfa2(:,3)=norm_vec_from(:,1).*hop_vec_u(:,2)-norm_vec_from(:,2).*hop_vec_u(:,1);
    
    tmpfa3(:,1)=norm_vec_to(:,2).*hop_vec_u(:,3)-norm_vec_to(:,3).*hop_vec_u(:,2);
    tmpfa3(:,2)=norm_vec_to(:,3).*hop_vec_u(:,1)-norm_vec_to(:,1).*hop_vec_u(:,3);
    tmpfa3(:,3)=norm_vec_to(:,1).*hop_vec_u(:,2)-norm_vec_to(:,2).*hop_vec_u(:,1);
    
    tmpfa4=sum(tmpfa2.*tmpfa3,2);
    
    hop_tmp1=tmp_sigma.*tmpfa1;
    hop_tmp2=tmp_pi.*tmpfa4;
    
    inplaneGAB_bottom_relax(:,3,indnn)=hop_tmp1+hop_tmp2;
    
    
end




sizetmp2=size(inplaneGAABB_table_bottom);

inplaneGAABB_bottom_relax=zeros(sizetmp2(1),5,sizetmp2(3));
inplaneGAABB_top_relax=zeros(sizetmp2(1),5,sizetmp2(3));


for indnn=1:sizetmp2(3)
    del_tmp=inplaneGAABB_vecs_bottom(indnn,:);
    pos_to_ind=squeeze(inplaneGAABB_table_bottom(:,1,indnn));
    pos_from_ind=squeeze(inplaneGAABB_table_bottom(:,2,indnn));
    
    %pos_to_tmp=sc_all_points_shift(pos_to_ind,1:2);
    %pos_from_tmp=sc_all_points_shift(pos_from_ind,1:2);
    
    pos_to_tmp_relax=Koshino_relax_shift(pos_to_ind,1:2);
    pos_from_tmp_relax=Koshino_relax_shift(pos_from_ind,1:2);
    
    pos_to_z_tmp=height_relax(pos_to_ind);
    pos_from_z_tmp=height_relax(pos_from_ind);
    
    % del_vec_relax=zeros(sizetmp(1),2);
    del_vec_relax=pos_to_tmp_relax-pos_from_tmp_relax;
    del_vec_relax(:,1)=del_vec_relax(:,1)+del_tmp(1);
    del_vec_relax(:,2)=del_vec_relax(:,2)+del_tmp(2);
    
    del_vec_relax_r=sqrt(del_vec_relax(:,1).^2+del_vec_relax(:,2).^2);
    
    %hop_t_tmp=Koshino_intra_t(del_vec_relax_r);
    
    
    inplaneGAABB_bottom_relax(:,1:2,indnn)=inplaneGAABB_table_bottom(:,1:2,indnn);
    inplaneGAABB_bottom_relax(:,4:5,indnn)=del_vec_relax(:,1:2)*lattice_a;
    inplaneGAABB_bottom_relax(:,6,indnn)=pos_to_z_tmp-pos_from_z_tmp;
    
    
    
    ind_to=squeeze(inplaneGAABB_bottom_relax(:,1,indnn));
    ind_from=squeeze(inplaneGAABB_bottom_relax(:,2,indnn));
    
    hop_vec=squeeze(inplaneGAABB_bottom_relax(:,4:6,indnn));
    hop_vec_r=sqrt(sum(hop_vec.^2,2));
    
    hop_vec_u=hop_vec;
    hop_vec_u(:,1)=hop_vec_u(:,1)./hop_vec_r;
    hop_vec_u(:,2)=hop_vec_u(:,2)./hop_vec_r;
    hop_vec_u(:,3)=hop_vec_u(:,3)./hop_vec_r;
    
    norm_vec_from=NN_norm_table(ind_from,:);
    norm_vec_to=NN_norm_table(ind_to,:);
    
    tmp_sigma=Koshino_sigma(hop_vec_r);
    tmp_pi=Koshino_pi(hop_vec_r);
    
    clear tmpfa2;
    clear tmpfa3;
    
    tmpfa1=sum(hop_vec_u.*norm_vec_from,2).*sum(hop_vec_u.*norm_vec_to,2);
    
    tmpfa2(:,1)=norm_vec_from(:,2).*hop_vec_u(:,3)-norm_vec_from(:,3).*hop_vec_u(:,2);
    tmpfa2(:,2)=norm_vec_from(:,3).*hop_vec_u(:,1)-norm_vec_from(:,1).*hop_vec_u(:,3);
    tmpfa2(:,3)=norm_vec_from(:,1).*hop_vec_u(:,2)-norm_vec_from(:,2).*hop_vec_u(:,1);
    
    tmpfa3(:,1)=norm_vec_to(:,2).*hop_vec_u(:,3)-norm_vec_to(:,3).*hop_vec_u(:,2);
    tmpfa3(:,2)=norm_vec_to(:,3).*hop_vec_u(:,1)-norm_vec_to(:,1).*hop_vec_u(:,3);
    tmpfa3(:,3)=norm_vec_to(:,1).*hop_vec_u(:,2)-norm_vec_to(:,2).*hop_vec_u(:,1);
    
    tmpfa4=sum(tmpfa2.*tmpfa3,2);
    
    hop_tmp1=tmp_sigma.*tmpfa1;
    hop_tmp2=tmp_pi.*tmpfa4;
    
    inplaneGAABB_bottom_relax(:,3,indnn)=hop_tmp1+hop_tmp2;
    
    
end

for indnn=1:sizetmp2(3)
    del_tmp=inplaneGAABB_vecs_top(indnn,:);
    pos_to_ind=squeeze(inplaneGAABB_table_top(:,1,indnn));
    pos_from_ind=squeeze(inplaneGAABB_table_top(:,2,indnn));
    
    %pos_to_tmp=sc_all_points_shift(pos_to_ind,1:2);
    %pos_from_tmp=sc_all_points_shift(pos_from_ind,1:2);
    
    pos_to_tmp_relax=Koshino_relax_shift(pos_to_ind,1:2);
    pos_from_tmp_relax=Koshino_relax_shift(pos_from_ind,1:2);
    
    pos_to_z_tmp=height_relax(pos_to_ind);
    pos_from_z_tmp=height_relax(pos_from_ind);
    
    % del_vec_relax=zeros(sizetmp(1),2);
    del_vec_relax=pos_to_tmp_relax-pos_from_tmp_relax;
    del_vec_relax(:,1)=del_vec_relax(:,1)+del_tmp(1);
    del_vec_relax(:,2)=del_vec_relax(:,2)+del_tmp(2);
    
    del_vec_relax_r=sqrt(del_vec_relax(:,1).^2+del_vec_relax(:,2).^2);
    
    %hop_t_tmp=Koshino_intra_t(del_vec_relax_r);
    
    
    inplaneGAABB_top_relax(:,1:2,indnn)=inplaneGAABB_table_top(:,1:2,indnn);
    % inplaneGAABB_top_relax(:,3,indnn)=hop_t_tmp(:);
    inplaneGAABB_top_relax(:,4:5,indnn)=del_vec_relax(:,1:2)*lattice_a;
    inplaneGAABB_top_relax(:,6,indnn)=pos_to_z_tmp-pos_from_z_tmp;
    
    
    ind_to=squeeze(inplaneGAABB_top_relax(:,1,indnn));
    ind_from=squeeze(inplaneGAABB_top_relax(:,2,indnn));
    
    hop_vec=squeeze(inplaneGAABB_top_relax(:,4:6,indnn));
    hop_vec_r=sqrt(sum(hop_vec.^2,2));
    
    hop_vec_u=hop_vec;
    hop_vec_u(:,1)=hop_vec_u(:,1)./hop_vec_r;
    hop_vec_u(:,2)=hop_vec_u(:,2)./hop_vec_r;
    hop_vec_u(:,3)=hop_vec_u(:,3)./hop_vec_r;
    
    norm_vec_from=NN_norm_table(ind_from,:);
    norm_vec_to=NN_norm_table(ind_to,:);
    
    tmp_sigma=Koshino_sigma(hop_vec_r);
    tmp_pi=Koshino_pi(hop_vec_r);
    
    tmpfa1=sum(hop_vec_u.*norm_vec_from,2).*sum(hop_vec_u.*norm_vec_to,2);
    
    clear tmpfa2;
    clear tmpfa3;
    
    tmpfa2(:,1)=norm_vec_from(:,2).*hop_vec_u(:,3)-norm_vec_from(:,3).*hop_vec_u(:,2);
    tmpfa2(:,2)=norm_vec_from(:,3).*hop_vec_u(:,1)-norm_vec_from(:,1).*hop_vec_u(:,3);
    tmpfa2(:,3)=norm_vec_from(:,1).*hop_vec_u(:,2)-norm_vec_from(:,2).*hop_vec_u(:,1);
    
    tmpfa3(:,1)=norm_vec_to(:,2).*hop_vec_u(:,3)-norm_vec_to(:,3).*hop_vec_u(:,2);
    tmpfa3(:,2)=norm_vec_to(:,3).*hop_vec_u(:,1)-norm_vec_to(:,1).*hop_vec_u(:,3);
    tmpfa3(:,3)=norm_vec_to(:,1).*hop_vec_u(:,2)-norm_vec_to(:,2).*hop_vec_u(:,1);
    
    tmpfa4=sum(tmpfa2.*tmpfa3,2);
    
    hop_tmp1=tmp_sigma.*tmpfa1;
    hop_tmp2=tmp_pi.*tmpfa4;
    
    inplaneGAABB_top_relax(:,3,indnn)=hop_tmp1+hop_tmp2;
    
    
end




% interlayer terms
% inter_vec, inter_table

sizetmp=size(inter_table);

% inter_relax: (number of interlayer terms) x 6
%             [1]:      index TO
%             [2]:      index FROM
%             [3]:      t (hopping term)
%             [4:6]:    3D hopping vector
inter_relax=zeros(sizetmp(1),6);
inter_relax(:,1:2)=inter_table(:,1:2);

pos_to_ind=squeeze(inter_table(:,1));
pos_from_ind=squeeze(inter_table(:,2));
pos_to_tmp=sc_all_points_shift(pos_to_ind,1:2);
pos_from_tmp=sc_all_points_shift(pos_from_ind,1:2);

pos_to_z_tmp=height_relax(pos_to_ind);
pos_from_z_tmp=height_relax(pos_from_ind);
pos_hop_zvec=pos_to_z_tmp-pos_from_z_tmp;

norm_vec_from=NN_norm_table(pos_from_ind,:);
norm_vec_to=NN_norm_table(pos_to_ind,:);


pos_to_tmp_relax=Koshino_relax_shift(pos_to_ind,1:2);
pos_from_tmp_relax=Koshino_relax_shift(pos_from_ind,1:2);
del_vec_relax=pos_to_tmp_relax-pos_from_tmp_relax;

inter_vec_relax=inter_vec;
inter_vec_relax=inter_vec_relax+del_vec_relax;

clear inter_vec_relax_3d;
inter_vec_relax_3d(:,1:2)=inter_vec_relax(:,1:2)*lattice_a;
inter_vec_relax_3d(:,3)=pos_hop_zvec(:);
inter_relax(:,4:6)=inter_vec_relax_3d(:,1:3);

%del_hh=height_relax(pos_to_ind)-height_relax(pos_from_ind);

inter_vec_relax_dd=sqrt(sum(inter_vec_relax_3d.^2,2));

hop_vec_u=inter_vec_relax_3d;
hop_vec_u(:,1)=hop_vec_u(:,1)./inter_vec_relax_dd;
hop_vec_u(:,2)=hop_vec_u(:,2)./inter_vec_relax_dd;
hop_vec_u(:,3)=hop_vec_u(:,3)./inter_vec_relax_dd;

% inter_vec_relax_dd=sqrt(inter_vec_relax(:,1).^2+inter_vec_relax(:,2).^2+(del_hh/lattice_a).^2);
inter_vec_relax_pi=Koshino_pi(inter_vec_relax_dd);
inter_vec_relax_sigma=Koshino_sigma(inter_vec_relax_dd);


tmpfa1=sum(hop_vec_u.*norm_vec_from,2).*sum(hop_vec_u.*norm_vec_to,2);


clear tmpfa2;
clear tmpfa3;
    

% tmpfa2 = norm_vec_from (x) hop_vec_u
tmpfa2(:,1)=norm_vec_from(:,2).*hop_vec_u(:,3)-norm_vec_from(:,3).*hop_vec_u(:,2);
tmpfa2(:,2)=norm_vec_from(:,3).*hop_vec_u(:,1)-norm_vec_from(:,1).*hop_vec_u(:,3);
tmpfa2(:,3)=norm_vec_from(:,1).*hop_vec_u(:,2)-norm_vec_from(:,2).*hop_vec_u(:,1);

% tmpfa2 = norm_vec_to (x) hop_vec_u
tmpfa3(:,1)=norm_vec_to(:,2).*hop_vec_u(:,3)-norm_vec_to(:,3).*hop_vec_u(:,2);
tmpfa3(:,2)=norm_vec_to(:,3).*hop_vec_u(:,1)-norm_vec_to(:,1).*hop_vec_u(:,3);
tmpfa3(:,3)=norm_vec_to(:,1).*hop_vec_u(:,2)-norm_vec_to(:,2).*hop_vec_u(:,1);

% tmpfa4 = ?
tmpfa4=sum(tmpfa2.*tmpfa3,2);

hop_tmp1=inter_vec_relax_sigma.*tmpfa1;
hop_tmp2=inter_vec_relax_pi.*tmpfa4;


% interlayer terms
if (inter_type == 0) % Koshino model
    inter_relax(:,3)=hop_tmp1+hop_tmp2;
elseif (inter_type == 1) % DFT model
    
    
    inplaneGAB_top_relax(:,3,:) = 0;
    inplaneGAABB_top_relax(:,3,:) = 0;
    
    inplaneGAB_bottom_relax(:,3,:) = 0;
    inplaneGAABB_bottom_relax(:,3,:) = 0;
    
    % graphene TBH from DFT
    %{
    t1 = -2.8922;
    t2 =  0.2425;
    t3 = -0.2656;
    t4 =  0.0235;
    t5 =  0.0524;
    t6 = -0.0209;
    t7 = -0.0148;
    t8 = -0.0211;
    
    inplaneGAB_top_relax(:,3,1:3)       = t1;
    inplaneGAABB_top_relax(:,3,1:6)     = t2;
    inplaneGAB_top_relax(:,3,4:6)       = t3;
    inplaneGAB_top_relax(:,3,7:12)      = t4;
    inplaneGAABB_top_relax(:,3,7:12)    = t5;
    inplaneGAABB_top_relax(:,3,13:18)   = t6;
    inplaneGAB_top_relax(:,3,13:18)     = t7;
    inplaneGAB_top_relax(:,3,19:21)     = t8;

    inplaneGAB_bottom_relax(:,3,1:3)       = t1;
    inplaneGAABB_bottom_relax(:,3,1:6)     = t2;
    inplaneGAB_bottom_relax(:,3,4:6)       = t3;
    inplaneGAB_bottom_relax(:,3,7:12)      = t4;
    inplaneGAABB_bottom_relax(:,3,7:12)    = t5;
    inplaneGAABB_bottom_relax(:,3,13:18)   = t6;
    inplaneGAB_bottom_relax(:,3,13:18)     = t7;
    inplaneGAB_bottom_relax(:,3,19:21)     = t8;
    %}
    
    % graphene TBH from DFT with strain corrections
    alph0 = -4.878; % currently unused
    t1 = -2.822;
    alph1 = 4.007;
    t2 =  0.254;
    alph2 = -0.463;
    t3 = -0.180;
    alph3 = 0.624;
    
    % turns off iso strain
    %alph1 = 0;
    %alph2 = 0;
    %alph3 = 0;
    
    avg_nn_bond_lengths_AB_top = sqrt(  sum(inplaneGAB_top_relax(:,4:6,1:3).^2  , 2));
    iso_strain_AB_top = (avg_nn_bond_lengths_AB_top - lattice_a/sqrt(3)) / (lattice_a/sqrt(3));

    avg_nn_bond_lengths_AB2_top = sqrt(  sum(inplaneGAB_top_relax(:,4:6,4:6).^2  , 2));
    iso_strain_AB2_top = (avg_nn_bond_lengths_AB2_top - 2*lattice_a/sqrt(3)) / (2*lattice_a/sqrt(3));
    
    %avg_nn_bond_lengths_AABB_top = mean(sqrt( inplaneGAABB_top_relax(:,4,1:6).^2 + inplaneGAABB_top_relax(:,5,1:6).^2) ,3);
    avg_nn_bond_lengths_AABB_top =  sqrt(  sum(inplaneGAABB_top_relax(:,4:6,1:6).^2  , 2));
    iso_strain_AABB_top = (avg_nn_bond_lengths_AABB_top - lattice_a) / (lattice_a);
    
    avg_nn_bond_lengths_AB_bot = sqrt(  sum(inplaneGAB_bottom_relax(:,4:6,1:3).^2  , 2));
    iso_strain_AB_bot = (avg_nn_bond_lengths_AB_bot - lattice_a/sqrt(3)) / (lattice_a/sqrt(3));
    
    avg_nn_bond_lengths_AB2_bot = sqrt(  sum(inplaneGAB_bottom_relax(:,4:6,4:6).^2  , 2));
    iso_strain_AB2_bot = (avg_nn_bond_lengths_AB2_bot - 2*lattice_a/sqrt(3)) / (2*lattice_a/sqrt(3));
        
    avg_nn_bond_lengths_AABB_bot =  sqrt(  sum(inplaneGAABB_bottom_relax(:,4:6,1:6).^2  , 2));
    iso_strain_AABB_bot = (avg_nn_bond_lengths_AABB_bot - lattice_a) / (lattice_a);

    
    inplaneGAB_top_relax(:,3,1:3)       =  t1 + alph1*iso_strain_AB_top;
    inplaneGAABB_top_relax(:,3,1:6)     =  t2 + alph2*iso_strain_AABB_top; %mtimes(t2 + alph2*iso_strain_AABB_top, [1 1 1 1 1 1]);
    inplaneGAB_top_relax(:,3,4:6)       =  t3 + alph3*iso_strain_AB2_top;
    
    inplaneGAB_bottom_relax(:,3,1:3)       = t1 + alph1*iso_strain_AB_bot;
    inplaneGAABB_bottom_relax(:,3,1:6)     = t2 + alph2*iso_strain_AABB_bot;
    inplaneGAB_bottom_relax(:,3,4:6)       = t3 + alph3*iso_strain_AB2_bot;

    % Get the list of all indices
    from_indices = inter_relax(:,1);
    to_indices = inter_relax(:,2);
    
    theta1 = zeros(length(inter_relax),1);
    theta2 = theta1;
    
    for nn_dir = 1:3
        
        % NN_pair_table is not used anymore, can remove from lines [168 to 274] of this file
        %from_bond_pairs = NN_pair_table(from_indices,nn_dir); % index of NN bonding direction
        %to_bond_pairs = NN_pair_table(to_indices,nn_dir); % index of NN bonding direction
        
        % nearest-neighbor lookup table:
        % first index tells which row of table to use
        % second index tells sign of bonding vector
        
        nn_lookup_table = zeros(tot_num,2);
        for table_idx = 1:size(inplaneGAB_bottom_relax(:,:,nn_dir),1)
            tmp_idx_1 = inplaneGAB_bottom_relax(table_idx,1,nn_dir);
            tmp_idx_2 = inplaneGAB_bottom_relax(table_idx,2,nn_dir);
            nn_lookup_table(tmp_idx_1,:) = [table_idx, -1];
            nn_lookup_table(tmp_idx_2,:) = [table_idx, 1];
        end
    
        for table_idx = 1:size(inplaneGAB_bottom_relax(:,:,nn_dir),1)
            tmp_idx_1 = inplaneGAB_top_relax(table_idx,1,nn_dir);
            tmp_idx_2 = inplaneGAB_top_relax(table_idx,2,nn_dir);
            nn_lookup_table(tmp_idx_1,:) = [table_idx, -1];
            nn_lookup_table(tmp_idx_2,:) = [table_idx, 1];
        end
        
        %from_bond_dirs = sc_all_points_shift(from_bond_pairs,:) - sc_all_points_shift(from_indices,:);
        %to_bond_dirs = sc_all_points_shift(to_bond_pairs,:) - sc_all_points_shift(to_indices,:);
                
        from_bond_dirs  = inplaneGAB_top_relax(      nn_lookup_table(from_indices,  1),     4:5, nn_dir)  .*    nn_lookup_table(from_indices,   2);
        to_bond_dirs    = inplaneGAB_bottom_relax(   nn_lookup_table(to_indices,    1),     4:5, nn_dir)  .*    nn_lookup_table(to_indices,     2);
 
        theta1_here = angle(from_bond_dirs(:,1) + 1i*from_bond_dirs(:,2));
        theta2_here = angle(to_bond_dirs(:,1) + 1i*to_bond_dirs(:,2));
        
        %fprintf("nn = %d, theta1_here(2) = %f \n",nn_dir, theta1_here(2)*180/pi);
        %fprintf("nn = %d, theta2_here(2) = %f \n",nn_dir, theta2_here(2)*180/pi);
        
        theta1 = theta1 + theta1_here;
        theta2 = theta2 + theta2_here;
        
    end
    
    theta1 = theta1/3;
    theta2 = theta2/3;
    
    % inter_vec_relax_3d points from TO to FROM so reverse it
    dft_inter_vec = -inter_vec_relax_3d;

    inter_relax(:,3) = dft_interlayer_coupling(dft_inter_vec,theta1,theta2,lattice_a);
    %inter_relax(:,3) = 0;
end












