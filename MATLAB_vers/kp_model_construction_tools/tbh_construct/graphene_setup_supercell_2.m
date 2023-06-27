% this is the file to establish atomic sites and the link between them
% both intra layer and inter layer

% input: sc_m, sc_n

code_check=false;


% expect sc_m >=2    sc_n>=1

% graphene_setup_supercell_env


% rescale now?

% inter_cutoff=3.99;

inter_cutoff=2.2975;



% inplane_hops_wan_mono_N8_5band=[-2.8217, 0.2540, -0.1797,-0.0143,0.0324,-0.0037,-0.0109,-0.0078];

inplane_hops_wan_mono_N8_pz=[-2.8922  ,  0.2425  , -0.2656  ,  0.0235, 0.0524 ,  -0.0209,   -0.0148  , -0.0211,0.3504];


inplane_hops=inplane_hops_wan_mono_N8_pz;




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





% build connections

% bottom layer---------------------------------
% actual boundary:     -M+1  to  N-1     (N+M-1)
% actual boundary:     0     to  N+2M-1  (N+2M)

% 6 directions: a1,a2,a2-a1,-a1,-a2,a1-a2
actual_index_mat_b=zeros(sc_m+sc_n-1,sc_n+2*sc_m,7);

for ind1=(-sc_m+1):(sc_n-1)
    for ind2=0:(sc_n+2*sc_m-1)
        
        tmp=index_mat_b(ind1-(-sc_m+1)+2,ind2+2);
        
        if tmp~=0
            actual_index_mat_b(ind1-(-sc_m+1)+1,ind2+1,1)=tmp;
            
            
            % [1,0]
            tmp1=index_mat_b(ind1+1-(-sc_m+1)+2,ind2+2);
            if tmp1~=0
                
                actual_index_mat_b(ind1-(-sc_m+1)+1,ind2+1,2)=tmp1;
            else
                
                [ newvec ] = moveto_primitive_method2( [sc_n,sc_m],[-sc_m,sc_n+sc_m],[ind1+1,ind2]' );
                actual_index_mat_b(ind1-(-sc_m+1)+1,ind2+1,2)=index_mat_b(newvec(1)-(-sc_m+1)+2,newvec(2)+2);
                
            end
            
            
            % [0,1]
            tmp1=index_mat_b(ind1-(-sc_m+1)+2,ind2+1+2);
            if tmp1~=0
                actual_index_mat_b(ind1-(-sc_m+1)+1,ind2+1,3)=tmp1;
            else
                [ newvec ] = moveto_primitive_method2( [sc_n,sc_m],[-sc_m,sc_n+sc_m],[ind1,ind2+1]' );
                actual_index_mat_b(ind1-(-sc_m+1)+1,ind2+1,3)=index_mat_b(newvec(1)-(-sc_m+1)+2,newvec(2)+2);
                
            end
            
            % [-1,1]
            tmp1=index_mat_b(ind1-1-(-sc_m+1)+2,ind2+1+2);
            if tmp1~=0
                actual_index_mat_b(ind1-(-sc_m+1)+1,ind2+1,4)=tmp1;
            else
                [ newvec ] = moveto_primitive_method2( [sc_n,sc_m],[-sc_m,sc_n+sc_m],[ind1-1,ind2+1]' );
                actual_index_mat_b(ind1-(-sc_m+1)+1,ind2+1,4)=index_mat_b(newvec(1)-(-sc_m+1)+2,newvec(2)+2);
                
            end
            
            % [-1,0]
            tmp1=index_mat_b(ind1-1-(-sc_m+1)+2,ind2+2);
            if tmp1~=0
                actual_index_mat_b(ind1-(-sc_m+1)+1,ind2+1,5)=tmp1;
            else
                [ newvec ] = moveto_primitive_method2( [sc_n,sc_m],[-sc_m,sc_n+sc_m],[ind1-1,ind2]' );
                actual_index_mat_b(ind1-(-sc_m+1)+1,ind2+1,5)=index_mat_b(newvec(1)-(-sc_m+1)+2,newvec(2)+2);
                
            end
            
            
            % [0,-1]
            tmp1=index_mat_b(ind1-(-sc_m+1)+2,ind2-1+2);
            if tmp1~=0
                actual_index_mat_b(ind1-(-sc_m+1)+1,ind2+1,6)=tmp1;
            else
                [ newvec ] = moveto_primitive_method2( [sc_n,sc_m],[-sc_m,sc_n+sc_m],[ind1,ind2-1]' );
                actual_index_mat_b(ind1-(-sc_m+1)+1,ind2+1,6)=index_mat_b(newvec(1)-(-sc_m+1)+2,newvec(2)+2);
                
            end
            
            % [1,-1]
            tmp1=index_mat_b(ind1+1-(-sc_m+1)+2,ind2-1+2);
            if tmp1~=0
                actual_index_mat_b(ind1-(-sc_m+1)+1,ind2+1,7)=tmp1;
            else
                [ newvec ] = moveto_primitive_method2( [sc_n,sc_m],[-sc_m,sc_n+sc_m],[ind1+1,ind2-1]' );
                actual_index_mat_b(ind1-(-sc_m+1)+1,ind2+1,7)=index_mat_b(newvec(1)-(-sc_m+1)+2,newvec(2)+2);
                
            end
            
            
        end
        
        
    end
end


if code_check
    for ind1=(-sc_m+1):(sc_n-1)
        for ind2=0:(sc_n+2*sc_m-1)
            
            tmp=actual_index_mat_b(ind1-(-sc_m+1)+1,ind2+1,1);
            
            if tmp~=0
                arr1=actual_index_mat_b(ind1-(-sc_m+1)+1,ind2+1,2:7);
                tmp2=sum(arr1==0);
                
                if tmp2>0
                    'warning bottom wrong!'
                end
            end
        end
    end
end



% top layer---------------------------------
% actual boundary:     -N+1  to  M-1     (N+M-1)
% actual boundary:     0     to  2N+M-1  (2N+M)

% 6 directions: a1,a2,a2-a1,-a1,-a2,a1-a2


actual_index_mat_t=zeros(sc_m+sc_n-1,2*sc_n+sc_m,7);


for ind1=(-sc_n+1):(sc_m-1)
    for ind2=0:(2*sc_n+sc_m-1)
        
        tmp=index_mat_t(ind1-(-sc_n+1)+2,ind2+2);
        
        if tmp~=0
            actual_index_mat_t(ind1-(-sc_n+1)+1,ind2+1,1)=tmp;
            
            % [1,0]
            tmp1=index_mat_t(ind1+1-(-sc_n+1)+2,ind2+2);
            if tmp1~=0
                
                actual_index_mat_t(ind1-(-sc_n+1)+1,ind2+1,2)=tmp1;
            else
                
                [ newvec ] = moveto_primitive_method2( [sc_m,sc_n],[-sc_n,sc_n+sc_m],[ind1+1,ind2]' );
                actual_index_mat_t(ind1-(-sc_n+1)+1,ind2+1,2)=index_mat_t(newvec(1)-(-sc_n+1)+2,newvec(2)+2);
                
            end
            
            
            % [0,1]
            tmp1=index_mat_t(ind1-(-sc_n+1)+2,ind2+1+2);
            if tmp1~=0
                actual_index_mat_t(ind1-(-sc_n+1)+1,ind2+1,3)=tmp1;
            else
                [ newvec ] = moveto_primitive_method2( [sc_m,sc_n],[-sc_n,sc_n+sc_m],[ind1,ind2+1]' );
                actual_index_mat_t(ind1-(-sc_n+1)+1,ind2+1,3)=index_mat_t(newvec(1)-(-sc_n+1)+2,newvec(2)+2);
                
            end
            
            % [-1,1]
            tmp1=index_mat_t(ind1-1-(-sc_n+1)+2,ind2+1+2);
            if tmp1~=0
                actual_index_mat_t(ind1-(-sc_n+1)+1,ind2+1,4)=tmp1;
            else
                [ newvec ] = moveto_primitive_method2([sc_m,sc_n],[-sc_n,sc_n+sc_m],[ind1-1,ind2+1]' );
                actual_index_mat_t(ind1-(-sc_n+1)+1,ind2+1,4)=index_mat_t(newvec(1)-(-sc_n+1)+2,newvec(2)+2);
                
            end
            
            % [-1,0]
            tmp1=index_mat_t(ind1-1-(-sc_n+1)+2,ind2+2);
            if tmp1~=0
                actual_index_mat_t(ind1-(-sc_n+1)+1,ind2+1,5)=tmp1;
            else
                [ newvec ] = moveto_primitive_method2([sc_m,sc_n],[-sc_n,sc_n+sc_m],[ind1-1,ind2]' );
                actual_index_mat_t(ind1-(-sc_n+1)+1,ind2+1,5)=index_mat_t(newvec(1)-(-sc_n+1)+2,newvec(2)+2);
                
            end
            
            
            % [0,-1]
            tmp1=index_mat_t(ind1-(-sc_n+1)+2,ind2-1+2);
            if tmp1~=0
                actual_index_mat_t(ind1-(-sc_n+1)+1,ind2+1,6)=tmp1;
            else
                [ newvec ] = moveto_primitive_method2([sc_m,sc_n],[-sc_n,sc_n+sc_m],[ind1,ind2-1]' );
                actual_index_mat_t(ind1-(-sc_n+1)+1,ind2+1,6)=index_mat_t(newvec(1)-(-sc_n+1)+2,newvec(2)+2);
                
            end
            
            % [1,-1]
            tmp1=index_mat_t(ind1+1-(-sc_n+1)+2,ind2-1+2);
            if tmp1~=0
                actual_index_mat_t(ind1-(-sc_n+1)+1,ind2+1,7)=tmp1;
            else
                [ newvec ] = moveto_primitive_method2([sc_m,sc_n],[-sc_n,sc_n+sc_m],[ind1+1,ind2-1]' );
                actual_index_mat_t(ind1-(-sc_n+1)+1,ind2+1,7)=index_mat_t(newvec(1)-(-sc_n+1)+2,newvec(2)+2);
                
            end
            
            
        end
        
        
    end
end

if code_check
    for ind1=(-sc_n+1):(sc_m-1)
        for ind2=0:(2*sc_n+sc_m-1)
            
            tmp=actual_index_mat_t(ind1-(-sc_n+1)+1,ind2+1,1);
            
            if tmp~=0
                arr1=actual_index_mat_t(ind1-(-sc_n+1)+1,ind2+1,2:7);
                tmp2=sum(arr1==0);
                
                if tmp2>0
                    'warning top wrong!'
                end
            end
        end
    end
end


% inverse lookup table

% bottom layer---------------------------------
% actual boundary:     -M+1  to  N-1     (N+M-1)
% actual boundary:     0     to  N+2M-1  (N+2M)

% inverse_index_mat_b(:,1)=inverse_index_mat_b(:,1)-(-sc_m+1)+1;
% inverse_index_mat_b(:,2)=inverse_index_mat_b(:,2)+1;

access_neighbor_bottom = @(ii,jj,kk) actual_index_mat_b(ii-(-sc_m+1)+1,jj+1,kk);


% top layer---------------------------------
% actual boundary:     -N+1  to  M-1     (N+M-1)
% actual boundary:     0     to  2N+M-1  (2N+M)

% inverse_index_mat_t(:,1)=inverse_index_mat_t(:,1)-(-sc_n+1)+1;
% inverse_index_mat_t(:,2)=inverse_index_mat_t(:,2)+1;

access_neighbor_top = @(ii,jj,kk) actual_index_mat_t(ii-(-sc_n+1)+1,jj+1,kk);




tot_num=4*num_pos_t;

index_bot_atoms=1:(2*num_pos_b);
bot_A_atoms=1:2:(2*num_pos_b);
bot_B_atoms=2:2:(2*num_pos_b);
index_top_atoms=(2*num_pos_b+1):(4*num_pos_b);
top_A_atoms=(2*num_pos_b+1):2:(4*num_pos_b);
top_B_atoms=(2*num_pos_b+2):2:(4*num_pos_b);

% convert to basis number index
conv_A_bot = @(indd) bot_A_atoms(indd);
conv_B_bot = @(indd) bot_B_atoms(indd);
conv_A_top = @(indd) top_A_atoms(indd);
conv_B_top = @(indd) top_B_atoms(indd);

% nearest interlayer bond between Bravvis lattice points from bottom to
% top layer
inter_table=0;
inter_vec=0;
inter_angles=0;
inter_dist=0;

bottom_angle_a=angle(deltav_b(1)+i*deltav_b(2));
top_angle_a=angle(deltav_t(1)+i*deltav_t(2));

bottom_angle_b=bottom_angle_a+pi;
top_angle_b=top_angle_a+pi;

indnn=1;


% Now find interlayer coupling pairs
% THIS IS THE SLOW PART!!!
% Should use the inversion trick to find the 9 "search areas" for each b.
% Currently this scales N^2, but with the inversion trick it scales like N

% new_inter_search: type of interlayer search
% 0: Old version: exhaustive search of every lattice site
% 1: New version: checks only a local region for coupled lattice sites
new_inter_search = 1;

% New version
if (new_inter_search)

    sc_search_width = 1; % number of supercell-equivalent sites to check
    index_search_width = 6;
    for indb=1:num_pos_b

        % This computes the positon of lattice number indb on bottom layer
        Rvec_b_int=inverse_index_mat_b(indb,:);
        Rvec_b=Rvec_b_int(1)*a1+Rvec_b_int(2)*a2;

        % Check all supercell equivalent sites
        for sc_x = -sc_search_width:sc_search_width
            for sc_y = -sc_search_width:sc_search_width

                % compute the new position
                sc_shift = sc_t1*sc_x + sc_t2*sc_y;
                pos_here = Rvec_b + sc_shift;
                % Rotate this into the top-layer's basis (transpose of rot mat)
                rot_pos = rot_mat'*pos_here;
                % find the L2 lattice closest to this site!
                top_grid = floor(inv([a1 a2])*rot_pos)';

                % Now loop over all posible grid indexing of the lattice
                for d_i = -index_search_width:index_search_width
                    for d_j = -index_search_width:index_search_width

                        top_grid_here = top_grid + [d_i d_j];
                        % And correct for how index_mat_t is arranged
                        top_grid_here_mat = top_grid_here + [sc_n+1 2];

                        % check if grid is within size of index_mat_t
                        if (top_grid_here_mat(1) < 1 || top_grid_here_mat(1) > sc_n+sc_m+1)
                            continue
                        end
                        if (top_grid_here_mat(2) < 1 || top_grid_here_mat(2) > 2*sc_n+sc_m+2)
                            continue
                        end

                        index_here = index_mat_t(top_grid_here_mat(1),top_grid_here_mat(2));
                        indt = index_here;

                        % If this index exists, we add it!
                        if index_here ~= 0
                            

                            Rvec_t = top_grid_here(1)*ra1 + top_grid_here(2)*ra2;
                            
                            % debug stuff
                            %{
                            if (abs(pos_here(1) - 0) < .1)
                                fprintf("------------------------------------------\n");
                                fprintf("pos_bot = [%f, %f], pos_top = [%f, %f] \n",pos_here(1),pos_here(2),Rvec_t(1),Rvec_t(2));
                                fprintf("top_grid = [%d, %d] \n",top_grid(1)+d_i,top_grid(2)+d_j);
                            end
                            %}
                            
                            % 0 = A, 1 = B
                            for o1 = 0:1
                                for o2 = 0:1
                                    minvec = (Rvec_t + shift_top_vec + o2*deltav_t) - (pos_here + shift_bot_vec + o1*deltav_b);
                                    mindist = norm(minvec);

                                    if mindist<=inter_cutoff

                                        if     (o1 == 0 && o2 == 0) % A to A
                                            inter_table(indnn,1:2)=[conv_A_top(indt),conv_A_bot(indb)];
                                        elseif (o1 == 1 && o2 == 1) % B to B
                                            inter_table(indnn,1:2)=[conv_B_top(indt),conv_B_bot(indb)];                                            
                                        elseif (o1 == 0 && o2 == 1) % A to B
                                            inter_table(indnn,1:2)=[conv_B_top(indt),conv_A_bot(indb)]; 
                                        elseif (o1 == 1 && o2 == 0) % B to A
                                            inter_table(indnn,1:2)=[conv_A_top(indt),conv_B_bot(indb)];  
                                        end
                                        
                                        inter_vec(indnn,1:2)=minvec(:);
                                        angle3=angle(minvec(1)+1j*minvec(2));
                                        inter_dist(indnn)=mindist;
                                        inter_angles(indnn,1:3)=[top_angle_a,bottom_angle_a,angle3];
                                        indnn=indnn+1;
                                    end

                                end
                            end % orbit loop

                        end


                    end
                end % grid index loop

            end
        end % supercell loop


    end % done with interlayer couplings!!

else
% Old version (very slow!)

    num_list=9;
    for indb=1:num_pos_b

        Rvec_b_int=inverse_index_mat_b(indb,:);
        Rvec_b=Rvec_b_int(1)*a1+Rvec_b_int(2)*a2;

        for indt=1:num_pos_t
            Rvec_t_int=inverse_index_mat_t(indt,:);
            Rvec_t=Rvec_t_int(1)*ra1+Rvec_t_int(2)*ra2;

            % A to A

            [listdist,period_vecs] = list_vec( sc_t1,sc_t2,Rvec_b+shift_bot_vec,Rvec_t+shift_top_vec );

            for indl=1:num_list
                mindist=listdist(indl);
                minvec=period_vecs(indl,:);
                if mindist<=inter_cutoff
                    inter_table(indnn,1:2)=[conv_A_top(indt),conv_A_bot(indb)];
                    inter_vec(indnn,1:2)=minvec(:);
                    angle3=angle(minvec(1)+i*minvec(2));
                    inter_dist(indnn)=mindist;
                    inter_angles(indnn,1:3)=[top_angle_a,bottom_angle_a,angle3];
                    indnn=indnn+1;
                end
            end

            % B to B
            [listdist,period_vecs] = list_vec( sc_t1,sc_t2,Rvec_b+deltav_b+shift_bot_vec,Rvec_t+deltav_t+shift_top_vec );

            for indl=1:num_list
                mindist=listdist(indl);
                minvec=period_vecs(indl,:);
                if mindist<=inter_cutoff
                    inter_table(indnn,1:2)=[conv_B_top(indt),conv_B_bot(indb)];
                    inter_vec(indnn,1:2)=minvec(:);
                    angle3=angle(minvec(1)+i*minvec(2));
                    inter_dist(indnn)=mindist;
                    inter_angles(indnn,1:3)=[top_angle_b,bottom_angle_b,angle3];
                    indnn=indnn+1;
                end
            end

            % A to B
            [listdist,period_vecs] = list_vec( sc_t1,sc_t2,Rvec_b+shift_bot_vec,Rvec_t+deltav_t+shift_top_vec );

            for indl=1:num_list
                mindist=listdist(indl);
                minvec=period_vecs(indl,:);
                if mindist<=inter_cutoff
                    inter_table(indnn,1:2)=[conv_B_top(indt),conv_A_bot(indb)];
                    inter_vec(indnn,1:2)=minvec(:);
                    angle3=angle(minvec(1)+i*minvec(2));
                    inter_dist(indnn)=mindist;
                    inter_angles(indnn,1:3)=[top_angle_b,bottom_angle_a,angle3];
                    indnn=indnn+1;
                end
            end

            % B to A
            [listdist,period_vecs] = list_vec( sc_t1,sc_t2,Rvec_b+deltav_b+shift_bot_vec,Rvec_t+shift_top_vec );

            for indl=1:num_list
                mindist=listdist(indl);
                minvec=period_vecs(indl,:);
                if mindist<=inter_cutoff
                    inter_table(indnn,1:2)=[conv_A_top(indt),conv_B_bot(indb)];
                    inter_vec(indnn,1:2)=minvec(:);
                    angle3=angle(minvec(1)+i*minvec(2));
                    inter_dist(indnn)=mindist;
                    inter_angles(indnn,1:3)=[top_angle_a,bottom_angle_b,angle3];
                    indnn=indnn+1;
                end
            end


        end
    end

end

fprintf("Average number of inter connections per lattice site: %f \n",(indnn-1)/(num_pos_b));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% all lattice points

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
