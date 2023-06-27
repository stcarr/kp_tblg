% build supercell hamiltonian

% variables to be used:
% sc_m, sc_n,



% goal is to set up the parameters to run supercell_hamiltonian

% total_num,
% inplane_table_bottom,
% inplane_vecs_bottom,
% inplane_table_top,
% inplane_vecs_top,
% inter_table,
% inter_vec,
% inter_hop


% setup variable
if num_pos_t~=num_pos_b
    'warning! something is wrong...'
end

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

% inplane nearest neighbor interaction
inplane_table_bottom=zeros(num_pos_b,2,3);
inplane_table_top=zeros(num_pos_t,2,3);

% bottom layer
inplane_vecs_bottom13=zeros(6,2);
inplane_vecs_top13=zeros(6,2);
inplane_vecs_bottom2=zeros(6,2);
inplane_vecs_top2=zeros(6,2);

deltav=(a1+a2)/3;
% N1
inplaneGAB_vecs_bottom(1,:)=deltav;
inplaneGAB_vecs_bottom(2,:)=deltav-a1;
inplaneGAB_vecs_bottom(3,:)=deltav-a2;

% N3
inplaneGAB_vecs_bottom(4,:)=deltav+a2-a1;
inplaneGAB_vecs_bottom(5,:)=deltav+a1-a2;
inplaneGAB_vecs_bottom(6,:)=deltav-a1-a2;

% N4
inplaneGAB_vecs_bottom(7,:)=deltav+a1;
inplaneGAB_vecs_bottom(8,:)=deltav+a2;
inplaneGAB_vecs_bottom(9,:)=deltav+a2-2*a1;
inplaneGAB_vecs_bottom(10,:)=deltav-2*a1;
inplaneGAB_vecs_bottom(11,:)=deltav-2*a2;
inplaneGAB_vecs_bottom(12,:)=deltav+a1-2*a2;

% N7
inplaneGAB_vecs_bottom(13,:)=deltav+2*a1-2*a2;
inplaneGAB_vecs_bottom(14,:)=deltav+2*a1-a2;
inplaneGAB_vecs_bottom(15,:)=deltav-a1+2*a2;
inplaneGAB_vecs_bottom(16,:)=deltav-2*a1+2*a2;
inplaneGAB_vecs_bottom(17,:)=deltav-2*a1-a2;
inplaneGAB_vecs_bottom(18,:)=deltav-a1-2*a2;

% N8
inplaneGAB_vecs_bottom(19,:)=deltav+a1+a2;
inplaneGAB_vecs_bottom(20,:)=deltav-3*a1+a2;
inplaneGAB_vecs_bottom(21,:)=deltav+a1-3*a2;


% N2
inplaneGAABB_vecs_bottom(1,:)=a1;
inplaneGAABB_vecs_bottom(2,:)=a2;
inplaneGAABB_vecs_bottom(3,:)=a2-a1;
inplaneGAABB_vecs_bottom(4,:)=-a1;
inplaneGAABB_vecs_bottom(5,:)=-a2;
inplaneGAABB_vecs_bottom(6,:)=a1-a2;


% N5
inplaneGAABB_vecs_bottom(7,:)=a1+a2;
inplaneGAABB_vecs_bottom(8,:)=2*a2-a1;
inplaneGAABB_vecs_bottom(9,:)=a2-2*a1;
inplaneGAABB_vecs_bottom(10,:)=-a2-a1;
inplaneGAABB_vecs_bottom(11,:)=-2*a2+a1;
inplaneGAABB_vecs_bottom(12,:)=2*a1-a2;

% N6
inplaneGAABB_vecs_bottom(13,:)=2*a1;
inplaneGAABB_vecs_bottom(14,:)=2*a2;
inplaneGAABB_vecs_bottom(15,:)=2*a2-2*a1;
inplaneGAABB_vecs_bottom(16,:)=-2*a1;
inplaneGAABB_vecs_bottom(17,:)=-2*a2;
inplaneGAABB_vecs_bottom(18,:)=-2*a2+2*a1;


% top layer
deltav=(ra1+ra2)/3;
inplaneGAB_vecs_top(1,:)=deltav;
inplaneGAB_vecs_top(2,:)=deltav-ra1;
inplaneGAB_vecs_top(3,:)=deltav-ra2;
inplaneGAB_vecs_top(4,:)=deltav+ra2-ra1;
inplaneGAB_vecs_top(5,:)=deltav+ra1-ra2;
inplaneGAB_vecs_top(6,:)=deltav-ra1-ra2;
inplaneGAB_vecs_top(7,:)=deltav+ra1;
inplaneGAB_vecs_top(8,:)=deltav+ra2;
inplaneGAB_vecs_top(9,:)=deltav+ra2-2*ra1;
inplaneGAB_vecs_top(10,:)=deltav-2*ra1;
inplaneGAB_vecs_top(11,:)=deltav-2*ra2;
inplaneGAB_vecs_top(12,:)=deltav+ra1-2*ra2;

% N7
inplaneGAB_vecs_top(13,:)=deltav+2*ra1-2*ra2;
inplaneGAB_vecs_top(14,:)=deltav+2*ra1-ra2;
inplaneGAB_vecs_top(15,:)=deltav-ra1+2*ra2;
inplaneGAB_vecs_top(16,:)=deltav-2*ra1+2*ra2;
inplaneGAB_vecs_top(17,:)=deltav-2*ra1-ra2;
inplaneGAB_vecs_top(18,:)=deltav-ra1-2*ra2;

% N8
inplaneGAB_vecs_top(19,:)=deltav+ra1+ra2;
inplaneGAB_vecs_top(20,:)=deltav-3*ra1+ra2;
inplaneGAB_vecs_top(21,:)=deltav+ra1-3*ra2;

inplaneGAABB_vecs_top(1,:)=ra1;
inplaneGAABB_vecs_top(2,:)=ra2;
inplaneGAABB_vecs_top(3,:)=ra2-ra1;
inplaneGAABB_vecs_top(4,:)=-ra1;
inplaneGAABB_vecs_top(5,:)=-ra2;
inplaneGAABB_vecs_top(6,:)=ra1-ra2;
inplaneGAABB_vecs_top(7,:)=ra1+ra2;
inplaneGAABB_vecs_top(8,:)=2*ra2-ra1;
inplaneGAABB_vecs_top(9,:)=ra2-2*ra1;
inplaneGAABB_vecs_top(10,:)=-ra2-ra1;
inplaneGAABB_vecs_top(11,:)=-2*ra2+ra1;
inplaneGAABB_vecs_top(12,:)=2*ra1-ra2;
inplaneGAABB_vecs_top(13,:)=2*ra1;
inplaneGAABB_vecs_top(14,:)=2*ra2;
inplaneGAABB_vecs_top(15,:)=2*ra2-2*ra1;
inplaneGAABB_vecs_top(16,:)=-2*ra1;
inplaneGAABB_vecs_top(17,:)=-2*ra2;
inplaneGAABB_vecs_top(18,:)=-2*ra2+2*ra1;


% iterate over number of primitive cells

% Nearest neighbor coupling: inplane_table_bottom,inplane_table_top

for indp=1:num_pos_b
    
    bottom_i=inverse_index_mat_b(indp,1);
    bottom_j=inverse_index_mat_b(indp,2);
    
    top_i=inverse_index_mat_t(indp,1);
    top_j=inverse_index_mat_t(indp,2);
    
    % type 1
    % from
    inplane_table_bottom(indp,2,1)=conv_A_bot(indp);
    % to
    inplane_table_bottom(indp,1,1)=conv_B_bot(indp);
    % from
    inplane_table_top(indp,2,1)=conv_A_top(indp);
    % to
    inplane_table_top(indp,1,1)=conv_B_top(indp);
    
    
    % type 2 (to unit [-1,0])
    % from
    inplane_table_bottom(indp,2,2)=conv_A_bot(indp);
    tt=access_neighbor_bottom(bottom_i,bottom_j,5);
    % to
    inplane_table_bottom(indp,1,2)=conv_B_bot(tt);
    
    % from
    inplane_table_top(indp,2,2)=conv_A_top(indp);
    tt=access_neighbor_top(top_i,top_j,5);
    % to
    inplane_table_top(indp,1,2)=conv_B_top(tt);
    
    % type 3 (to unit [0,-1])
    % from
    inplane_table_bottom(indp,2,3)=conv_A_bot(indp);
    tt=access_neighbor_bottom(bottom_i,bottom_j,6);
    % to
    inplane_table_bottom(indp,1,3)=conv_B_bot(tt);
    
    % from
    inplane_table_top(indp,2,3)=conv_A_top(indp);
    tt=access_neighbor_top(top_i,top_j,6);
    % to
    inplane_table_top(indp,1,3)=conv_B_top(tt);
    
    
end

% in-plane 2nd nearest neighbors
inplane_NN_table_bottom=zeros(2*num_pos_b,2,6);
inplane_NN_table_top=zeros(2*num_pos_t,2,6);

% 4th neighbor
inplane_N4_table_bottom=zeros(num_pos_b,2,6);
inplane_N4_table_top=zeros(num_pos_t,2,6);

% 5th neighbor
inplane_N5_table_bottom=zeros(2*num_pos_b,2,6);
inplane_N5_table_top=zeros(2*num_pos_t,2,6);

% 6th neighbor
inplane_N6_table_bottom=zeros(2*num_pos_b,2,6);
inplane_N6_table_top=zeros(2*num_pos_t,2,6);



N4_table=[0,0;
    3,7;
    2,4;
    4,5;
    5,5;
    6,6;
    6,7];


for indp=1:num_pos_b
    %     bottom_i=inverse_index_mat_b(indp,1);
    %     bottom_j=inverse_index_mat_b(indp,2);
    %
    %     top_i=inverse_index_mat_t(indp,1);
    %     top_j=inverse_index_mat_t(indp,2);
    
    indd1=2*(indp-1)+1;
    indd2=2*(indp);
    
    for indu=2:7
        
        % bottom/top
        [newbotu]=jump_pcs( indp,actual_index_mat_b,offset_inds_bottom,inverse_index_mat_b,[indu] );
        [newtopu]=jump_pcs( indp,actual_index_mat_t,offset_inds_top,inverse_index_mat_t,[indu] );
        
        % from
        inplane_NN_table_bottom(indd1,2,indu-1)=conv_A_bot(indp);
        % to
        inplane_NN_table_bottom(indd1,1,indu-1)=conv_A_bot(newbotu);
        
        % from
        inplane_NN_table_top(indd1,2,indu-1)=conv_A_top(indp);
        % to
        inplane_NN_table_top(indd1,1,indu-1)=conv_A_top(newtopu);
        
        
        % from
        inplane_NN_table_bottom(indd2,2,indu-1)=conv_B_bot(indp);
        % to
        inplane_NN_table_bottom(indd2,1,indu-1)=conv_B_bot(newbotu);
        
        % from
        inplane_NN_table_top(indd2,2,indu-1)=conv_B_top(indp);
        % to
        inplane_NN_table_top(indd2,1,indu-1)=conv_B_top(newtopu);
        
        % 6th NN
        % bottom/top
        [newbotu]=jump_pcs( indp,actual_index_mat_b,offset_inds_bottom,inverse_index_mat_b,[indu,indu] );
        [newtopu]=jump_pcs( indp,actual_index_mat_t,offset_inds_top,inverse_index_mat_t,[indu,indu] );
        
        % from
        inplane_N6_table_bottom(indd1,2,indu-1)=conv_A_bot(indp);
        % to
        inplane_N6_table_bottom(indd1,1,indu-1)=conv_A_bot(newbotu);
        
        % from
        inplane_N6_table_top(indd1,2,indu-1)=conv_A_top(indp);
        % to
        inplane_N6_table_top(indd1,1,indu-1)=conv_A_top(newtopu);
        
        
        % from
        inplane_N6_table_bottom(indd2,2,indu-1)=conv_B_bot(indp);
        % to
        inplane_N6_table_bottom(indd2,1,indu-1)=conv_B_bot(newbotu);
        
        % from
        inplane_N6_table_top(indd2,2,indu-1)=conv_B_top(indp);
        % to
        inplane_N6_table_top(indd2,1,indu-1)=conv_B_top(newtopu);
        
        
        % 5th NN
        % bottom/top
        if indu<7
            [newbotu]=jump_pcs( indp,actual_index_mat_b,offset_inds_bottom,inverse_index_mat_b,[indu,indu+1] );
            [newtopu]=jump_pcs( indp,actual_index_mat_t,offset_inds_top,inverse_index_mat_t,[indu,indu+1] );
        else
            [newbotu]=jump_pcs( indp,actual_index_mat_b,offset_inds_bottom,inverse_index_mat_b,[7,2] );
            [newtopu]=jump_pcs( indp,actual_index_mat_t,offset_inds_top,inverse_index_mat_t,[7,2] );
        end
        % from
        inplane_N5_table_bottom(indd1,2,indu-1)=conv_A_bot(indp);
        % to
        inplane_N5_table_bottom(indd1,1,indu-1)=conv_A_bot(newbotu);
        
        % from
        inplane_N5_table_top(indd1,2,indu-1)=conv_A_top(indp);
        % to
        inplane_N5_table_top(indd1,1,indu-1)=conv_A_top(newtopu);
        
        
        % from
        inplane_N5_table_bottom(indd2,2,indu-1)=conv_B_bot(indp);
        % to
        inplane_N5_table_bottom(indd2,1,indu-1)=conv_B_bot(newbotu);
        
        % from
        inplane_N5_table_top(indd2,2,indu-1)=conv_B_top(indp);
        % to
        inplane_N5_table_top(indd2,1,indu-1)=conv_B_top(newtopu);
        
        
        
        % 4th NN   NOTE: this is an AB type!
        % bottom/top
        
        [newbotu]=jump_pcs( indp,actual_index_mat_b,offset_inds_bottom,inverse_index_mat_b,[N4_table(indu,1),N4_table(indu,2)] );
        [newtopu]=jump_pcs( indp,actual_index_mat_t,offset_inds_top,inverse_index_mat_t,[N4_table(indu,1),N4_table(indu,2)] );
        
        % from
        inplane_N4_table_bottom(indp,2,indu-1)=conv_A_bot(indp);
        % to
        inplane_N4_table_bottom(indp,1,indu-1)=conv_B_bot(newbotu);
        
        % from
        inplane_N4_table_top(indp,2,indu-1)=conv_A_top(indp);
        % to
        inplane_N4_table_top(indp,1,indu-1)=conv_B_top(newtopu);
        
        
        
        
        
    end
end


% in-plane 3rd nearest neighbors % can use hermitian conjugate trick! take
% care only A to B as the nearest neighbor terms inplane
inplane_NNN_table_bottom=zeros(num_pos_b,2,3);
inplane_NNN_table_top=zeros(num_pos_t,2,3);


NNN_table=[3,5;
    2,6;
    5,6];


% 7th neighbor
inplane_N7_table_bottom=zeros(num_pos_b,2,6);
inplane_N7_table_top=zeros(num_pos_t,2,6);

% 8th neighbor
inplane_N8_table_bottom=zeros(num_pos_b,2,3);
inplane_N8_table_top=zeros(num_pos_t,2,3);


% 7/8 AB type
N7_table=[7,2,6;
    2,2,6;
    5,3,3;
    4,5,3;
    5,5,6;
    5,6,6];

N8_table=[2,2,4;
    5,5,4;
    6,6,7];


for indp=1:num_pos_b
    
    % N3 type
    for indu=1:3
        % bottom/top
        [newbotu]=jump_pcs( indp,actual_index_mat_b,offset_inds_bottom,inverse_index_mat_b,[NNN_table(indu,1),NNN_table(indu,2)] );
        [newtopu]=jump_pcs( indp,actual_index_mat_t,offset_inds_top,inverse_index_mat_t,[NNN_table(indu,1),NNN_table(indu,2)] );
        
        % from
        inplane_NNN_table_bottom(indp,2,indu)=conv_A_bot(indp);
        % to
        inplane_NNN_table_bottom(indp,1,indu)=conv_B_bot(newbotu);
        
        % from
        inplane_NNN_table_top(indp,2,indu)=conv_A_top(indp);
        % to
        inplane_NNN_table_top(indp,1,indu)=conv_B_top(newtopu);
        
    end
    
    % N7 type
    for indu=1:6
        % bottom/top
        [newbotu]=jump_pcs( indp,actual_index_mat_b,offset_inds_bottom,inverse_index_mat_b,[N7_table(indu,1),N7_table(indu,2),N7_table(indu,3)] );
        [newtopu]=jump_pcs( indp,actual_index_mat_t,offset_inds_top,inverse_index_mat_t,[N7_table(indu,1),N7_table(indu,2),N7_table(indu,3)] );
        
        % from
        inplane_N7_table_bottom(indp,2,indu)=conv_A_bot(indp);
        % to
        inplane_N7_table_bottom(indp,1,indu)=conv_B_bot(newbotu);
        
        % from
        inplane_N7_table_top(indp,2,indu)=conv_A_top(indp);
        % to
        inplane_N7_table_top(indp,1,indu)=conv_B_top(newtopu);
        
    end
    
    % N8 type
    for indu=1:3
        % bottom/top
        [newbotu]=jump_pcs( indp,actual_index_mat_b,offset_inds_bottom,inverse_index_mat_b,[N8_table(indu,1),N8_table(indu,2),N8_table(indu,3)] );
        [newtopu]=jump_pcs( indp,actual_index_mat_t,offset_inds_top,inverse_index_mat_t,[N8_table(indu,1),N8_table(indu,2),N8_table(indu,3)] );
        
        % from
        inplane_N8_table_bottom(indp,2,indu)=conv_A_bot(indp);
        % to
        inplane_N8_table_bottom(indp,1,indu)=conv_B_bot(newbotu);
        
        % from
        inplane_N8_table_top(indp,2,indu)=conv_A_top(indp);
        % to
        inplane_N8_table_top(indp,1,indu)=conv_B_top(newtopu);
        
    end
    
    
    
    
    
end


% combine!
inplaneGAB_table_bottom=zeros(num_pos_b,2,21);
inplaneGAB_table_top=zeros(num_pos_t,2,21);

inplaneGAB_table_bottom(:,:,1:3)=inplane_table_bottom(:,:,:);
inplaneGAB_table_top(:,:,1:3)=inplane_table_top(:,:,:);
inplaneGAB_table_bottom(:,:,4:6)=inplane_NNN_table_bottom(:,:,:);
inplaneGAB_table_top(:,:,4:6)=inplane_NNN_table_top(:,:,:);
inplaneGAB_table_bottom(:,:,7:12)=inplane_N4_table_bottom(:,:,:);
inplaneGAB_table_top(:,:,7:12)=inplane_N4_table_top(:,:,:);
inplaneGAB_table_bottom(:,:,13:18)=inplane_N7_table_bottom(:,:,:);
inplaneGAB_table_top(:,:,13:18)=inplane_N7_table_top(:,:,:);
inplaneGAB_table_bottom(:,:,19:21)=inplane_N8_table_bottom(:,:,:);
inplaneGAB_table_top(:,:,19:21)=inplane_N8_table_top(:,:,:);


inplaneGAABB_table_bottom=zeros(2*num_pos_b,2,6*3);
inplaneGAABB_table_top=zeros(2*num_pos_t,2,6*3);

inplaneGAABB_table_bottom(:,:,1:6)=inplane_NN_table_bottom(:,:,:);
inplaneGAABB_table_top(:,:,1:6)=inplane_NN_table_top(:,:,:);
inplaneGAABB_table_bottom(:,:,7:12)=inplane_N5_table_bottom(:,:,:);
inplaneGAABB_table_top(:,:,7:12)=inplane_N5_table_top(:,:,:);
inplaneGAABB_table_bottom(:,:,13:18)=inplane_N6_table_bottom(:,:,:);
inplaneGAABB_table_top(:,:,13:18)=inplane_N6_table_top(:,:,:);


%
% inter_table=0;
% inter_vec=0;
% inter_angles=0;
% inter_dist=0; 


% goal: inter_hop
tmpsize=size(inter_table);
num_nn_int=tmpsize(1);


tmp_nn_vecs=zeros(num_nn_int,2);
tmp_nn_vecs(:,1:2)=inter_vec(:,1:2)*lattice_a;
% tmp_nn_vecs(:,3)=3.35;


%inter_hop=graphene_interlayer_hopping(tmp_nn_vecs ,inter_angles(:,2),inter_angles(:,1));



