function [ Hmat] = ham_supercell_N8_Koshino_relax_height( know, total_num, inplaneGAB_bottom_relax,inplaneGAB_top_relax, inplaneGAABB_bottom_relax,inplaneGAABB_top_relax,inter_relax ,graphene_onsite_e)

% index_to, index_from, hop_t, vec_x, vec_y


% inplane A to B only!
% interaction, bottom to top layer only



Hmat=zeros(total_num);
Hmat_inplane=zeros(total_num);
Hmat_interaction=zeros(total_num);

% index_bottom=1:(total_num/2);
% index_top=(total_num/2+1):total_num;

% inplane part

% 3 nearest bonds and for 2 layers and 3 3rd nearest bonds

sizetmp=size(inplaneGAB_bottom_relax);

for ind=1:sizetmp(3)
    delvec_b=squeeze(inplaneGAB_bottom_relax(:,4:5,ind));
    delvec_t=squeeze(inplaneGAB_top_relax(:,4:5,ind));
    
    tmp_table_b=squeeze(inplaneGAB_bottom_relax(:,1:2,ind));
    tmp_table_t=squeeze(inplaneGAB_top_relax(:,1:2,ind));
    
    intra_coupling_b=squeeze(inplaneGAB_bottom_relax(:,3,ind));
    intra_coupling_t=squeeze(inplaneGAB_top_relax(:,3,ind));
    
    
    % convert index
    tmp_index_b=total_num*(tmp_table_b(:,2)-1)+tmp_table_b(:,1);
    tmp_index_t=total_num*(tmp_table_t(:,2)-1)+tmp_table_t(:,1);
    
    %size(intra_coupling_b)
    %size(delvec_b)
    %size(know)
    
    Hmat_inplane(tmp_index_b)=Hmat_inplane(tmp_index_b)+intra_coupling_b.*exp(-i*delvec_b*(know'));
    Hmat_inplane(tmp_index_t)=Hmat_inplane(tmp_index_t)+intra_coupling_t.*exp(-i*delvec_t*(know'));
    
end
Hmat_inplane=(Hmat_inplane+Hmat_inplane');



% inplane 2nd nearest bonds to be added!

sizetmp=size(inplaneGAABB_bottom_relax);

for ind=1:sizetmp(3)
    delvec_b=squeeze(inplaneGAABB_bottom_relax(:,4:5,ind));
    delvec_t=squeeze(inplaneGAABB_top_relax(:,4:5,ind));
    
    tmp_table_b=squeeze(inplaneGAABB_bottom_relax(:,1:2,ind));
    tmp_table_t=squeeze(inplaneGAABB_top_relax(:,1:2,ind));
    
    intra_coupling_b=squeeze(inplaneGAABB_bottom_relax(:,3,ind));
    intra_coupling_t=squeeze(inplaneGAABB_top_relax(:,3,ind));

    % convert index
    tmp_index_b=total_num*(tmp_table_b(:,2)-1)+tmp_table_b(:,1);
    tmp_index_t=total_num*(tmp_table_t(:,2)-1)+tmp_table_t(:,1);
    
    Hmat_inplane(tmp_index_b)=Hmat_inplane(tmp_index_b)+intra_coupling_b.*exp(-i*delvec_b*(know'));
    Hmat_inplane(tmp_index_t)=Hmat_inplane(tmp_index_t)+intra_coupling_t.*exp(-i*delvec_t*(know'));
    
end

% diagonal part? to be added later!
Hmat_inplane=Hmat_inplane+eye(total_num)*graphene_onsite_e;


% interaction part

tmp_index=total_num*(inter_relax(:,2)-1)+inter_relax(:,1);

%Hmat_interaction(tmp_index)=Hmat_interaction(tmp_index)+inter_hop.*exp(-i*inter_vec*know');
intsize=size(inter_relax);
for inds=1:intsize(1)
	Hmat_interaction(tmp_index(inds))=Hmat_interaction(tmp_index(inds))+inter_relax(inds,3)*exp(-i*inter_relax(inds,4:5)*know');
end


Hmat_interaction=Hmat_interaction+Hmat_interaction';

 Hmat=Hmat_inplane+Hmat_interaction;
 
 

end

