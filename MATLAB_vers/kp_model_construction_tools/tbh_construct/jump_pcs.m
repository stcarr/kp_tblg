function [ new_ind ] = jump_pcs( start_ind,actual_index_mat,offset_inds,inverse_index_mat,jump_seq )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

new_ind=start_ind;
new_ij=inverse_index_mat(start_ind,:);

tmp=length(jump_seq);


for inds=1:tmp
    now_ind=new_ind;
    now_ij=new_ij;
    
    new_ind=actual_index_mat(now_ij(1)-offset_inds(1)+1,now_ij(2)-offset_inds(2)+1,jump_seq(inds));
    new_ij=inverse_index_mat(new_ind,:);
    
    
end



end

