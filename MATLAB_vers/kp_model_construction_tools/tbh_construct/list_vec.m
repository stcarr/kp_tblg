function [listdist,period_vecs] = list_vec( sc1,sc2,vecb,vect )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

delvec=vect-vecb;

period_vecs=zeros(9,2);
period_vecs(1,:)=delvec;
period_vecs(2,:)=sc1+delvec;
period_vecs(3,:)=-sc1+delvec;
period_vecs(4,:)=sc2+delvec;
period_vecs(5,:)=-sc2+delvec;
period_vecs(6,:)=sc1+sc2+delvec;
period_vecs(7,:)=-sc1-sc2+delvec;
period_vecs(8,:)=sc1-sc2+delvec;
period_vecs(9,:)=-sc1+sc2+delvec;

% !
listdist=sqrt(sum(period_vecs.^2,2));

end

