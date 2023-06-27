function [ newvec ] = moveto_primitive_method2( vec1,vec2,targetvec )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% vec1,vec2 components are "integers"

n1vec=[-vec1(2),vec1(1)];

n2vec=[vec2(2),-vec2(1)];

bvd1=dot(n1vec,vec2);
bvd2=dot(n2vec,vec1);

proj1=dot(n1vec,targetvec);
proj2=dot(n2vec,targetvec);

mod1part=mod(proj1,bvd1);
mod2part=mod(proj2,bvd2);

int1part=(proj1-mod1part)/bvd1;
int2part=(proj2-mod2part)/bvd2;

newvec=targetvec-int1part*vec2'-int2part*vec1';


end

