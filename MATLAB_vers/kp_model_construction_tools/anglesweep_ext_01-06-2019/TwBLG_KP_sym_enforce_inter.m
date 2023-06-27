function [All_Eff_inter_sym,All_Eff_inter_kplus_sym,All_Eff_inter_kminus_sym] = TwBLG_KP_sym_enforce_inter(All_Eff_inter_init,All_Eff_inter_kplus_init,All_Eff_inter_kminus_init)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


% basic variables
sigx=[0,1;1,0];
sigy=[0,-i;i,0];
sigz=[1,0;0,-1];

% enforce the symmetry

All_Eff_inter_C2T=0*All_Eff_inter_init;

All_Eff_inter_C2TM=0*All_Eff_inter_init;

All_Eff_inter_C2TMR3=0*All_Eff_inter_init;

for indq=1:12
    tmpM=squeeze(All_Eff_inter_init(:,:,indq));
    All_Eff_inter_C2T(:,:,indq)=(tmpM+sigx*conj(tmpM)*(sigx))/2;
    
end

Mirror_ind=[1,3,2,5,4,6,8,7,11,12,9,10];

for indq=1:12
    tmpM=squeeze(All_Eff_inter_init(:,:,indq));
    tmpM2=squeeze(All_Eff_inter_init(:,:,Mirror_ind(indq)));
    
    All_Eff_inter_C2TM(:,:,indq)=(tmpM+sigx*(tmpM2')*sigx)/2;
    
end


rot_inter=[3,1,2,5,6,4,9,10,12,11,8,7];

rot_inter_inv=[2,3,1,6,4,5,12,11,7,8,10,9];

rot_mat=diag([exp(-i*2*pi/3),exp(i*2*pi/3)]);

for indq=1:12
    Mtmp=squeeze(All_Eff_inter_C2TM(:,:,indq));
    Mtmp1=squeeze(All_Eff_inter_C2TM(:,:,rot_inter(indq)));
    Mtmp1b=squeeze(All_Eff_inter_C2TM(:,:,rot_inter_inv(indq)));
    All_Eff_inter_C2TMR3(:,:,indq)=(Mtmp+(rot_mat')*Mtmp1*rot_mat+(rot_mat)*Mtmp1b*(rot_mat'))/3;
    
end

%All_Eff_inter_sym=All_Eff_inter_C2TMR3; 
% constant terms are done during projection now
All_Eff_inter_sym = All_Eff_inter_init;

% for 4:12 momentum scattering, they will be removed.
All_Eff_inter_kplus_sym=0*All_Eff_inter_kplus_init;
All_Eff_inter_kminus_sym=0*All_Eff_inter_kminus_init;

% miller relabeling
miller_idx = [2 3 1];

Mplus1=squeeze(All_Eff_inter_kplus_init(:,:,miller_idx(1)));
Mplus2=squeeze(All_Eff_inter_kplus_init(:,:,miller_idx(2)));
Mplus3=squeeze(All_Eff_inter_kplus_init(:,:,miller_idx(3)));

Mminus1=squeeze(All_Eff_inter_kminus_init(:,:,miller_idx(1)));
Mminus2=squeeze(All_Eff_inter_kminus_init(:,:,miller_idx(2)));
Mminus3=squeeze(All_Eff_inter_kminus_init(:,:,miller_idx(3)));

phi=exp(i*2*pi/3);

Mplus_phase1=[1,1;1,1];
Mplus_phase2=[phi,1;conj(phi),phi];
Mplus_phase3=[conj(phi),1;(phi),conj(phi)];

Mminus_phase1=[1,1;1,1];
Mminus_phase2=[conj(phi),phi;1,conj(phi)];
Mminus_phase3=[(phi),conj(phi);1,(phi)];

Mplus_avg=(Mplus1.*conj(Mplus_phase1)+Mplus2.*conj(Mplus_phase2)+Mplus3.*conj(Mplus_phase3))/3;
Mminus_avg=(Mminus1.*conj(Mminus_phase1)+Mminus2.*conj(Mminus_phase2)+Mminus3.*conj(Mminus_phase3))/3;

Kdep_parm1=(Mplus_avg(1,1)-conj(Mplus_avg(2,2)-Mminus_avg(1,1)+conj(Mminus_avg(2,2))))/4;
Kdep_parm2=(Mplus_avg(2,1)-Mminus_avg(1,2))/2;
Kdep_parm3=(Mplus_avg(1,2)-Mminus_avg(2,1))/2;

Kdep_parm1=i*imag(Kdep_parm1);
Kdep_parm2=i*imag(Kdep_parm2);
Kdep_parm3=i*imag(Kdep_parm3);

Msym_plus=[Kdep_parm1,Kdep_parm3;Kdep_parm2,Kdep_parm1];
Msym_minus=[-Kdep_parm1,-Kdep_parm2;-Kdep_parm3,-Kdep_parm1];

All_Eff_inter_kplus_sym(:,:,miller_idx(1))=Msym_plus.*Mplus_phase1;
All_Eff_inter_kplus_sym(:,:,miller_idx(2))=Msym_plus.*Mplus_phase2;
All_Eff_inter_kplus_sym(:,:,miller_idx(3))=Msym_plus.*Mplus_phase3;

All_Eff_inter_kminus_sym(:,:,miller_idx(1))=Msym_minus.*Mminus_phase1;
All_Eff_inter_kminus_sym(:,:,miller_idx(2))=Msym_minus.*Mminus_phase2;
All_Eff_inter_kminus_sym(:,:,miller_idx(3))=Msym_minus.*Mminus_phase3;


end

