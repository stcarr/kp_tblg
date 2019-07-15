clear all;

side = 0; % 0: top, 1: bot

% have to load the correct WF file here (should be REWAN output)
if (side == 0)
    load('run_folder_5band/fiveWANtopnew_wf_1p10_16k.mat')
elseif (side == 1)
    load('run_folder_5band/fiveWANbotnew_wf_1p10_16k.mat')
end
WF_phase = [1 1 1 1 1];

% use to see if pre-REWAN WFs need phase fixing to get proper symmetry
% this is fixed during the REWAN step, so not needed if using REWAN WFs
%{
    load('run_folder_5band/fiveWANtopnew_wf_1p10_16k.mat')
    if (side == 0)
        WF_phase = [exp(1j*pi/4),exp(1j*3*pi/4),1,1j,1j];
    elseif (side == 1)
        WF_phase = [exp(1j*pi/4),exp(1j*3*pi/4),1,-1j,1j];        
    end
%}


wf1 = 1;
wf2 = 2;

ax_m = 100;
c_max = 0.05;

nx = size(all_wfAL1,1);
all_WF = zeros(nx,nx,5,4);
all_WF(:,:,:,1) = all_wfAL1;
all_WF(:,:,:,2) = all_wfBL1;
all_WF(:,:,:,3) = all_wfAL2;
all_WF(:,:,:,4) = all_wfBL2;

%
% Rules for how the sublattice/layer DoFs change under symm operation
%

% in-plane 3-fold rotation
% C3: All orbitals to themselves...

% x to -x mirrory sym
% Mx: AL1 <-> BL2
%     BL1 <-> AL2

% in-plane 2-fold rotation w/ conjugate
% C2T: AL1 <-> BL1
%      AL2 <-> BL2

tar_wf1 = squeeze(WF_phase(wf1)*all_WF(:,:,wf1,:));
tar_wf2 = squeeze(WF_phase(wf2)*all_WF(:,:,wf2,:));

titles = {'$A L_1$','$B L_1$', '$A L_2$', '$B L_2$'};

clf
for orb = 1:4
    subplot(4,2,orb)
    surf(mmX,mmY,real(tar_wf1(:,:,orb)))
    view(2)
    shading interp
    axis equal
    axis([-ax_m ax_m -ax_m ax_m])
    colorbar
    caxis([-c_max c_max])
    title(['WF1 ' titles{orb}])
    set(gca,'XTick',[])

end

for orb = 1:4
    subplot(4,2,4+orb)
    surf(mmX,mmY,real(tar_wf2(:,:,orb)))
    view(2)
    shading interp
    axis equal
    axis([-ax_m ax_m -ax_m ax_m])
    colorbar
    caxis([-c_max c_max])
    title(['WF2 ' titles{orb}])
    set(gca,'XTick',[])

end

colormap default
