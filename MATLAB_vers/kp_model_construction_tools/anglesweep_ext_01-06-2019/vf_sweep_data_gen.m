theta_range = [1.0 0.19];
theta_0 = 1.1;
alpha_0 = sqrt(1/3);
alpha = (theta_0./theta_range)*alpha_0;
alpha2 = alpha.^2;


alpha2_list = linspace(alpha2(1),alpha2(2),150);
tar_theta_list = [5:-.1:1.6, 1.5:-0.025:1.025, theta_0*(alpha_0./sqrt(alpha2_list))];

[vals,scaleaxis,kpts] = tblg_kp_calc_ext('vf_only',1','theta_list',tar_theta_list,'knum',61);

save('dft_full_inter_3_intra_relax_vf_02-09-2019.mat','tar_theta_list','vals','kpts')
%%
[bands,scaleaxis,allkpts] = tblg_kp_calc_ext('theta_list',[0.44]);
%%
clf
hold on

%b_c = bands_cutoff{1};
%plot(scaleaxis{1},b_c(:,1:123),'r')

b = bands{1};
plot(scaleaxis{1},b(:,1:123),'-k')
plot(scaleaxis{1},b(:,124:end),'-k')

ax_m = 0.03;
axis([0 1 -ax_m ax_m])
%%

%close all
clear all

f_size = 16;

set(groot, 'DefaultTextInterpreter', 'latex')
set(groot, 'DefaultLegendInterpreter', 'latex')
set(groot, 'DefaultAxesTickLabelInterpreter', 'latex')
set(0,'DefaultAxesFontSize',f_size)

fig = figure('pos',[50 300 600 400]);

f_size = 16;

set(groot, 'DefaultTextInterpreter', 'Latex')
set(groot, 'DefaultLegendInterpreter', 'Latex')
set(groot, 'DefaultAxesTickLabelInterpreter', 'Latex')
set(0,'DefaultAxesFontSize',f_size)

relax_color = [0 0 0];
norelax_color = [1 0 0];
bmd_color = [0 150 0]./255;
elec_type = '-';
hole_type = '--';

dk_idx = 2;
theta_axis = [0.245 1.5];


% hbar vf del_k = del_E
% vf = del_E / (hbar del_k)
hbar = 6.58211951*(10^-16); % hbar in eV seconds
Dirac_v= 5.2275;
v0 = Dirac_v / (hbar*10^10); % fermi velocity for graphene
% if v_f = 1e6 m/s, then we would have Dirac_v = 6.582

%load dft_full_relax_fullsweep_vf
%load vf_data_larger_cutoff_02-06-2019.mat
%load dft_full_relax_vf_02-08-2019.mat
load dft_full_relax_vf_02-09-2019.mat
%load dft_3_inter_3_intra_relax_vf_02-09-2019.mat
%load dft_3_inter_full_intra_relax_vf_02-09-2019.mat
%load dft_full_inter_3_intra_relax_vf_02-09-2019.mat

%load dft_full_relax_fullsweep_vf_10-19-2018_6inter
%load dft_full_relax_fullsweep_vf_10-19-2018_3inter
sweep_kpts1 = kpts;
sweep_bands1 = vals;

offsets = [0.5 0.25 0];
vf_max = 0.9;
%alpha_max = 15;

for theta_idx = 1:length(tar_theta_list)
    kpts = sweep_kpts1{theta_idx};
    bands = sweep_bands1{theta_idx};
    nb = size(bands,1);
    
    start_idx = nb/2;
    for trial_idx = start_idx-1:start_idx+1
        if ( (bands(trial_idx+1,1) - bands(trial_idx,1)) < 1e-10 )
           start_idx = trial_idx;
           break
        end
    end
    
    d_k(theta_idx) = norm(kpts(1,:) - kpts(2,:));
    d_E_hole(theta_idx) = bands(start_idx,1) - bands(start_idx,2);
    d_E_elec(theta_idx) = bands(start_idx+1,1) - bands(start_idx+1,2);
    
end

vf_elec_relax = d_E_elec./(hbar*d_k)*10^-10; % angstroms to meters
vf_hole_relax = d_E_hole./(hbar*d_k)*10^-10;


theta_0 = 1.1;
alpha_0 = sqrt(1/3);
alpha = (theta_0./tar_theta_list)*alpha_0;
alpha2 = alpha.^2;
alpha2_relax = alpha2;

l_theta = 1.43*sqrt(3)./(2*sind(tar_theta_list/2));
k_theta = (2/3)*2*pi./l_theta;
alpha_theta = (0.11/Dirac_v)./k_theta;
alpha2_relax = alpha_theta.^2;

l_axis = 1.43*sqrt(3)./(2*sind(theta_axis/2));
k_axis = (2/3)*2*pi./l_axis;
alpha_axis = (0.11/Dirac_v)./k_axis;
alpha2_axis = alpha_axis.^2;
%alpha_max = alpha2_axis(1);
alpha_max = 8.5;

plot(alpha2_relax,vf_elec_relax,'Color',relax_color,'LineWidth',2)

%{
alpha2_tick_labels = [0:2:12];
alpha_tick_labels = sqrt(alpha2_tick_labels);
k_tick_labels = (0.11/Dirac_v)./alpha_tick_labels;
l_tick_labels = (2/3)*2*pi./k_tick_labels;
theta_tick_labels = 2*asind((1.43*sqrt(3)./l_tick_labels)/2);
%}

load vf_data/dft_no_relax_fullsweep_vf_11-02-2018
%load dft_no_relax_fullsweep
clear d_k
clear d_E_hole
clear d_E_elec

for theta_idx = 1:length(tar_theta_list)
    kpts = sweep_kpts1{theta_idx};
    bands = sweep_bands1{theta_idx};
    nb = size(bands,1);
    
    start_idx = nb/2;
    for trial_idx = start_idx-1:start_idx+1
        if ( (bands(trial_idx+1,1) - bands(trial_idx)) < 1e-7 )
           start_idx = trial_idx;
           break
        end
    end
    
    d_k(theta_idx) = norm(kpts(end,:) - kpts(end-dk_idx,:));
    d_E_hole(theta_idx) = bands(start_idx,end) - bands(start_idx,end-dk_idx);
    d_E_elec(theta_idx) = bands(start_idx+1,end) - bands(start_idx+1,end-dk_idx);
    
    bands2 = sweep_bands2{theta_idx};
    
    n4_elec = max([bands(start_idx+1,:) bands2(start_idx+1,:)]);
    min_elec = min([bands(start_idx+2,:) bands2(start_idx+2,:)]);
    n4_hole = min([bands(start_idx,:) bands2(start_idx,:)]);
    max_hole = max([bands(start_idx-1,:) bands2(start_idx-1,:)]);
    
    gap_elec(theta_idx) = min_elec - n4_elec;
    gap_hole(theta_idx) = n4_hole - max_hole;
    
    disp_elec(theta_idx) = range(bands(start_idx+1,:));
    disp_hole(theta_idx) = range(bands(start_idx,:));
    
end

vf_elec_norelax = d_E_elec./(hbar*d_k)*10^-10; % angstroms to meters
vf_hole_norelax = d_E_hole./(hbar*d_k)*10^-10;

l_theta = 1.43*sqrt(3)./(2*sind(tar_theta_list/2));
k_theta = (2/3)*2*pi./l_theta;
alpha_theta = (0.11/Dirac_v)./k_theta;
alpha2_norelax = alpha_theta.^2;


load vf_data/bmd_fullsweep_vf_10-27-2018
%load bmd_fullsweep_vf
clear d_k
clear d_E_hole
clear d_E_elec

for theta_idx = 1:length(tar_theta_list)
    kpts = sweep_kpts1{theta_idx};
    bands = sweep_bands1{theta_idx};
    nb = size(bands,1);
    
    start_idx = nb/2;
    for trial_idx = start_idx-1:start_idx+1
        if ( (bands(trial_idx+1,1) - bands(trial_idx)) < 1e-7 )
           start_idx = trial_idx;
           break
        end
    end
    
    d_k(theta_idx) = norm(kpts(end,:) - kpts(end-dk_idx,:));
    d_E_hole(theta_idx) = bands(start_idx,end) - bands(start_idx,end-dk_idx);
    d_E_elec(theta_idx) = bands(start_idx+1,end) - bands(start_idx+1,end-dk_idx);
    
    bands2 = sweep_bands2{theta_idx};
    
    n4_elec = max([bands(start_idx+1,:) bands2(start_idx+1,:)]);
    min_elec = min([bands(start_idx+2,:) bands2(start_idx+2,:)]);
    n4_hole = min([bands(start_idx,:) bands2(start_idx,:)]);
    max_hole = max([bands(start_idx-1,:) bands2(start_idx-1,:)]);
    
    gap_elec(theta_idx) = min_elec - n4_elec;
    gap_hole(theta_idx) = n4_hole - max_hole;
    
    disp_elec(theta_idx) = range(bands(start_idx+1,:));
    disp_hole(theta_idx) = range(bands(start_idx,:));
end

vf_elec_bmd = d_E_elec./(hbar*d_k)*10^-10; % angstroms to meters
vf_hole_bmd = d_E_hole./(hbar*d_k)*10^-10;

l_theta = 1.43*sqrt(3)./(2*sind(tar_theta_list/2));
k_theta = (2/3)*2*pi./l_theta;
alpha_theta = (0.11/Dirac_v)./k_theta;

%alpha_axis = (theta_0./theta_axis)*alpha_0;
alpha_axis = alpha_theta;
alpha2_axis = alpha_axis.^2;

clf
%subaxis(2,1,1,'sh', 0.01, 'sv', 0.05, 'padding', 0, 'margintop', .16, 'marginbot',.15, 'marginleft', .125);
hold on

plot([0 100],offsets(1)+[0 0],'-','Color',bmd_color,'HandleVisibility','off')
plot([0 100],offsets(2)+[0 0],'-','Color',norelax_color,'HandleVisibility','off')

bmd_elec_set = [1:90, 94:length(vf_elec_bmd)];

plot(alpha2_axis(bmd_elec_set),abs(vf_elec_bmd(bmd_elec_set))/v0 + offsets(1),elec_type,'Color',bmd_color,'LineWidth',2)
%plot(alpha2_axis(bmd_elec_set),abs(vf_elec_bmd(bmd_elec_set))/v0 + offsets(1),elec_type,'Color',bmd_color,'LineWidth',2)
%plot(alpha2_norelax,abs(vf_hole_bmd),hole_type,'Color',bmd_color,'LineWidth',2)

noelec_set = [1:88,92:123,125:143, 146:length(vf_elec_norelax)];
plot(alpha2_norelax(noelec_set),abs(vf_elec_norelax(noelec_set))/v0 + offsets(2),elec_type,'Color',norelax_color,'LineWidth',2)
%plot(alpha2_norelax,abs(vf_hole_norelax),hole_type,'Color',norelax_color,'LineWidth',2)
%plot(alpha2_norelax,abs(vf_hole_norelax)/v0 + offsets(2),hole_type,'Color',norelax_color,'LineWidth',2)

%relax_split = 83;

%plot(alpha2_relax(1:relax_split),abs(vf_elec_relax(1:relax_split))/v0 + offsets(3),elec_type,'Color',relax_color,'LineWidth',2)
%plot(alpha2_relax(relax_split+1:end),abs(vf_elec_relax(relax_split+1:end))/v0 + offsets(3),'--','Color',relax_color,'LineWidth',2)
plot(alpha2_relax,abs(vf_elec_relax)/v0 + offsets(3),elec_type,'Color',relax_color,'LineWidth',2)

%plot(alpha2_relax,abs(vf_hole_relax),hole_type,'Color',relax_color,'LineWidth',2)
box on

theta_tick_labels = [1.00 .5 .35 .26 .21];
l_theta = 1.43*sqrt(3)./(2*sind(theta_tick_labels/2));
k_theta = (2/3)*2*pi./l_theta;
alpha_theta = (0.11/Dirac_v)./k_theta;
alpha2_tick_labels = alpha_theta.^2;

axis([0 alpha_max 0 vf_max])
set(gca,'xaxislocation','top')
set(gca,'box','off')
ylabel('$v_F/v_0$')
set(gca,'XTick',alpha2_tick_labels)
theta_tick_list = {['$' num2str(theta_tick_labels(1),'%.2f') '^\circ$']
                   ['$' num2str(theta_tick_labels(2),'%.2f') '^\circ$']
                   ['$' num2str(theta_tick_labels(3),'%.2f') '^\circ$']
                   ['$' num2str(theta_tick_labels(4),'%.2f') '^\circ$']
                   ['$' num2str(theta_tick_labels(5),'%.2f') '^\circ$']

                        };
    
xticklabels(theta_tick_list)
xlabel('$\theta$')

%{
lgd = legend('BMD','Unrelaxed','Relaxed');
lgd.FontSize = 11;
lgd.Position = [0.732 0.68 .15 .1];
%}

text(alpha_max-.1,offsets(1)+.075,'BMD','color',bmd_color,'FontSize',f_size,'HorizontalAlignment','right')
text(alpha_max-.1,offsets(2)+.100,'Unrelaxed','color',norelax_color,'FontSize',f_size,'HorizontalAlignment','right')
text(alpha_max-.1,offsets(3)+.050,'Full Relax','color',relax_color,'FontSize',f_size,'HorizontalAlignment','right')

%text(.3*alpha2_axis(end),.9*vf_max,'electron','FontSize',f_size)

set(gca,'YTick',[0 .8])
set(gca,'TickDir','out')
text(-.05,offsets(1),'0','Color',bmd_color,'FontSize',f_size,'HorizontalAlignment','right')
text(-.05,offsets(2),'0','Color',norelax_color,'FontSize',f_size,'HorizontalAlignment','right')
set(gca,'pos',[.125 .155 .8 .675])

ax1 = gca;
ax1_pos = ax1.Position; % position of first axes
ax2 = axes('Position',ax1_pos,...
    'XAxisLocation','bottom',...
    'YAxisLocation','right',...
    'Color','none');
set(gca,'YTick',[])
set(gca,'TickDir','out')
axis([0 alpha_max 0 vf_max])
set(gca,'Xtick',[0:2:12])
xlabel('$\alpha^2$')

%{
subaxis(2,1,2)
hold on

plot([0 100],offsets(1)+[0 0],'-','Color',bmd_color)
plot([0 100],offsets(2)+[0 0],'-','Color',norelax_color)

%bmd_hole_set = [1:42, 45:113, 117:120];

%plot(alpha2_norelax,abs(vf_elec_bmd),elec_type,'Color',bmd_color,'LineWidth',2)
%plot(alpha2_axis(bmd_hole_set),abs(vf_hole_bmd(bmd_hole_set))/v0 + offsets(1),hole_type,'Color',bmd_color,'LineWidth',2)
plot(alpha2_axis,abs(vf_hole_bmd)/v0 + offsets(1),hole_type,'Color',bmd_color,'LineWidth',2)

%plot(alpha2_norelax,abs(vf_elec_norelax),elec_type,'Color',norelax_color,'LineWidth',2)
plot(alpha2_norelax,abs(vf_hole_norelax)/v0 + offsets(2),hole_type,'Color',norelax_color,'LineWidth',2)

%plot(alpha2_relax,abs(vf_elec_relax),elec_type,'Color',relax_color,'LineWidth',2)
plot(alpha2_relax,abs(vf_hole_relax)/v0 + offsets(3),hole_type,'Color',relax_color,'LineWidth',2)


text(.3*alpha2_axis(end),.9*vf_max,'hole','FontSize',f_size)
set(gca,'YTick',[0 .8])
text(-.05,offsets(1),'0','Color',bmd_color,'FontSize',f_size,'HorizontalAlignment','right')
text(-.05,offsets(2),'0','Color',norelax_color,'FontSize',f_size,'HorizontalAlignment','right')

axis([0 alpha_max 0 vf_max])
set(gca,'Xtick',[0:.5:12])
ylabel('$v_F/v_0$')
xlabel('$\alpha^2$')
%}





