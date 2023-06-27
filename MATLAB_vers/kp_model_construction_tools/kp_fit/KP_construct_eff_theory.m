%clear all;

% here the length scale adopts a=2.46A = 1

sc_m=60;
sc_n=sc_m-1;

[sc_m,sc_n]

%load(['TwBLG_EffKP_relax_sym_BZscan_',num2str(sc_m),'_',num2str(sc_n),'_no']);

% savename=['TwBLG_EffKP_symparms1_',num2str(sc_m),'_',num2str(sc_n),'_no'];

% note the sc_bi are in units of 1/lattice_a

% data needed from TBH construction:
%{
sc_b1
sc_b2
pc_vec_b1
pc_vec_b2
pc_vec_rb1
pc_vec_rb2
BZ_allham
allkxy
sc_all_points_shift
bot_B_atoms
bot_A_atoms
top_B_atoms
top_A_atoms
%}
%%
save( ['reduced_dat_BZscan_',num2str(sc_m),'_',num2str(sc_n)], ...
        'sc_b1','sc_b2','pc_vec_b1','pc_vec_b2','pc_vec_rb1','pc_vec_rb2', ...
        'BZ_allham','allkxy','sc_all_points_shift','tot_num', ...
        'bot_A_atoms','bot_B_atoms','top_A_atoms','top_B_atoms');
    
%%

clear all;

% here the length scale adopts a=2.46A = 1

sc_m=30;
sc_n=sc_m-1;

[sc_m,sc_n]

load(['reduced_dat_BZscan_',num2str(sc_m),'_',num2str(sc_n)]);

    
%% constructs k-p degrees of freedom


hex_cut=2.51;
hex_cut=hex_cut*sqrt(dot(sc_b1,sc_b1));

hex_shift=(-sc_b1+sc_b2)/3;

hex_M=40;

hex_table=zeros(2*hex_M+1);

hex_index=0;
hex_coor=0;

ind=1;
for ind1=(-hex_M):hex_M
    for ind2=(-hex_M):hex_M
        vec=sc_b1*ind1+sc_b2*ind2+hex_shift;
        
        if sqrt(dot(vec,vec))<hex_cut
            %hex_index(ind,1:2)=[ind1,ind2];
            hex_table(ind1+hex_M+1,ind2+hex_M+1)=ind;
            hex_coor(ind,1:2)=vec(1:2);
            
            ind=ind+1;
        end
        
    end
end

num_hex=ind-1;


% find all sym qs

bot_K_point=(-pc_vec_b1+pc_vec_b2)/3
top_K_point=(-pc_vec_rb1+pc_vec_rb2)/3
hex_shift=(-sc_b1+sc_b2)/3;
%
% all_qs(1,:)=hex_shift;
% all_qs(2,:)=hex_shift+sc_b2;
% all_qs(3,:)=hex_shift-sc_b2;
% all_qs(4,:)=hex_shift-2*sc_b2;
% all_qs(5,:)=hex_shift-sc_b1;
% all_qs(6,:)=hex_shift-sc_b1-sc_b2;
% all_qs(7,:)=hex_shift-sc_b1-2*sc_b2;
% all_qs(8,:)=hex_shift+sc_b1;
% all_qs(9,:)=hex_shift+sc_b1+sc_b2;
% all_qs(10,:)=hex_shift+sc_b1-sc_b2;
% all_qs(11,:)=hex_shift+2*sc_b1;
% all_qs(12,:)=hex_shift+2*sc_b1+sc_b2;

% create momenta values for bottom and top layers by shifting
% the sampling by the momenta of the respective K points

all_bot_qs=hex_coor;
all_bot_qs(:,1)=all_bot_qs(:,1)+bot_K_point(1);
all_bot_qs(:,2)=all_bot_qs(:,2)+bot_K_point(2);
all_top_qs=-hex_coor;
all_top_qs(:,1)=all_top_qs(:,1)+top_K_point(1);
all_top_qs(:,2)=all_top_qs(:,2)+top_K_point(2);


% supercell reciprical lattice
qqmat=[sc_b1(1),sc_b2(1);sc_b1(2),sc_b2(2)];
invqqmat=inv(qqmat);

% checks if the new momenta lie on the lattice?
for inds=1:size(all_bot_qs,1)
    qqnow=all_bot_qs(inds,1:2)';
    
    qqconv=invqqmat*qqnow;
    
    qqconv
    
    
    qqnow=all_top_qs(inds,1:2)';
    
    qqconv2=invqqmat*qqnow;
    
    qqconv2
end


%% intra layer terms

% sample intralyer terms over the hex samples
% this is done for every unique k-point for which we have a Hamiltonian


bot_intra_list=zeros(2,2,num_hex*length(BZ_allham));
bot_all_intra_kk=zeros(num_hex*length(BZ_allham),2);

top_intra_list=zeros(2,2,num_hex*length(BZ_allham));
top_all_intra_kk=zeros(num_hex*length(BZ_allham),2);

indcc=1;
for indk=1:length(BZ_allham)
    % get the Hamiltonian at this k point
    Hmatnow=BZ_allham{indk};
    % get the k point
    k_scnow=allkxy(indk,1:2);
    
    for indh=1:num_hex
        
        % get k points of each layer
        kkbot=all_bot_qs(indh,1:2);
        kktop=all_top_qs(indh,1:2);
        
        
        bkpt_test=kkbot;
        tkpt_test=kktop;
        
        % construct bottom wavefuction
        wf_set1_b=zeros(tot_num,2);
        
        tmp_k=bkpt_test(1:2);
        tmp_k=reshape(tmp_k,2,1);
        tmp_pos=sc_all_points_shift(bot_A_atoms,:);
        wf_set1_b(bot_A_atoms,1)=exp(i*tmp_pos*tmp_k);
        tmp_pos=sc_all_points_shift(bot_B_atoms,:);
        wf_set1_b(bot_B_atoms,2)=exp(i*tmp_pos*tmp_k);
        wf_set1_b=wf_set1_b/sqrt(tot_num/4);
        
        %  construct top wavefunction
        wf_set1_t=zeros(tot_num,2);
        
        tmp_k=tkpt_test(1:2);
        tmp_k=reshape(tmp_k,2,1);
        tmp_pos=sc_all_points_shift(top_A_atoms,:);
        wf_set1_t(top_A_atoms,1)=exp(i*tmp_pos*tmp_k);
        tmp_pos=sc_all_points_shift(top_B_atoms,:);
        wf_set1_t(top_B_atoms,2)=exp(i*tmp_pos*tmp_k);
        wf_set1_t=wf_set1_t/sqrt(tot_num/4);
        
        % get the matrix elements
        Eff_intra_bot=wf_set1_b'*Hmatnow*wf_set1_b;
        Eff_intra_top=wf_set1_t'*Hmatnow*wf_set1_t;
        
        % save matrix elements to list
        bot_intra_list(:,:,indcc)=Eff_intra_bot(:,:);
        top_intra_list(:,:,indcc)=Eff_intra_top(:,:);
        
        % get relative momenta of this hex k-point
        kkbot_rel=hex_coor(indh,1:2);
        kktop_rel=-hex_coor(indh,1:2);
        
        % save global momenta to list
        bot_all_intra_kk(indcc,:)=kkbot_rel(1:2)+k_scnow(1:2);
        top_all_intra_kk(indcc,:)=kktop_rel(1:2)+k_scnow(1:2);
        
        
        indcc=indcc+1;
    end
end

%% Check the intralyer coupling values (should be top half of a Dirac cone!)
clf
hold on
scatter3(bot_all_intra_kk(:,1),bot_all_intra_kk(:,2),abs(bot_intra_list(1,2,:)))
%scatter3(top_all_intra_kk(:,1),top_all_intra_kk(:,2),abs(top_intra_list(1,2,:)))
view(3)
%%

close all
clear splot;

figure('Position', [100 100 800 800]);
f_size = 16;
set(groot, 'DefaultTextInterpreter', 'Latex')
set(groot, 'DefaultLegendInterpreter', 'Latex')
set(groot, 'DefaultAxesTickLabelInterpreter', 'Latex')
set(0,'DefaultAxesFontSize',f_size)

abs_cmap = colorcet('L3');
phase_cmap = colorcet('C2');

max_grid_k = .2;
grid_dk = 0.002;


k_grid = -max_grid_k:grid_dk:max_grid_k;
[kx_grid, ky_grid] = meshgrid(k_grid,k_grid);
H_ab_interp = griddata(squeeze(bot_all_intra_kk(:,1)),squeeze(bot_all_intra_kk(:,2)),squeeze(bot_intra_list(1,2,:)),kx_grid,ky_grid);
H_ba_interp = griddata(squeeze(bot_all_intra_kk(:,1)),squeeze(bot_all_intra_kk(:,2)),squeeze(bot_intra_list(2,1,:)),kx_grid,ky_grid);
plot_data1 = H_ab_interp;
plot_data2 = H_ba_interp;

clf

splot(1) = subplot(2,2,1);
surf(kx_grid,ky_grid,abs(plot_data1),'EdgeColor','none')
view(2)
axis equal
colormap(gca,abs_cmap);
title('$|H_{AB}(k)|$')
ylabel('$k_y$')

splot(2) = subplot(2,2,2);
surf(kx_grid,ky_grid,angle(plot_data1),'EdgeColor','none')
view(2)
axis equal
colormap(gca,phase_cmap);
xlabel('$k_x$')
ylabel('$k_y$')
%title('$H_{AB}$ phase')

splot(3) = subplot(2,2,3);
surf(kx_grid,ky_grid,abs(plot_data2),'EdgeColor','none')
view(2)
axis equal
colormap(gca,abs_cmap);
title('$|H_{BA}(k)|$')
cbar_abs = colorbar(gca);

splot(4) = subplot(2,2,4);
surf(kx_grid,ky_grid,mod(angle(plot_data2),2*pi),'EdgeColor','none')
view(2)
axis equal
colormap(gca,phase_cmap);
%title('$H_{BA}$ phase')
xlabel('$k_x$')
cbar_phase = colorbar(gca);

ax_w = .375;
ax_sx = .02;
ax_sy = .05;
b_s = .1;

set(splot(1),'pos',[b_s             ax_w+ax_sy+b_s    ax_w ax_w])
set(splot(3),'pos',[ax_w+ax_sx+b_s   ax_w+ax_sy+b_s   ax_w ax_w])
set(splot(2),'pos',[b_s             b_s              ax_w ax_w])
set(splot(4),'pos',[ax_w+ax_sx+b_s   b_s             ax_w ax_w])
set(cbar_phase,'YTick',[0.0001:pi/2:3.5*pi/2, 2*pi-.0001])
set(cbar_phase,'yticklabel',{'$0^\circ$','$90^\circ$','$180^\circ$','$270^\circ$','$360^\circ$'})
cbar_phase.TickLabelInterpreter = 'latex';
cbar_abs.TickLabelInterpreter = 'latex';

%%

% data_kx=top_all_intra_kk(:,1);
% data_ky=top_all_intra_kk(:,2);
% data_fit=squeeze(top_intra_list(2,1,:));


data_kx=bot_all_intra_kk(:,1);
data_ky=bot_all_intra_kk(:,2);
data_fit=squeeze(bot_intra_list(2,1,:));

data_kplus=data_kx+1i*data_ky;
data_kminus=data_kx-1i*data_ky;

% CC is a 3x3 matrix which represents:
% sum[<(k+, k-, 1)| (x) |(k+, k-, 1)>]
% over every unique k-point

CC_mat=zeros(3);
CC_mat(1,1)=sum(conj(data_kplus).*data_kplus);
CC_mat(1,2)=sum(conj(data_kplus).*data_kminus);
CC_mat(1,3)=sum(conj(data_kplus));

CC_mat(2,1)=sum(conj(data_kminus).*data_kplus);
CC_mat(2,2)=sum(conj(data_kminus).*data_kminus);
CC_mat(2,3)=sum(conj(data_kminus));

CC_mat(3,1)=sum(data_kplus);
CC_mat(3,2)=sum(data_kminus);
CC_mat(3,3)=length(data_fit);

% DDvec is a 3x1 vector which represents:
% sum[H_ab(intra)|k+, k-, 1>] over every kpoint

DDvec(1,1)=sum(conj(data_kplus).*data_fit);
DDvec(2,1)=sum(conj(data_kminus).*data_fit);
DDvec(3,1)=sum(data_fit);


% CC optimal gives the prefactors for the k+, k-, and identity terms of the
% intralayer Hamiltonian
CC_optimal=inv(CC_mat)*DDvec;

CC_optimal



% top (1,2):
%   -0.0448 - 2.1198i
%    0.0415 - 0.0225i
%   -0.0010 + 0.0020i

% top (2,1)
%    0.0415 + 0.0225i
%   -0.0448 + 2.1198i
%   -0.0010 - 0.0020i

% bot (1,2)
%    0.0033 - 2.1203i
%    0.0407 - 0.0239i
%   -0.0013 + 0.0019i

% bot (2,1)
%    0.0407 + 0.0239i
%    0.0033 + 2.1203i
%   -0.0013 - 0.0019i
  

%%

% Expand the original data in terms of k+, k-, and 1
data_optimal=CC_optimal(1)*data_kplus+CC_optimal(2)*data_kminus+CC_optimal(3);

% Plot the data: fit in black, error in red
clf
hold on
scatter3(data_kx,data_ky,abs(data_fit-data_optimal),'r')
scatter3(data_kx,data_ky,abs(data_optimal),'k')
view(3)

%%


q_resolution=1E-6;

% inter plane terms

q_resolution=1E-6;

interall_given_qs(1,:)=(-2*hex_shift+sc_b2);
interall_given_qs(2,:)=-hex_shift-(hex_shift+sc_b1);
interall_given_qs(3,:)=(-hex_shift-sc_b1)-(hex_shift-sc_b2);
%
% interall_given_qs(4,:)=-2*hex_shift;
% interall_given_qs(5,:)=-2*hex_shift+2*sc_b2;
% interall_given_qs(6,:)=-2*hex_shift-2*sc_b1;
%
% interall_given_qs(7,:)=(-2*hex_shift+sc_b2)+sc_b1;
% interall_given_qs(8,:)=(-2*hex_shift+sc_b2)+sc_b1+sc_b2;
%
% interall_given_qs(9,:)=(-2*hex_shift+sc_b2)-sc_b1+sc_b2;
% interall_given_qs(10,:)=(-2*hex_shift+sc_b2)-2*sc_b1;
%
% interall_given_qs(11,:)=(-2*hex_shift+sc_b2)-(sc_b1+sc_b2)-sc_b2;
% interall_given_qs(12,:)=(-2*hex_shift+sc_b2)-(sc_b1+sc_b2)*2;

numq_inter=size(interall_given_qs);
numq_inter=numq_inter(1);
numqs=12;

q_resolution=1E-6;

All_inter_list={};

for indqq=1:numq_inter
    given_q=interall_given_qs(indqq,1:2);
    
    clear Eff_inter;
    
    
    tmp_list=[];
    indc=0;
    for indq1=1:num_hex
        qtmptop=-hex_coor(indq1,:);
        for indq2=1:num_hex
            qtmpbot=hex_coor(indq2,:);
            
            qtmpdiff=qtmptop-qtmpbot-given_q;
            if sqrt(dot(qtmpdiff,qtmpdiff))<q_resolution
                indc=indc+1;
                
                tmp_list(indc,1:2)=[indq1,indq2]; %,qtmptop(1)-qtmpbot(1),qtmptop(2)-qtmpbot(2)];
                
                
            end
            
            
        end
    end
    indc
    All_inter_list{indqq}=tmp_list;
    
end


%%

indqq=3;

list_now=All_inter_list{indqq};
sizetmp=size(list_now,1);

indcc=1;

collect_list=zeros(2,2,sizetmp*length(BZ_allham));
all_collect_kk=zeros(sizetmp*length(BZ_allham),2);

for indk=1:length(BZ_allham)
    
    indk/length(BZ_allham)
    %{
    if (indk == 1)
        
    else
       for i = 1:num_backspace-1
          fprintf('\b')
       end
    end
    
    percent_string = ['Completed: ' num2str(100*indk/length(BZ_allham),'%1.f') '%%'];
    fprintf(percent_string);
    %pause(.01)
    num_backspace = length(percent_string);
    %}
    
    Hmatnow=BZ_allham{indk};
    k_scnow=allkxy(indk,1:2);
    
    for indl=1:sizetmp
        tmpnn=list_now(indl,1:2);
        
        kkbot=all_bot_qs(tmpnn(2),1:2);
        kktop=all_top_qs(tmpnn(1),1:2);
        
        
        bkpt_test=kkbot;
        tkpt_test=kktop;
        
        wf_set1_b=zeros(tot_num,2);
        
        tmp_k=bkpt_test(1:2);
        tmp_k=reshape(tmp_k,2,1);
        tmp_pos=sc_all_points_shift(bot_A_atoms,:);
        wf_set1_b(bot_A_atoms,1)=exp(i*tmp_pos*tmp_k);
        tmp_pos=sc_all_points_shift(bot_B_atoms,:);
        wf_set1_b(bot_B_atoms,2)=exp(i*tmp_pos*tmp_k);
        wf_set1_b=wf_set1_b/sqrt(tot_num/4);
        
        wf_set1_t=zeros(tot_num,2);
        
        tmp_k=tkpt_test(1:2);
        tmp_k=reshape(tmp_k,2,1);
        tmp_pos=sc_all_points_shift(top_A_atoms,:);
        wf_set1_t(top_A_atoms,1)=exp(i*tmp_pos*tmp_k);
        tmp_pos=sc_all_points_shift(top_B_atoms,:);
        wf_set1_t(top_B_atoms,2)=exp(i*tmp_pos*tmp_k);
        wf_set1_t=wf_set1_t/sqrt(tot_num/4);
        
        Eff_inter_tmp=wf_set1_t'*Hmatnow*wf_set1_b;
        
        
        kkbot_rel=hex_coor(tmpnn(2),1:2);
        
        
        all_collect_kk(indcc,:)=kkbot_rel(1:2)+k_scnow(1:2);
        collect_list(:,:,indcc)=Eff_inter_tmp(:,:);
        
        indcc=indcc+1;
    end
    
    
    
    
    
    
end

All_CCopt=zeros(3,2,2);

for ind1=1:2
    for ind2=1:2
        
        sel_data=squeeze(collect_list(ind1,ind2,:));
        
        % k_\pm = kx \pm iky
        kk_plus=all_collect_kk(:,1)+i*all_collect_kk(:,2);
        kk_minus=all_collect_kk(:,1)-i*all_collect_kk(:,2);
        
        
        CC_mat=zeros(3);
        CC_mat(1,1)=sum(conj(kk_plus).*kk_plus);
        CC_mat(1,2)=sum(conj(kk_plus).*kk_minus);
        CC_mat(1,3)=sum(conj(kk_plus));
        
        CC_mat(2,1)=sum(conj(kk_minus).*kk_plus);
        CC_mat(2,2)=sum(conj(kk_minus).*kk_minus);
        CC_mat(2,3)=sum(conj(kk_minus));
        
        CC_mat(3,1)=sum(kk_plus);
        CC_mat(3,2)=sum(kk_minus);
        CC_mat(3,3)=length(kk_minus);
        
        DDvec=zeros(3,1);
        DDvec(1)=sum(conj(kk_plus).*sel_data);
        DDvec(2)=sum(conj(kk_minus).*sel_data);
        DDvec(3)=sum(sel_data);
        
        CCoptimal=inv(CC_mat)*DDvec;
        All_CCopt(:,ind1,ind2)=CCoptimal(:);
        
    end
    
    
end


All_CCopt(:,1,1)
All_CCopt(:,1,2)
All_CCopt(:,2,1)
All_CCopt(:,2,2)

% const term 

'done'

%%
clf
scatter3(all_collect_kk(:,1),all_collect_kk(:,2),abs(collect_list(1,2,:)),'k')


%%
close all
clear splot;

figure('Position', [100 100 800 800]);
f_size = 16;
set(groot, 'DefaultTextInterpreter', 'Latex')
set(groot, 'DefaultLegendInterpreter', 'Latex')
set(groot, 'DefaultAxesTickLabelInterpreter', 'Latex')
set(0,'DefaultAxesFontSize',f_size)

abs_cmap = colorcet('L3');
phase_cmap = colorcet('C2');

max_grid_k = .2;
grid_dk = 0.002;

raw_data1 = squeeze(collect_list(1,1,:));
raw_data2 = squeeze(collect_list(1,2,:));
raw_k_list = all_collect_kk;

k_grid = -max_grid_k:grid_dk:max_grid_k;
[kx_grid, ky_grid] = meshgrid(k_grid,k_grid);
H_ab_interp = griddata(squeeze(raw_k_list(:,1)),squeeze(raw_k_list(:,2)),raw_data1,kx_grid,ky_grid);
H_ba_interp = griddata(squeeze(raw_k_list(:,1)),squeeze(raw_k_list(:,2)),raw_data2,kx_grid,ky_grid);
plot_data1 = H_ab_interp;
plot_data2 = H_ba_interp;

clf

phase_max = [180]*pi/180;
amp_max = .15;


splot(1) = subplot(2,2,1);
surf(kx_grid,ky_grid,abs(plot_data1),'EdgeColor','none')
view(2)
axis equal
colormap(gca,abs_cmap);
title('$|H_{AB}(k)|$')
ylabel('$k_y$')
caxis([0 amp_max])

splot(2) = subplot(2,2,2);
surf(kx_grid,ky_grid,angle(plot_data1),'EdgeColor','none')
view(2)
axis equal
colormap(gca,phase_cmap);
xlabel('$k_x$')
ylabel('$k_y$')
%title('$H_{AB}$ phase')
caxis([-phase_max phase_max])

splot(3) = subplot(2,2,3);
surf(kx_grid,ky_grid,abs(plot_data2),'EdgeColor','none')
view(2)
axis equal
colormap(gca,abs_cmap);
title('$|H_{BA}(k)|$')
cbar_abs = colorbar(gca);
caxis([0 amp_max])

splot(4) = subplot(2,2,4);
surf(kx_grid,ky_grid,angle(plot_data2),'EdgeColor','none')
view(2)
axis equal
colormap(gca,phase_cmap);
%title('$H_{BA}$ phase')
xlabel('$k_x$')
cbar_phase = colorbar(gca);
caxis([-phase_max phase_max])

ax_w = .375;
ax_sx = .02;
ax_sy = .05;
b_s = .1;

set(splot(1),'pos',[b_s             ax_w+ax_sy+b_s    ax_w ax_w])
set(splot(3),'pos',[ax_w+ax_sx+b_s   ax_w+ax_sy+b_s   ax_w ax_w])
set(splot(2),'pos',[b_s             b_s              ax_w ax_w])
set(splot(4),'pos',[ax_w+ax_sx+b_s   b_s             ax_w ax_w])
%set(cbar_phase,'YTick',[0.0001:pi/2:3.5*pi/2, 2*pi-.0001])
%set(cbar_phase,'yticklabel',{'$0^\circ$','$90^\circ$','$180^\circ$','$270^\circ$','$360^\circ$'})
cbar_phase.TickLabelInterpreter = 'latex';
cbar_abs.TickLabelInterpreter = 'latex';

%%

interall_given_qs(indqq,:)
All_CCopt_const=squeeze(All_CCopt(3,:,:));
All_CCopt_const
% angle(All_CCopt_const)/pi*180

% q=1
%
%    0.1116 + 0.0000i   0.1116 + 0.0000i
%    0.1116 + 0.0000i   0.1116 - 0.0000i

% q=2
%    0.1117 + 0.0000i  -0.0559 - 0.0967i
%   -0.0559 + 0.0967i   0.1117 - 0.0000i

% q=3
%    0.1115 + 0.0000i  -0.0557 + 0.0966i
%   -0.0557 - 0.0966i   0.1115 - 0.0000i


All_CCopt_kplus=squeeze(All_CCopt(1,:,:));
All_CCopt_kplus
% q=1
%   -0.0002 + 0.0451i  -0.0002 + 0.0451i
%   -0.0002 + 0.0451i  -0.0002 + 0.0451i

% q=2
%   -0.0395 - 0.0216i   0.0010 + 0.0450i
%    0.0385 - 0.0234i  -0.0395 - 0.0216i

% q=3
%    0.0413 - 0.0244i   0.0005 + 0.0480i
%   -0.0418 - 0.0236i   0.0413 - 0.0244i



All_CCopt_kminus=squeeze(All_CCopt(2,:,:));
All_CCopt_kminus
% q=1
%   -0.0002 - 0.0451i  -0.0002 - 0.0451i
%   -0.0002 - 0.0451i  -0.0002 - 0.0451i

% q=2
%   -0.0395 + 0.0216i   0.0385 + 0.0234i
%    0.0010 - 0.0450i  -0.0395 + 0.0216i

% q=3
%    0.0413 + 0.0244i  -0.0418 + 0.0236i
%    0.0005 - 0.0480i   0.0413 + 0.0244i


www=exp(i*2*pi/3);














