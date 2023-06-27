%clear all
%load dos_sweep_0p5_01-09-2019.mat;

f_size = 20;

clf
hold on

theta_list = [0.52, 0.48, 0.44];

plot((idos_sweep{1}-0.0)*1.0,dos_sweep{1},'r')
plot((idos_sweep{2}-0.0)*1.0,dos_sweep{2},'k')
plot((idos_sweep{3}-0.0)*1.0,dos_sweep{3},'b')

%plot(idos{1},dos_sweep{1},'r')
%plot(idos{2},dos_sweep{2},'k')
%plot(idos{3},dos_sweep{3},'b')

dos_max_0p5_r = 4;
axis([-10 10 0 dos_max_0p5_r])

%axis([7.9 8.1 0 0.015])
%axis([-8.05 -7.95 0 0.015])

for x = -8:4:8
   plot(x+[0 0],[-10 10],'--k') 
end
set(gca,'XTick',[-28:4:28])
xlabel('$n/n_0$')
ylabel('DoS $(eV^{-1} nm^{-2})$')

% title
text(-9.8,dos_max_0p5_r*.9,['relaxed'],'Color','k','FontSize',f_size);

% legend
text(8.5,dos_max_0p5_r*.85,'$0.52^\circ$','Color','r','FontSize',f_size-4);
text(8.5,dos_max_0p5_r*.7,'$0.48^\circ$','Color','k','FontSize',f_size-4);
text(8.5,dos_max_0p5_r*.55,'$0.44^\circ$','Color','b','FontSize',f_size-4);


%%
dos_max = 60;
m = 10;

idos_h = idos{1};
theta_here = theta_list(1);
dos_h = dos_sweep{1};

b_size = 3;
tot_bands = 2*(b_size+1);
tot_bands = 4*tot_bands; % 2 for valley, 2 for spin
idos_rescale = tot_bands/idos_h(end);
alpha = 2.47;
sc_alpha = alpha/(2*sind(theta_here/2));
sc_area = sc_alpha^2*sind(60)*1e-2; %area in nm^2
n0 = 1/sc_area;
dos_rescale = idos_rescale*n0;

clf
plot(idos_h,dos_rescale*dos_h,'k','LineWidth',2);
%text(-14,dos_max*.9,'unrelaxed 0.48^\circ')

axis([-m m 0 dos_max])
set(gca,'XTick',[-m+rem(m,4):4:m])
%xticklabels({})
set(gca,'YTick',[])
ylabel('DoS')
xlabel('n/n0')
%%

data_dos_r_0p52 = dos_sweep{1};
data_idos_r_0p52 = idos_sweep{1};

data_dos_r_0p48 = dos_sweep{2};
data_idos_r_0p48 = idos_sweep{2};

data_dos_r_0p44 = dos_sweep{3};
data_idos_r_0p44 = idos_sweep{3};

save('hyobin_paper_data_extended_relaxed_dos_v2.mat',...
    'data_dos_r_0p52','data_idos_r_0p52',...
    'data_dos_r_0p48','data_idos_r_0p48',...
    'data_dos_r_0p44','data_idos_r_0p44')

    

