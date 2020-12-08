%% Get full BZ data (LDOS)
clear all;

theta_list = 2.0;

% realspace sampling
nsamps = 20;
dx = linspace(0,1,nsamps+1);
dx = dx(1:end-1);
[x,y] = meshgrid(dx,dx);

ldos_pos = [x(:),y(:)];

tic
[sweep_vals, scaleaxis, sweep_kpts, sweep_weights] = tblg_kp_calc_ext('theta_list',theta_list,'knum',20,'full_bz',1,'ldos_positions',ldos_pos);

[dos_sweep, ldos_sweep, E_list] = interp_kp_ldos(theta_list, sweep_vals, sweep_weights, sweep_kpts);
toc

%% check the weights (used for debugging, can skip)

wh = sweep_weights{1};
tar_k = 1;
tar_b = size(wh,3)/2;
tar_o = 4;
w_r = squeeze(mean(wh(:,:,tar_b,tar_k),1));
sum(w_r(:))

clf
subplot(2,1,1)
plot(w_r)
subplot(2,1,2)
A = [sqrt(3)/2 sqrt(3)/2; -0.5 0.5];
r_sc = A(:,1).*x(:)' + A(:,2).*y(:)';
rx = reshape(r_sc(1,:),nsamps,nsamps);
ry = reshape(r_sc(2,:),nsamps,nsamps);
surf(rx,ry,reshape(w_r,nsamps,nsamps))
axis equal
view(2)
colorbar

%% nice plot of DOS, energy selection, and real-space LDOS
ldos = ldos_sweep{1};
dos = dos_sweep{1};

tar_E = 1250;

clf
subplot(2,1,1)
title(['DOS for $' num2str(theta_list(1)) '^\circ$ TBG'])
hold on
box on
plot(E_list,dos,'-k');
%plot(mean(ldos,2),'--r')
plot(E_list(tar_E)+[0 0],[0 max(dos)*1.2],'-r')
axis([-inf inf 0 max(dos)*1.2])
xlabel('Energy (eV)')
ylabel('DoS')

subplot(2,1,2)
hold on
box on
%{
titles = {'AA','AB','BA','DW1','DW2'};
%%{
sscl = max(ldos(:)); % plot splitting scale
nr = size(ldos,2);
plot(E_list(tar_E)+[0 0],999*[-1 1],'--k')

pidx = 1;
for idx = [1,100,200]
    plot(E_list,abs(ldos(:,idx))+sscl*(pidx-1),'LineWidth',.5)
    %text(E_list(100),sscl*(pidx-1)+abs(ldos(1,idx))+.01,titles{idx});
    plot([E_list(1) E_list(end)],sscl*(pidx-1)+[0 0],'--k')
    pidx = pidx+1;
end
%
%plot(E_list, sum(abs(ldos),2));

axis([E_list(1) E_list(end) 0 3*sscl])
%}

% LDOS surface at given energy
A = [sqrt(3)/2 sqrt(3)/2; -0.5 0.5];
r_sc = A(:,1).*x(:)' + A(:,2).*y(:)';
nr = sqrt(size(ldos,2));

rx = reshape(r_sc(1,:),nr,nr);
ry = reshape(r_sc(2,:),nr,nr);

scg = 1; % supercell grid scaling
nrsc = (2*scg+1)*nr;
rscx = zeros(nrsc,nrsc);
rscy = rscx;
ldos_sc = rscx;

for scx = -scg:scg
    for scy = -scg:scg
        tar_idx_x = [1:nr]+nr*(scx+scg);
        tar_idx_y = [1:nr]+nr*(scy+scg);
        sc_h = A(:,1)*scy + A(:,2)*scx;
        
        rscx(tar_idx_x,tar_idx_y) = rx + sc_h(1);
        rscy(tar_idx_x,tar_idx_y) = ry + sc_h(2);
        ldos_sc(tar_idx_x,tar_idx_y) = reshape(ldos(tar_E,:),nr,nr);
        

    end
end

title('LDOS at red line')
surf(rscx, rscy, ldos_sc)
shading flat
view(2)
axis equal
colorbar
xlabel('$x$ (SC units)')
ylabel('$y$ (SC units)')