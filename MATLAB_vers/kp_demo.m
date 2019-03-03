%% Get bandstructure data

theta_list = [0.9:.05:1.1];
[band_vals, bands_scaleaxis, band_kpts] = tblg_kp_calc_ext('theta_list',theta_list,'knum',10);
%% Plot bandstructure data
clf

ax_m = 0.1;

for idx = 1:length(theta_list)
   subplot(1,length(theta_list),idx)
   hold on
   
   bands = band_vals{idx};
   nb = size(bands,1);
   nk = size(bands,2)/2;
   E_f = bands(nb/2,1);
   
   plot(bands_scaleaxis{idx},bands(nb/2+[-30:30],1:nk)-E_f,'k')
   plot(bands_scaleaxis{idx},bands(nb/2+[-30:30],nk+1:end)-E_f,'k')

   axis([0 1 -ax_m ax_m])
   title(['$' num2str(theta_list(idx)) '^\circ$' ]);
   
   if(idx == 1)
       ylabel('Energy (eV)')
   else
      yticklabels({}) 
   end
   
end

%% Get full BZ data
clear all;

theta_list = [0.9:.05:1.1];

[sweep_vals, scaleaxis, sweep_kpts] = tblg_kp_calc_ext('theta_list',theta_list,'knum',10,'full_bz',1);

[dos_sweep, idos_sweep, E_list, half_filling_hole_E] = interp_kp_dos(theta_list, sweep_vals, sweep_kpts);

%% Plot DoS/IDoS data
% not great resolution, need higher knum for good results on the flat-bands!

clf

ax_m = .14;

for idx = 1:length(theta_list)
   subplot(length(theta_list),1,length(theta_list)-idx+1)
   hold on
   
   bands = band_vals{idx};
   nb = size(bands,1);
   nk = size(bands,2)/2;
   E_f = bands(nb/2,1);
   
   dos_here = dos_sweep{idx};
   idos_here = idos_sweep{idx};
   plot(E_list-E_f,dos_here)
   
   axis([-ax_m ax_m 0 10])
   text(-ax_m*0.9,50,['$' num2str(theta_list(idx)) '^\circ$'],'FontSize',16)
 
   if (idx == 1)
       xlabel('Energy (eV)')
   else
      xticklabels({}) 
   end
   
   ylabel('DoS')

   % DoS vs IDoS (carrier filling)
   % need to manually set Fermi-level
   % as the k-dot-p model can not predict it
   %plot(idos_here,dos_here)

   
end

