% README:
% The following three code sections generate band structures for TBLG
%   The first is a standard angle-sweep
%   The second adds symmetry-breaking terms
%   The third modifies the effective AA and AB interlayer couplings
% After this is a band-plotting utility, followed by code for DoS/IDoS

%% Get bandstructure data

theta_list = [0.9:.05:1.1];
[band_vals, bands_scaleaxis, band_kpts] = tblg_kp_calc_ext('theta_list',theta_list,'knum',10);


%% Play with layer (e.g. E-field) or sublattice (e.g. hBN) symmetry breaking terms
% see paper: Phys. Rev. Research 1, 033072 (2019).
%       url: https://journals.aps.org/prresearch/abstract/10.1103/PhysRevResearch.1.033072

theta_list = [1.1];

% this term breaks layer symmetry (lifts degeneracy at supercell K points)
displacement_strength =  0.075;% onsite energy (positive for layer 1 and negative for layer 2, in eV)

% these terms breaks sublattice symmetry (opens gap at CNP)
sublattice_strength_sym = 0.0; % sublattice symmetry breaking term, in eV. (identical on L1 and L2)
sublattice_strength_asym = 0.0; % sublattice symmetry breaking term, in eV. (opposite sign for L1 vs L2)

[band_vals, bands_scaleaxis, band_kpts] = tblg_kp_calc_ext('theta_list',theta_list,'knum',20,...
                                                           'displacement_strength',displacement_strength,...
                                                           'sublattice_strength_sym',sublattice_strength_sym,...
                                                           'sublattice_strength_asym',sublattice_strength_asym);

%% Play with interlayer coupling strengths
% (see: https://arxiv.org/abs/1910.07893)

theta_list = [1.0];
inter_AA_fac = 0.5; % screening can greatly reduce AA coupling
inter_AB_fac = 1.0; % but won't change AB much
[band_vals, bands_scaleaxis, band_kpts] = tblg_kp_calc_ext('theta_list',theta_list,'knum',20,'inter_aa_fac',inter_AA_fac,'inter_ab_fac',inter_AB_fac);

%% Plot bandstructure data
clf

% latex > tex?
f_size = 14;
set(groot, 'DefaultTextInterpreter', 'latex')
set(groot, 'DefaultLegendInterpreter', 'latex')
set(groot, 'DefaultAxesTickLabelInterpreter', 'latex')
set(groot, 'DefaultLegendInterpreter','latex')
set(groot, 'DefaultAxesFontSize',f_size)


ax_m = 0.15;

for idx = 1:length(theta_list)
   subplot(1,length(theta_list),idx)
   hold on
   box on
   
   bands = band_vals{idx};
   nb = size(bands,1);
   nk = size(bands,2)/2;
   E_f = bands(nb/2,1);
   
   plot(bands_scaleaxis{idx},bands(nb/2+[-30:30],1:nk)-E_f,'k')
   plot(bands_scaleaxis{idx},bands(nb/2+[-30:30],nk+1:end)-E_f,'k')

   axis([0 1 -ax_m ax_m])
   title(['$\theta =' num2str(theta_list(idx)) '^\circ$' ]);
   
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

