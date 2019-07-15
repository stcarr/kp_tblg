
% MODEL_TYPE: Which Reduced model to use
% 1: 8-band
% 2: 5-band (top)
% 3: 5-band (bot)
model_type = 1; 

% TAR_THETA: twist angle in degrees, only supports
% [0.60 : 0.05 : 1.30] for 8-band model
% [0.85 : 0.05 : 1.30] for 5-band models
% TODO: Interpolate hamiltonians between the sampled points?
tar_theta = 1.05;

%
% SETUP AND CALCULATE BAND STRUCTURE
%

if (model_type == 1)
    num_orbs = 8;
    model_str = '8band';
elseif (model_type == 2)
    num_orbs = 5;
    model_str = '5band_top';
elseif (model_type == 3)
    num_orbs = 5;
    model_str = '5band_bot';
end

tar_folder = 'hmats_2019-07-11'; % Hamiltonians generated on July 11th, 2019

tar_theta_str = strrep(num2str(tar_theta,'%.2f'),'.','p');
clf

filename = [tar_folder '/hmat_' model_str '_' tar_theta_str '_rewan_RS.dat'];


[bands, scale_axis] =  bands_from_hmat_asci(filename, num_orbs);

Nk = length(scale_axis)/3;
Gamma_pt = scale_axis(Nk+1);
M_pt = scale_axis(2*Nk+1);
ax_m = 0.3;

plot(scale_axis,bands,'k');
axis([0 1 -ax_m ax_m]);
set(gca,'XTick',[0 Gamma_pt M_pt 1]);
xticklabels({'K','G','M','K'})
ylabel('Energy (eV)')
