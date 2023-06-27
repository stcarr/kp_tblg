%
main_dir = pwd;

% supercell choice

all_sc_m= 30; %[8 9 10 11];
all_rot_angle=zeros(length(all_sc_m),1);

% relaxation factor:
fft_relax_onoff_factor= 0.0;

% bandstructure sampling (total k points = knum*3)
knum = 20;


% two choices for interlayer coupling here:

inter_type = 1;

if (inter_type == 0)
    interlayer_model = 'koshino';
elseif (inter_type == 1)
    interlayer_model = 'dft';
end

kp_fit_range = 1;

% fit_range = 1 -> only nearest q point
% fit_range = 2 -> includes NN q point
% fit_range = 3 -> includes NNN q point
