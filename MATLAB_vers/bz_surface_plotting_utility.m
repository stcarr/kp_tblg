%% Get full BZ data (for plutting a surface band structure over IBZ)
clear all;

theta_list = 1.0;%[0.9:.05:1.1];
knum = 3*7;

[sweep_vals, scaleaxis, sweep_kpts] = tblg_kp_calc_ext('theta_list',theta_list,'knum',knum,'full_bz',2,'vf_fac',1.2);
%%
bands = sweep_vals{1};
kpts = sweep_kpts{1};

nk = sqrt(size(kpts,1));
nb = size(bands,1);
kpts_sq = reshape(kpts,nk,nk,3);
bands_sq = reshape(bands',nk,nk,nb);

b1 = squeeze(kpts_sq( (nk+1)/2,      end,  1:2));
b2 = squeeze(kpts_sq(      end, (nk+1)/2,  1:2));

K_pt = (b1+b2)/3;

kpts_x = squeeze(kpts_sq(:,:,1));
kpts_y = squeeze(kpts_sq(:,:,2));

clf

K_cut_max = 1.25;

rot_mat = [cosd(60) sind(60); -sind(60) cosd(60)];
for r_idx = 0:5
    
    K1 = rot_mat^r_idx*K_pt*K_cut_max;
    K2 = rot_mat^(r_idx+1)*K_pt*K_cut_max;
    line_seg = K2-K1;

    for idx_x = 1:nk
        for idx_y = 1:nk
            
            R_h = squeeze(kpts_sq(idx_x,idx_y,1:2));
            off_seg = R_h-K2;
            if (sign(det([line_seg,off_seg])) ~= -1)
               bands_sq(idx_x,idx_y,:) = NaN; 
            end
            
            
            %plot3(K1(1),K1(2),.3,'xr')
            hold on
            
        end
        
    end
end


for tar_b = nb/2+[-6:5]

    bands_h = squeeze(bands_sq(:,:,tar_b));
    surf(kpts_x,kpts_y,bands_h);
    hold on
end
shading interp
axis equal
lighting gouraud
material dull
lightangle(30,40)
axis([-0.035 0.035 -0.035 0.035 -.25 .25])
%view(2)