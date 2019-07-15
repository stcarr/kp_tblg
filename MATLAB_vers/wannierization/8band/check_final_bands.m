
num_orbs = 8;
tar_theta = 1;
tar_theta_str = strrep(num2str(tar_theta,'%.2f'),'.','p');
clf
tar_folder = 'run_folder_8band_final_07-10-2019';
filename_orig = [tar_folder '/hmat_8band_' tar_theta_str '_rewan.dat'];
filename_RS = [tar_folder '/hmat_8band_' tar_theta_str '_rewan_RS.dat'];


[bands_orig, scale_axis] =  bands_from_hmat_asci(filename_orig, num_orbs);
[bands_RS, scale_axis] =  bands_from_hmat_asci(filename_RS, num_orbs);

subplot(1,2,1)
plot(scale_axis,bands_orig,'k');
subplot(1,2,2)
plot(scale_axis,bands_RS,'r');