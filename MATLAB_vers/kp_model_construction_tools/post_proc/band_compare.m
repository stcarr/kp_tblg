
run(global_var_file);

sc_m=all_sc_m(1);
sc_n=sc_m-1;

koshino_bands_filename = ['TwBLG_tbh-bands_',num2str(sc_m),'_',num2str(sc_n),'_koshino-inter'];
dft_bands_filename = ['TwBLG_tbh-bands_',num2str(sc_m),'_',num2str(sc_n),'_dft-inter'];

 %{
fprintf(['loading file: ' koshino_bands_filename '.mat \n']);
load([bands_data_dir '/' koshino_bands_filename]);

clf
hold on

nb = size(allbands,1);
plot(scale_axis,allbands' - allbands(nb/2,1),'k')
 %}

% %{
fprintf(['loading file: ' dft_bands_filename '.mat \n']);
load([bands_data_dir '/' dft_bands_filename]);

nb = size(allbands,1);
plot(scale_axis,allbands' - allbands(nb/2,1),'r')
% %}

%axis([0 1 -2 2])
hold off