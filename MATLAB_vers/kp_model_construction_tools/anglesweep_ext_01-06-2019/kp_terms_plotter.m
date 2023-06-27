% 1 to 41

for tar_theta = 1:length(thetas)
    inter_here_AB = squeeze(inter_kp(tar_theta,3,1,2,:));
    inter_here_AA = squeeze(inter_kp(tar_theta,3,1,1,:));

    c_here = (length(thetas)-tar_theta)/length(thetas)*[0 0 1] + tar_theta/length(thetas)*[1 0 0];
    
    subplot(2,1,1)
    plot(abs(inter_here_AB),'Color',c_here)
    title(num2str(thetas(tar_theta)*180/pi))
    hold on

    subplot(2,1,2)
    plot(abs(inter_here_AA),'Color',c_here)
    hold on

end


%%
clf
% 1 8 23 72 106

plot((thetas*180/pi),abs(squeeze(inter_kp(:,3,1,2,[1 8 23 72 106]))))