clear all;

load dft_full_relax_data_01-06-2019;

thetas_deg = thetas*180/pi;

for theta_idx = 1:length(thetas_deg)

        All_Eff_inter_ext = squeeze(inter_kp(theta_idx,:,:,:,:));
        All_Eff_intra_bot_ext = squeeze(intra_bot_kp(theta_idx,:,:,:,:));
        All_Eff_intra_top_ext = squeeze(intra_top_kp(theta_idx,:,:,:,:));
        
        % clean up NAN errors during fitting (small terms -> singular)
        All_Eff_inter=squeeze(All_Eff_inter_ext(3,:,:,:));
        All_Eff_intra_bot=squeeze(All_Eff_intra_bot_ext(3,:,:,:))*1;
        All_Eff_intra_top=squeeze(All_Eff_intra_top_ext(3,:,:,:))*1;
        All_Eff_inter_kplus=squeeze(All_Eff_inter_ext(1,:,:,:));
        All_Eff_inter_kminus=squeeze(All_Eff_inter_ext(2,:,:,:));
        
        All_Eff_inter(isnan(All_Eff_inter)) = 0;
        All_Eff_intra_bot(isnan(All_Eff_intra_bot)) = 0;
        All_Eff_intra_top(isnan(All_Eff_intra_top)) = 0;
        All_Eff_inter_kplus(isnan(All_Eff_inter_kplus)) = 0;
        All_Eff_inter_kminus(isnan(All_Eff_inter_kminus)) = 0;
        
        All_Eff_intra_shell_indices = intra_shells;
        All_Eff_inter_shell_indices = inter_shells;
        
        % Symm enforcer, to make sure the numerically calculated terms
        % do not break any assumed symmetery of our k-p Hamiltonian
        [All_Eff_inter,All_Eff_inter_kplus,All_Eff_inter_kminus] = TwBLG_KP_sym_enforce_inter(All_Eff_inter,All_Eff_inter_kplus,All_Eff_inter_kminus);
         
        inter{theta_idx} = All_Eff_inter;
        inter_kplus{theta_idx} = All_Eff_inter_kplus;
        inter_kminus{theta_idx} = All_Eff_inter_kminus;
        intra_bot{theta_idx} = All_Eff_intra_bot;
        intra_top{theta_idx} = All_Eff_intra_top;
        
end

ntheta = length(thetas_deg);
n_inter_couplings = size(inter{1},3);
n_intra_couplings = size(intra_bot{1},3);


tar_filename = 'full_relax_kp_01-06-2019.dat';



fileID = fopen(tar_filename,'w');
fprintf(fileID,'%d %d %d %d %d \n',ntheta,n_inter_couplings,n_intra_couplings,length(inter_shells),length(intra_shells)-1);

fprintf(fileID,"thetas (deg): \n");
for t_idx = 1:ntheta
    fprintf(fileID,"%.10f ",thetas_deg(t_idx));
end
fprintf(fileID,"\n");

fprintf(fileID,"inter terms: \n");
for t_idx = 1:ntheta
    inter_here = inter{t_idx};
    for p_idx = 1:n_inter_couplings
        for o1 = 1:2
            for o2 = 1:2
                if imag(inter_here(o1,o2,p_idx)) < 0
                    fprintf(fileID,"%.10f%.10fi ",real(inter_here(o1,o2,p_idx)),imag(inter_here(o1,o2,p_idx)));
                else
                    fprintf(fileID,"%.10f+%.10fi ",real(inter_here(o1,o2,p_idx)),imag(inter_here(o1,o2,p_idx)));
                end
            end
        end
    end
    fprintf(fileID,"\n")

end

fprintf(fileID,"inter kplus terms: \n");
for t_idx = 1:ntheta
    inter_here = inter_kplus{t_idx};
    for p_idx = 1:n_inter_couplings
        for o1 = 1:2
            for o2 = 1:2
                if imag(inter_here(o1,o2,p_idx)) < 0
                    fprintf(fileID,"%.10f%.10fi ",real(inter_here(o1,o2,p_idx)),imag(inter_here(o1,o2,p_idx)));
                else
                    fprintf(fileID,"%.10f+%.10fi ",real(inter_here(o1,o2,p_idx)),imag(inter_here(o1,o2,p_idx)));
                end
            end
        end
    end
    fprintf(fileID,"\n")

end

fprintf(fileID,"inter kminus terms: \n");
for t_idx = 1:ntheta
    inter_here = inter_kminus{t_idx};
    for p_idx = 1:n_inter_couplings
        for o1 = 1:2
            for o2 = 1:2
                if imag(inter_here(o1,o2,p_idx)) < 0
                    fprintf(fileID,"%.10f%.10fi ",real(inter_here(o1,o2,p_idx)),imag(inter_here(o1,o2,p_idx)));
                else
                    fprintf(fileID,"%.10f+%.10fi ",real(inter_here(o1,o2,p_idx)),imag(inter_here(o1,o2,p_idx)));
                end
            end
        end    
    end
    fprintf(fileID,"\n")

end

fprintf(fileID,"intra bot terms: \n");
for t_idx = 1:ntheta
    intra_here = intra_bot{t_idx};
    for p_idx = 1:n_intra_couplings
        for o1 = 1:2
            for o2 = 1:2
                if imag(intra_here(o1,o2,p_idx)) < 0
                    fprintf(fileID,"%.10f%.10fi ",real(intra_here(o1,o2,p_idx)),imag(intra_here(o1,o2,p_idx)));
                else
                    fprintf(fileID,"%.10f+%.10fi ",real(intra_here(o1,o2,p_idx)),imag(intra_here(o1,o2,p_idx)));
                end
            end
        end
    end
    fprintf(fileID,"\n")

end

fprintf(fileID,"intra top terms: \n");
for t_idx = 1:ntheta
    intra_here = intra_top{t_idx};
    for p_idx = 1:n_intra_couplings
        for o1 = 1:2
            for o2 = 1:2
                if imag(intra_here(o1,o2,p_idx)) < 0
                    fprintf(fileID,"%.10f%.10fi ",real(intra_here(o1,o2,p_idx)),imag(intra_here(o1,o2,p_idx)));
                else
                    fprintf(fileID,"%.10f+%.10fi ",real(intra_here(o1,o2,p_idx)),imag(intra_here(o1,o2,p_idx)));
                end
            end
        end
    end
    fprintf(fileID,"\n")

end

fprintf(fileID,"inter shells: \n");
for s_idx = 1:length(inter_shells)
    fprintf(fileID,"%d ",size(inter_shells{s_idx},1));
end
fprintf(fileID,"\n");

for s_idx = 1:length(inter_shells)
    shell_here = inter_shells{s_idx};
    for idx = 1:size(shell_here,1)
        fprintf(fileID,"%d %d %d \n",shell_here(idx,:));
    end
end

fprintf(fileID,"intra shells: \n");
for s_idx = 1:length(intra_shells)-1
    fprintf(fileID,"%d ",size(intra_shells{s_idx+1},1));
end
fprintf(fileID,"\n");

for s_idx = 1:length(intra_shells)
    shell_here = intra_shells{s_idx};
    for idx = 1:size(shell_here,1)
        fprintf(fileID,"%d %d %d \n",shell_here(idx,:));
    end
end


fclose(fileID);
