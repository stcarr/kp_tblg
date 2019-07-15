function [hmat_filename] = generate_hmat_datafile_5band(filename, tar_folder, tar_theta, cutoff_scale, side)

    if (side == 0)
        side_text = 'top';
    elseif (side == 1)
        side_text = 'bot';
    end

    num_orbs = 5;

    %cutoff_scale = 5;

    load(filename);

    % WAN 1,2,3: at 0
    % WAN 4 at (moire_L_x1+moire_L_x2)/3;
    % WAN 5 at -(moire_L_x1+moire_L_x2)/3;

    %moire_L_x1
    %moire_L_x2

    all_wan_xyz=zeros(5,3);
    all_wan_xyz(4,:)=(moire_L_x1+moire_L_x2)/3;
    all_wan_xyz(5,:)=-(moire_L_x1+moire_L_x2)/3;

    all_wan_xyz(6,:)=(moire_L_x1-moire_L_x2)/2;
    all_wan_xyz(7,:)=(moire_L_x2)/2;
    all_wan_xyz(8,:)=(-moire_L_x1)/2;


    %all_wan_xyz



    %{
    figure(1);
    clf

    scatter3(all_kpts(:,1),all_kpts(:,2),imag(all_new_hmat(1,4,:)));

    axis([-inf,inf,-inf,inf,-inf,inf]);
    %}

    %

    hop_cut=norm(moire_L_x1)*9.9;
    %hop_cut

    MMM=20;

    clear hex_R;
    clear idx_R;
    indh=1;
    for ind1=(-MMM):MMM
        for ind2=(-MMM):MMM
            vectmp=ind1*moire_L_x1+ind2*moire_L_x2;


            if norm(vectmp)<=hop_cut
                hex_R(indh,:)=vectmp(:);
                idx_R(indh,:)=[ind1, ind2];
                indh=indh+1;

            end


        end
    end
    num_hex=indh-1;


    %scatter(hex_R(:,1),hex_R(:,2))

    BZ_area=cross(moire_k_vec1,moire_k_vec2);
    BZ_area=abs(BZ_area(3));

    dBZ_area=BZ_area/numk/numk;


    all_R_mat=zeros(num_orbs,num_orbs,num_hex);



    for indr=1:num_hex
        R_hop=hex_R(indr,1:2);
        for indk=1:knum_tot
            know=all_kpts(indk,1:2);

            Heff_now=squeeze(all_new_hmat(1:num_orbs,1:num_orbs,indk));

            all_R_mat(:,:,indr)=all_R_mat(:,:,indr)+Heff_now*exp(1j*dot(know,R_hop))/(numk^2);

        end

        tmp=squeeze(all_R_mat(:,:,indr));
        maxRR(indr)=sum(abs(tmp(:)));

    end


    % scatter3(hex_R(:,1),hex_R(:,2),abs(all_R_mat(5,4,:)));

    % symmetric cutoff
    TBH_Rcut=norm(moire_L_x1)*cutoff_scale;

    Rcut_mask=zeros(num_orbs,num_orbs,num_hex);

    for indr=1:num_hex
        bigR=hex_R(indr,1:2);
        for indp1=1:num_orbs
            pos1=all_wan_xyz(indp1,1:2);
            for indp2=1:num_orbs
                pos2=all_wan_xyz(indp2,1:2);


                r_connect=bigR+pos2-pos1;
                if norm(r_connect)<TBH_Rcut
                    Rcut_mask(indp2,indp1,indr)=1;
                end

            end

        end
    end

    all_R_mat_cut=zeros(num_orbs,num_orbs,num_hex);
    for indr=1:num_hex
        tmph=squeeze(all_R_mat(:,:,indr));
        tmpm=squeeze(Rcut_mask(:,:,indr));

        all_R_mat_cut(:,:,indr)=tmph.*tmpm;


    end

    %save_folder = 'hmats_5band';
    %mkdir(save_folder);
    save_folder = tar_folder;
    hmat_filename = [save_folder '/hmat_5band_' side_text '_' strrep(num2str(tar_theta,'%.2f'),'.','p') '.dat'];

    fileID = fopen(hmat_filename,'w');

    fprintf(fileID,"Moire lattice vectors \n");
    fprintf(fileID," %s   %s \n", print_float(moire_L_x1(1)), print_float(moire_L_x1(2)));
    fprintf(fileID," %s   %s \n", print_float(moire_L_x2(1)), print_float(moire_L_x2(2)));

    fprintf(fileID,"Moire reciprocal vectors \n");
    fprintf(fileID," %s   %s \n", print_float(moire_k_vec1(1)), print_float(moire_k_vec1(2)));
    fprintf(fileID," %s   %s \n", print_float(moire_k_vec2(1)), print_float(moire_k_vec2(2)));

    fprintf(fileID,"Orbital locations \n");
    for idx = 1:num_orbs
        fprintf(fileID," %s   %s \n", print_float(all_wan_xyz(idx,1)), print_float(all_wan_xyz(idx,2)));
    end

    fprintf(fileID," Hamiltonian \n");
    fprintf(fileID," %s  %s   %s    %s         %s           %s \n","R_x", "R_y", "m", "n", "t_real", "t_cpx");
    for idx = 1:size(idx_R,1)
        R_x = idx_R(idx,1);
        R_y = idx_R(idx,2);
        for m = 1:num_orbs
            for n = 1:num_orbs
                t = all_R_mat_cut(m,n,idx);
                if (t ~= 0)
                    fprintf(fileID," %s   %s   %s   %s   %s  %s \n",print_int(R_x), print_int(R_y), print_int(m), print_int(n), print_float(real(t)), print_float(imag(t)));
                end
            end
        end 
    end
    fclose(fileID);
end


% functions
function [str] = print_int(x)
    if (x < 0)
        str = string(x);
    else
        str = " "+string(x);
    end

end

function [str] = print_float(x)
    if (x < 0)
        str = num2str(x,'%.10f');
    else
        str = " "+num2str(x,'%.10f');
    end
    
    if abs(x) < 10
       str = " "+str; 
    end
    if abs(x) < 100
       str = " "+str; 
    end

end
