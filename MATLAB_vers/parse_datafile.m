function [thetas, inter_kp, intra_bot_kp, intra_top_kp, inter_shells, intra_shells] = parse_datafile(filename)
    
    %load dft_full_relax_data_02-04-2019.mat;
    %filename = '../data/full_relax_kp_01-06-2019.dat';
    
    fid = fopen(filename);

    sizes = str2num(fgetl(fid));  % get array sizes
    n_theta = sizes(1);
    n_inter = sizes(2);
    n_intra = sizes(3);
    num_inter_shells = sizes(4);
    num_intra_shells = sizes(5);

    fgetl(fid); % skip text
    thetas = pi/180*str2num(fgetl(fid));  % get theta list

    inter_kp = zeros(n_theta,3,2,2,n_inter);
    intra_bot_kp = zeros(n_theta,3,2,2,n_intra);
    intra_top_kp = zeros(n_theta,3,2,2,n_intra);

    fgetl(fid); % skip text
    for t_idx = 1:n_theta
        data_idx = 1;
        data_list = str2num(fgetl(fid));  % get theta list    
        for p_idx = 1:n_inter
            for o1 = 1:2
                for o2 = 1:2
                    inter_kp(t_idx,3,o1,o2,p_idx) = data_list(data_idx);
                    data_idx = data_idx+1;
                end
            end
        end
    end

    fgetl(fid); % skip text "inter kplus terms"
    for t_idx = 1:n_theta
        data_idx = 1;
        data_list = str2num(fgetl(fid));  % get theta list    
        for p_idx = 1:n_inter
            for o1 = 1:2
                for o2 = 1:2
                    inter_kp(t_idx,1,o1,o2,p_idx) = data_list(data_idx);
                    data_idx = data_idx+1;
                end
            end
        end
    end

    fgetl(fid); % skip text "inter kminus terms"
    for t_idx = 1:n_theta
        data_idx = 1;
        data_list = str2num(fgetl(fid));  % get theta list    
        for p_idx = 1:n_inter
            for o1 = 1:2
                for o2 = 1:2
                    inter_kp(t_idx,2,o1,o2,p_idx) = data_list(data_idx);
                    data_idx = data_idx+1;
                end
            end
        end
    end

    fgetl(fid); % skip text "intra bot terms"
    for t_idx = 1:n_theta
        data_idx = 1;
        data_list = str2num(fgetl(fid));  % get theta list    
        for p_idx = 1:n_intra
            for o1 = 1:2
                for o2 = 1:2
                    intra_bot_kp(t_idx,3,o1,o2,p_idx) = data_list(data_idx);
                    data_idx = data_idx+1;
                end
            end
        end
    end

    fgetl(fid); % skip text "intra top terms"
    for t_idx = 1:n_theta
        data_idx = 1;
        data_list = str2num(fgetl(fid));  % get theta list    
        for p_idx = 1:n_intra
            for o1 = 1:2
                for o2 = 1:2
                    intra_top_kp(t_idx,3,o1,o2,p_idx) = data_list(data_idx);
                    data_idx = data_idx+1;
                end
            end
        end
    end


    fgetl(fid);% skip text "inter shells:"
    shell_size_list = str2num(fgetl(fid));  % get theta list   
    for s_idx = 1:num_inter_shells
        temp_shells = zeros(shell_size_list(s_idx),3);
        for idx = 1:shell_size_list(s_idx)
            temp_shells(idx,:) = str2num(fgetl(fid));
        end
        inter_shells{s_idx} = temp_shells;
    end

    fgetl(fid);% skip text "intra shells:"
    shell_size_list = str2num(fgetl(fid));  % get theta list   
    for s_idx = 1:num_intra_shells
        temp_shells = zeros(shell_size_list(s_idx),3);
        for idx = 1:shell_size_list(s_idx)
            temp_shells(idx,:) = str2num(fgetl(fid));
        end
        intra_shells{s_idx+1} = temp_shells;    
    end

    fclose(fid);
    end