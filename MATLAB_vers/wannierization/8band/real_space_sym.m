function [h_filename] = real_space_sym(filename,tar_folder,tar_theta,rewan)
    % Real-space symmetrization for the 5/8 band models

    %% system paramters and load file (Stephen's)
    num_orbs = 8;
    %tar_theta = 1.1;
    %tar_folder = 'hmats_8band';
    %tar_filename = [tar_folder '/hmat_8band_' strrep(num2str(tar_theta,'%.2f'),'.','p') '.dat'];

    if (rewan == 0)

        h_filename = [tar_folder '/hmat_8band_' strrep(num2str(tar_theta,'%.2f'),'.','p') '_RS.dat'];


        fid = fopen(filename);

        temp = fgetl(fid);  % skip: Moire lattice vectors
        temp = str2num(fgetl(fid));  
        moire_L_x1 = [temp(1) temp(2) 0];
        temp = str2num(fgetl(fid));  
        moire_L_x2 = [temp(1) temp(2) 0];

        temp = fgetl(fid);  % skip: Moire reciprocal vectors
        temp = str2num(fgetl(fid));  
        moire_k_vec1 = [temp(1) temp(2) 0];
        temp = str2num(fgetl(fid));  
        moire_k_vec2 = [temp(1) temp(2) 0];

        temp = fgetl(fid);  % skip: Orbital locations
        for o = 1:num_orbs
            temp = str2num(fgetl(fid));  
            all_wan_xyz(o,:) = [temp(1) temp(2)];
        end

        temp = fgetl(fid);  % skip: Hamiltonian
        temp = fgetl(fid);  % skip: R_x R_y ... (etc)
        idx = 1;
        while (1)
            temp = fgetl(fid); % load the Hamiltonian data lines!
            if (temp == -1)
                break
            end
            H(idx,:) = str2num(temp);  
            idx = idx+1;
        end

        fclose(fid);
    elseif (rewan == 1)
        
        h_filename = [tar_folder '/hmat_8band_' strrep(num2str(tar_theta,'%.2f'),'.','p') '_rewan_RS.dat'];
        H = load(filename);
        
    end
    
    bonds_nosym = H(:,1:5);
    bonds_nosym(:,5) = bonds_nosym(:,5)+1j*H(:,6);

    
    %% Rephase WFs
    % transform c^d  -> xi c^d 
    WF_phase = [exp(1j*pi/4),exp(1j*3*pi/4),1,1j,1j,1,1,1];

    for brun = 1:size(bonds_nosym,1)
        xi_to = WF_phase(bonds_nosym(brun,3));
        xi_fr = WF_phase(bonds_nosym(brun,4));
        bonds_nosym(brun,5) = bonds_nosym(brun,5)*xi_to'*xi_fr;

    end

    %% Symmetry 
    phiphase = exp(1j*2*pi/3);


    if num_orbs == 5
        % to implement top vs bottom
        fprintf('5 band to be implemented. \n');
        return
    end


    sym_op(1).name = 'C3';
    sym_op(1).order = 3;
    sym_op(1).antiU = 0; % 0/1 for unitary/ anti-unitary
    sym_op(1).aMap = [-1,-1; 1,0]; % unit cell [l;m] maps to aMap*[l;m]
    sym_op(1).siteMap = [1,0,0,phiphase;...
                        2,0,0,phiphase';...
                        3,0,0,1;...
                        4,-1,0,1;...
                        5,1,0,1;...
                        7,0,0,1;...
                        8,0,0,1;...
                        6,0,0,1]; 
        % suppose c_n^d in the home cell is mapped to xi c_m^d at x a1 + y a_2,
        % then the n-th row of siteMap is [m,x,y,xi]



    sym_op(2).name = 'Mx';
    sym_op(2).order = 2;
    sym_op(2).antiU = 0; % 0/1 for unitary/ anti-unitary
    sym_op(2).aMap = [0,1; 1,0]; % unit cell [l;m] maps to aMap*[l;m]
    sym_op(2).siteMap = [2,0,0,-1;... % previously -1j
                        1,0,0,-1;... % previously 1j
                        3,0,0,1;...
                        4,0,0,-1;...
                        5,0,0,-1;...
                        6,-1,1,1;...
                        8,1,0,1;...
                        7,0,-1,1]; 

    sym_op(3).name = 'C2T';
    sym_op(3).order = 2;
    sym_op(3).antiU = 1; % 0/1 for unitary/ anti-unitary
    sym_op(3).aMap = [-1,0; 0,-1]; % unit cell [l;m] maps to aMap*[l;m]
    sym_op(3).siteMap = [2,0,0,-1;... % previously 1
                        1,0,0,-1;... % previously 1
                        3,0,0,1;...
                        5,0,0,1;... % previously -1
                        4,0,0,1;... % previously -1
                        6,-1,1,1;...
                        7,0,-1,1;...
                        8,1,0,1]; 







    %% Symmetrize
    bonds_sym = bonds_nosym;

    to_sym = [1,2,3];

    for s_now = to_sym
        bonds_sym = symmetrize(bonds_sym,sym_op(s_now));
    end

    H_sym = zeros(size(bonds_sym,1),6);
    H_sym(:,1:4) = bonds_sym(:,1:4);
    H_sym(:,5) = real(bonds_sym(:,5));
    H_sym(:,6) = imag(bonds_sym(:,5));

    non_zero_idx = find((abs(H_sym(:,5)) + abs(H_sym(:,6))) > 0);
    H_sym_nozeros = H_sym(non_zero_idx,:);

    fid = fopen(h_filename,'w');
    for idx = 1:size(H_sym_nozeros,1)
        R_x = H_sym_nozeros(idx,1);
        R_y = H_sym_nozeros(idx,2);
        m = H_sym_nozeros(idx,3);
        n = H_sym_nozeros(idx,4);
        t = H_sym_nozeros(idx,5) + 1j*H_sym_nozeros(idx,6);
        fprintf(fid," %s   %s   %s   %s   %s  %s \n",print_int(R_x), print_int(R_y), print_int(m), print_int(n), print_float(real(t)), print_float(imag(t)));
    end
    fclose(fid);
    %save(h_filename,'H_sym_nozeros','-ascii');

end

%% Symmetrize function
function [bonds_out] = symmetrize(bonds_in,sym_now)
    bond_tot = size(bonds_in,1);
    bonds_sym = zeros(bond_tot*sym_now.order,5);
    
    cnt = 0;
    for brun = 1:bond_tot
        bond = bonds_in(brun,:);
        cnt = cnt + 1;
        bonds_sym(cnt,:) = bond;
        
        
        for sym_run = 1:(sym_now.order-1)
            to = sym_now.siteMap(bond(3),:);
            fr = sym_now.siteMap(bond(4),:);

            r_vec = sym_now.aMap*[bond(1);bond(2)]+to(2:3).'-fr(2:3).';
            if sym_now.antiU==1
                t_now = bond(5)';
            else
                t_now = bond(5);
            end
            
            t_now = t_now*to(4)*fr(4)';
            
            bond = [r_vec(1),r_vec(2),to(1),fr(1),t_now];
            
            cnt = cnt + 1;
            bonds_sym(cnt,:) = bond;
        end
        
    end
    
    
    % clean up
    [bonds_data,~,remap] = unique(bonds_sym(:,1:4),'rows');
    bonds_out = zeros(size(bonds_data,1),5);
    bonds_out(:,1:4) = bonds_data;
    
    
    for brun = 1:length(remap)
        bonds_out(remap(brun),5) = bonds_out(remap(brun),5)+bonds_sym(brun,5);
    end
    
    bonds_out(:,5)=bonds_out(:,5)/sym_now.order;

end
    
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
    
    












