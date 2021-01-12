function [dos_sweep, idos_sweep, E_list, half_filling_hole_E] = interp_kp_dos(theta_list, sweep_vals, sweep_kpts)

    tar_theta_list = theta_list;

    % number of extra bands to include
    b_size = 4;


    for t_idx = 1:length(tar_theta_list)
        fprintf("Starting theta %d/%d: ",t_idx,length(tar_theta_list));
        tic
        
        tri_idx = 1;

        all_kpts1 = sweep_kpts{t_idx};
        allbands1 = sweep_vals{t_idx};

        kpoints = all_kpts1(:,[2 1]);
        eig_vals = allbands1';

        nk = sqrt(size(eig_vals,1));
        nb = size(eig_vals,2);

        n_tri = nk*nk*(2*b_size + 2);

        triangles1 = zeros(n_tri,3);
        triangles2 = zeros(n_tri,3);
        k_tris1 = zeros(n_tri,3,2);
        k_tris2 = zeros(n_tri,3,2);



        for tar_b = (nb/2)-b_size:(nb/2+1)+b_size
        %for tar_b = 1;
        %fprintf("on band %d / %d \n",tar_b - ((nb/2)-b_size) + 1, 2*b_size + 2);

            for i = 1:nk
                for j = 1:nk
                    tar_band(i,j) = eig_vals((i-1)*nk + j,tar_b);
                    k_mesh(i,j,:) = kpoints((i-1)*nk + j,:);
                    %k_mesh(i,j,1) = 0.1*(j-nk/2);
                    %k_mesh(i,j,2) = 0.1*(i-nk/2);

                    %tar_band(i,j) = norm( squeeze(k_mesh(i,j,:)) ).^2;
                end
            end
            %nk = sqrt(size(k_mesh,1)*size(k_mesh,2));

            % make triangles!
            for i = 1:nk
                for j = 1:nk

                    ip = i+1;
                    im = i-1;
                    jp = j+1;
                    jm = j-1;

                    if (ip > nk)
                       ip = 1; 
                    end
                    if (im < 1)
                       im = nk; 
                    end
                    if (jp > nk)

                       jp = 1; 
                    end
                    if (jm < 1)
                       jm = nk; 
                    end

                    triangles1(tri_idx,1) = tar_band(i,j);
                    triangles1(tri_idx,2) = tar_band(ip,j);
                    triangles1(tri_idx,3) = tar_band(i,jp);
                    k_tris1(tri_idx,1,:) = k_mesh(1,1,:);
                    k_tris1(tri_idx,2,:) = k_mesh(2,1,:);
                    k_tris1(tri_idx,3,:) = k_mesh(1,2,:);
                    k_0s1(tri_idx,:) = k_mesh(i,j,:);
                    %k_tris1(tri_idx,1,:) = k_tris1(tri_idx,1,:) + k_mesh(i,j,:);
                    %k_tris1(tri_idx,2,:) = k_tris1(tri_idx,2,:) + k_mesh(i,j,:);
                    %k_tris1(tri_idx,3,:) = k_tris1(tri_idx,3,:) + k_mesh(i,j,:);

                    triangles2(tri_idx,1) = tar_band(i,j);
                    triangles2(tri_idx,2) = tar_band(im,j);
                    triangles2(tri_idx,3) = tar_band(i,jm);
                    k_tris2(tri_idx,1,:) = k_mesh(end,end,:);
                    k_tris2(tri_idx,2,:) = k_mesh(end-1,end,:);
                    k_tris2(tri_idx,3,:) = k_mesh(end,end-1,:);
                    %k_tris2(tri_idx,1,:) = k_tris2(tri_idx,1,:) + k_mesh(i,j,:);
                    %k_tris2(tri_idx,2,:) = k_tris2(tri_idx,2,:) + k_mesh(i,j,:);
                    %k_tris2(tri_idx,3,:) = k_tris2(tri_idx,3,:) + k_mesh(i,j,:);


                    tri_idx = tri_idx+1;            

                end
            end
        end
       
        max_E = 0.3;
        dE = max_E/6000;
        
        E_list = [-max_E:dE:max_E];
        dos = zeros(length(E_list),1);

        for E_idx = 1:length(E_list)
            E = E_list(E_idx);
            for t = 1:length(triangles1)
               if E > min(triangles1(t,:)) && E < max(triangles1(t,:))  
                %dos(E_idx) = dos(E_idx)+1;
                % replace this with proper gradient computation function!
                % need slope |b| and cross sectional area (length) f!
                
                % our triangular mesh element is spanned by two vectors, v,w
                
                % v is dk1
                v(1) = k_tris1(t,2,1) - k_tris1(t,1,1);
                v(2) = k_tris1(t,2,2) - k_tris1(t,1,2);
                % w is dk2
                w(1) = k_tris1(t,3,1) - k_tris1(t,1,1);
                w(2) = k_tris1(t,3,2) - k_tris1(t,1,2);
                % NOTE: here w(2) == 0!

                % E0 is the difference between the sampled Energy and the
                % vertex of the triangular mesh element
                E0 = E - triangles1(t,1);
                % Ev,Ew are the differences between the vertex and the v,w
                % points of the triangular mesh element
                Ev = triangles1(t,2) - triangles1(t,1);
                Ew = triangles1(t,3) - triangles1(t,1);

                if (w(2) ~= 0)
                   fprintf('WARNING: dk2 is not purely along x dir! \n');
                   pause(10) 
                end

                % b is the gradient, we use w ~ x-hat to quickly compute it
                b(1) = Ew/w(1);
                b(2) = (Ev - b(1)*v(1)) / v(2);

                % For w(2) not equal to zero:
                %{ 
                b(1) = Ew*v(1) - Ev*w(1);
                b(2) = Ev*w(2) - Ew*v(2);
                b = b/( w(2)*v(1) - w(1)*v(2));
                %}

                % The "t"s tell us where the "E0 energy cross section" 
                % crosses the boundaries of the triangular mesh element.
                % I.e: one of these has to fall outside of the finite range
                % [0,1], and then we can calculate the cross-sectional area
                % by using the other two (which will both be within [0,1]!)
                t1 = E0/dot(b,v);
                t2 = E0/dot(b,w);
                t3 = (E0 - dot(b,v)) / dot(b,w-v);

                f = 0;

                if (t3 < 0 || t3 > 1)
                    p1 = t1*v;
                    p2 = t2*w;
                    f = sqrt( sum((p1 - p2).^2) );
                end
                if (t1 < 0 || t1 > 1)
                    p2 = t2*w;
                    p3 = v + t3*(w - v);
                    f = sqrt( sum((p2 - p3).^2) );
                end
                if (t2 < 0 || t2 > 1)
                    p1 = t1*v;
                    p3 = v + t3*(w - v);
                    f = sqrt( sum((p1 - p3).^2) );
                end

                if(norm(b) == 0)
                    fprintf('WARNING: b = 0! \n');
                    pause(10) 
                end
                if (f/norm(b) > 10)
                    fprintf('WARNING!! \n');
                    pause(10)                
                end

                dos(E_idx) = dos(E_idx) + f/norm(b);
                %dos_tots1(E_idx,t) = f/norm(b);

               end
               if E > min(triangles2(t,:)) && E < max(triangles2(t,:))  
                % replace this with proper gradient computation function!
                % need slope |b| and cross sectional area (length) 
                % v is dk1
                v(1) = k_tris2(t,2,1) - k_tris2(t,1,1);
                v(2) = k_tris2(t,2,2) - k_tris2(t,1,2);
                % w is dk2

                w(1) = k_tris2(t,3,1) - k_tris2(t,1,1);
                w(2) = k_tris2(t,3,2) - k_tris2(t,1,2);

                E0 = E - triangles2(t,1);
                Ev = triangles2(t,2) - triangles2(t,1);
                Ew = triangles2(t,3) - triangles2(t,1);

                if (w(2) ~= 0)
                   fprintf('WARNING: dk2 is not purely along x dir! \n');
                   pause(10) 
                end

                b(1) = Ew/w(1);
                b(2) = (Ev - b(1)*v(1)) / v(2);

                % For w(2) not equal to zero:
                %{ 
                b(1) = Ew*v(1) - Ev*w(1);
                b(2) = Ev*w(2) - Ew*v(2);
                b = b/( w(2)*v(1) - w(1)*v(2));
                %}

                t1 = E0/dot(b,v);
                t2 = E0/dot(b,w);
                t3 = (E0 - dot(b,v)) / dot(b,w-v);

                f = 0;

                if (t3 < 0 || t3 > 1)
                    p1 = t1*v;
                    p2 = t2*w;
                    f = sqrt( sum((p1 - p2).^2) );
                end
                if (t1 < 0 || t1 > 1)
                    p2 = t2*w;
                    p3 = v + t3*(w - v);
                    f = sqrt( sum((p2 - p3).^2) );
                end
                if (t2 < 0 || t2 > 1)
                    p1 = t1*v;
                    p3 = v + t3*(w - v);
                    f = sqrt( sum((p1 - p3).^2) );
                end

                if(norm(b) == 0)
                    fprintf('WARNING: b = 0! \n');
                    pause(10) 
                end
                if (f/norm(b) > 10)
                    fprintf('WARNING!! \n');
                    pause(10)                
                end

                dos(E_idx) = dos(E_idx) + f/norm(b);
                %dos_tots2(E_idx,t) = f/norm(b);

               end
            end
            %E_idx/length(E_list)
            
        end

        %dos_sweep{t_idx} = dos;

        idos = zeros(size(dos));
        for x = 1:length(E_list)
            if (x > 1)
                idos(x) = trapz(E_list(1:x),dos(1:x));
            end
        end

        %tot_bands = 4*2*(b_size+1);% 2 for valley, 2 for spin
        
        alpha = 2.47;
        sc_alpha = alpha/(2*sind(tar_theta_list(t_idx)/2));
        sc_area = sc_alpha^2*sind(60)*1e-2; %area in nm^2
        %n0 = 1/sc_area;
        
        %idos_rescale = tot_bands/idos(end);
        %dos_rescale = idos_rescale*n0
        
        % 100 for A^2 -> nm^2, 4 for valley/spin, 4pi^2 for 
        dos_rescale = 100*4/(2*pi)^2;
        idos_rescale = dos_rescale*sc_area;

        idos(:) = idos_rescale*(idos(:) - 0.5*idos(end));
        [val, idx] = min(abs(idos - (-2)));
        
        
        idos_sweep{t_idx} = idos;
        dos_sweep{t_idx} = dos_rescale*dos;
        half_filling_hole_E(t_idx) = E_list(idx);
    
        toc
    end
    
    % old plotting utilities
    %{
    clf
    hold on
    plot(idos,dos_rescale*dos,'k','LineWidth',2);
    %text(-14,dos_max*.9,'unrelaxed 0.48^\circ')

    dos_max = 60;
    m = 6;
    axis([-m m 0 dos_max])
    set(gca,'XTick',[-m+rem(m,4):4:m])
    plot([idos(idx) idos(idx)],[0 dos_max],'--k')
    %xticklabels({})
    set(gca,'YTick',[])
    ylabel('DoS')
    xlabel('n/n0')
    %}

    %{
    tar_b = (nb/2);
    for i = 1:nk
        for j = 1:nk
            tar_band(i,j) = eig_vals((i-1)*nk + j,tar_b);
            k_mesh(i,j,:) = kpoints((i-1)*nk + j,:);
            %k_mesh(i,j,1) = 0.1*(j-nk/2);
            %k_mesh(i,j,2) = 0.1*(i-nk/2);

            %tar_band(i,j) = norm( squeeze(k_mesh(i,j,:)) ).^2;
        end
    end
    surf(k_mesh(:,:,1),k_mesh(:,:,2),tar_band,'EdgeColor','none','FaceColor','r')
    axis([-inf inf -inf inf -inf half_filling_hole_E(end)])
    view(2)
    %}
    
end