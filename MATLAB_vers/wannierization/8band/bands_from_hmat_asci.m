function [bands, scale_axis] =  bands_from_hmat_asci(filename, num_orbs)

    %ax_m = 0.3;

    %num_orbs = 8;

    kk1 = [0,0,0];
    kk2 = [-1/2,0,0];
    kk3 = [-1/3,1/3,0];

    H_in = load(filename);
    bonds = zeros(size(H_in,1),5);
    bonds(:,1:5) = H_in(:,1:5);
    bonds(:,5) = bonds(:,5) + 1j*H_in(:,6);

    [Har_list,~,Har_ind] = unique(bonds(:,1:2),'rows');

    H_Har = zeros(num_orbs,num_orbs,max(Har_ind));
    for brun = 1:length(Har_ind)
        bond_now = bonds(brun,3:5);

        H_Har(bond_now(1),bond_now(2),Har_ind(brun)) = ...
            H_Har(bond_now(1),bond_now(2),Har_ind(brun)) + bond_now(3);
        % provision for repeated bonds although it shouldn't happen
    end

    Nk = 31;

    Hermitianize = @(h) (h+h')/2;
    Hk = @(g) Hermitianize(sum(...
                bsxfun(@times,H_Har,reshape(exp(-1j*2*pi*Har_list*g.'),1,1,[])),...
                    3)); 
    kscan_list=[kk3;kk1;kk2;kk3];
    %kscan_list=[kk1;kk2;kk3;kk1];

    [ all_kpts, scale_axis] = generate_k_line( Nk, kscan_list);

    all_kpts = all_kpts(:,1:2);

    bands = zeros(num_orbs,size(all_kpts,1));
    weights = zeros(num_orbs,num_orbs,size(all_kpts,1));
    for krun = 1:size(all_kpts,1)

        Hk_now = Hk(all_kpts(krun,:));
        [vecs, vals] = eig(Hk_now);
        [bands(:,krun), ord] = sort(diag(vals));
        %weights(:,:,krun) = abs(vecs(:,ord)).^2;
    end

    %E_f = bands(4,1);
    %G_pt = scale_axis(32);
    %M_pt = scale_axis(63);

    %plot(scale_axis_proj,bands_proj'-E_f,'Color','k')

end
        