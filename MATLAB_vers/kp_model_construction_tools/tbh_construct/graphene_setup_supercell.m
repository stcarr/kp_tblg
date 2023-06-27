% setup basic environment for supercell

% input: sc_m, sc_n

% expect sc_m >=2    sc_n>=1


lattice_a=1.42*sqrt(3);

layer_d=[0,0,3.35];

% lattice_a*
a1=[sqrt(3)/2,-1/2]';
a2=[sqrt(3)/2,+1/2]';


sc_int_bottom_a1=[sc_n,sc_m];
sc_int_bottom_a2=[-sc_m,sc_n+sc_m];

sc_int_top_a1=[sc_m,sc_n];
sc_int_top_a2=[-sc_n,sc_n+sc_m];


num_pc=(sc_m^2+sc_m*sc_n+sc_n^2);
total_num=4*num_pc;

rot_theta=acos((sc_m^2+sc_n^2+4*sc_m*sc_n)/2/(sc_m^2+sc_m*sc_n+sc_n^2));
% show in degrees
fprintf("constructing %.2f degree TBH (%d atoms) \n", rot_theta/pi*180, total_num);

rot_mat=[cos(rot_theta),-sin(rot_theta);sin(rot_theta),cos(rot_theta)];

ra1=rot_mat*a1;
ra2=rot_mat*a2;

ra_mat=zeros(2);
ra_mat(:,1)=ra1;
ra_mat(:,2)=ra2;


% for both bottom and top unit
sc_t1=sc_n*a1+sc_m*a2;
sc_t2=-sc_m*a1+(sc_m+sc_n)*a2;


sc_ft1=[sc_t1(1),sc_t1(2),0];
sc_ft2=[sc_t2(1),sc_t2(2),0];
sc_ft3=[0,0,1];

sc_v=abs(dot(sc_ft1,cross(sc_ft2,sc_ft3)));
sc_b1=2*pi*cross(sc_ft2,sc_ft3)/sc_v;
sc_b2=2*pi*cross(sc_ft3,sc_ft1)/sc_v;
sc_b3=2*pi*cross(sc_ft1,sc_ft2)/sc_v;

sc_vec1=sc_b1(1:2);
sc_vec2=sc_b2(1:2)+sc_b1(1:2);
sc_gamma=sc_vec1*0;
sc_kpoint=(sc_vec1+sc_vec2)/3;
sc_mpoint=sc_vec1*0.5;

% primitive cell coordinate

pc_vec_a1=[a1',0];
pc_vec_a2=[a2',0];
pc_vec_a3=[0,0,1];

pc_vec_ra1=[ra1',0];
pc_vec_ra2=[ra2',0];
pc_vec_ra3=[0,0,1];


pc_vec_v=abs(dot(pc_vec_a1,cross(pc_vec_a2,pc_vec_a3)));
pc_vec_b1=2*pi*cross(pc_vec_a2,pc_vec_a3)/pc_vec_v;
pc_vec_b2=2*pi*cross(pc_vec_a3,pc_vec_a1)/pc_vec_v;
pc_vec_b3=2*pi*cross(pc_vec_a1,pc_vec_a2)/pc_vec_v;

pc_vec_rv=abs(dot(pc_vec_ra1,cross(pc_vec_ra2,pc_vec_ra3)));
pc_vec_rb1=2*pi*cross(pc_vec_ra2,pc_vec_ra3)/pc_vec_rv;
pc_vec_rb2=2*pi*cross(pc_vec_ra3,pc_vec_ra1)/pc_vec_rv;
pc_vec_rb3=2*pi*cross(pc_vec_ra1,pc_vec_ra2)/pc_vec_rv;

pc_bot_Mpoint=pc_vec_b1/2;
pc_bot_Kpoint=(2*pc_vec_b1+pc_vec_b2)/3;


pc_top_Mpoint=pc_vec_rb1/2;
pc_top_Kpoint=(2*pc_vec_rb1+pc_vec_rb2)/3;



% shift the layer

shift_top_vec=(-1/3)*(2*ra2-ra1);
shift_bot_vec=(-1/3)*(2*a2-a1);






