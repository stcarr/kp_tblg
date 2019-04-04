/*
 * File:   kp_tblg_construct.cpp
 * Author: Stephen Carr
 *
 * Created on January 13, 2019, 3:41 PM
 */

#include "kp_tblg_construct.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdexcept>
#include <math.h>     // sqrt, M_PI
#include <algorithm>  // max
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>

// --------------------------------------
// Old constructor, no use of class Sdata
// --------------------------------------
// Empty constructor
Kp_tblg_construct::Kp_tblg_construct(){

  inter_fac = 1.0;
  strain_fac = 1.0;
  full_mono_ham = 0;

}

// Copy constructor
Kp_tblg_construct::Kp_tblg_construct(const Kp_tblg_construct& orig){


}

// Destructor
Kp_tblg_construct::~Kp_tblg_construct(){


}

// constructor with angle supplied
Kp_tblg_construct::Kp_tblg_construct(double theta_in){

    theta = theta_in;
    inter_fac = 1.0;
    strain_fac = 1.0;
    full_mono_ham = 0;
}

// define twisting angle, can be changed before call to prepare
void Kp_tblg_construct::setTwist(double theta_in){

  theta = theta_in;


}

void Kp_tblg_construct::setInterFac(double interfac_in){

  inter_fac = interfac_in;

}

void Kp_tblg_construct::setStrainFac(double strainfac_in){

  strain_fac = strainfac_in;

}

void Kp_tblg_construct::setFullMonoHam(int fullmonoham_in){

  full_mono_ham = fullmonoham_in;

}


void Kp_tblg_construct::loadFiles(string filename){

  string line;
  //string filename = "full_relax_kp_01-06-2019.dat";
  ifstream dat_file (filename);

  int file_section = 0;
  int kp_term = 0;
  int shell_term = 0;

  if (dat_file.is_open())
  {
    while ( getline (dat_file,line) )
    {
      stringstream ss(line);

      switch (file_section){

        case 0: // header

          ss >> ntheta;
          ss >> n_inter_couplings;
          ss >> n_intra_couplings;
          ss >> n_inter_shells;
          ss >> n_intra_shells;

          //cout << ntheta << ", " << n_inter_couplings << ", " << n_intra_couplings << ", " << n_inter_shells << ", " << n_intra_shells << '\n';
          file_section++;
          getline (dat_file,line); // skip "theta (deg):" line
          break;

        case 1: // thetas
          thetas_deg.resize(ntheta);
          for(int i = 0; i < ntheta; ++i){
            ss >> thetas_deg[i];
          }
          file_section++;
          // skip the getline since case 2 will read per theta
          //getline (dat_file,line); // skip "inter terms:" line
          break;

        case 2: {  // kp terms
          vector< vector< vector< vector< complex<double> > > > > temp_kp_data;
          temp_kp_data.resize(ntheta);

          for (int t = 0; t < ntheta; ++t){
            getline (dat_file,line); // get data for this theta
            stringstream loop_ss(line);

            int n_couplings;

            if (kp_term < 3)
              n_couplings = n_inter_couplings;
            else
              n_couplings = n_intra_couplings;

            temp_kp_data[t].resize(n_couplings);

            for (int c = 0; c < n_couplings; ++c){
              temp_kp_data[t][c].resize(2);
              temp_kp_data[t][c][0].resize(2);
              temp_kp_data[t][c][1].resize(2);
              for (int o1 = 0; o1 < 2; ++o1){
                for (int o2 = 0; o2 < 2; ++o2){

                  string temp_string;
                  loop_ss >> temp_string;
                  //istringstream stm(temp_string);
		  double real, imag;

		  // keep track of if the complex # read works or not
		  int str_failure = 0;
		  // sign of imaginary part
		  double imag_sign = 1.0;

		  // manual read of complex #
		  size_t str_idx = temp_string.find_first_of("+-");

		  // if real part is negative, we need to use this work around
		  if(temp_string[0] == '-'){
		    str_idx = temp_string.substr(str_idx+1).find_first_of("+-") + 1;
		  }

		  if (temp_string[str_idx] == '+'){

		  }
	 	  else if (temp_string[str_idx] == '-'){
		    imag_sign = -1.0;
		  } else {
		    str_failure = 1;
		  }

		  if (str_failure == 0){
		    real = stod(temp_string);
		    imag = imag_sign*stod(temp_string.substr(str_idx+1));
		    temp_kp_data[t][c][o1][o2] = complex<double>(real,imag);
                  } else{
                    cerr << "badly formed number '" << temp_string << "'\n" ;
                    cerr << "exiting... \n";
                    return;
                  }


                }
              }
            }

          }

          if (kp_term == 0)
            All_Eff_inter = temp_kp_data;
          else if (kp_term == 1)
            All_Eff_inter_kplus = temp_kp_data;
          else if (kp_term == 2)
            All_Eff_inter_kminus = temp_kp_data;
          else if (kp_term == 3)
            All_Eff_intra_bot = temp_kp_data;
          else if (kp_term == 4) {
            All_Eff_intra_top = temp_kp_data;
            file_section++;
          }

          kp_term++;
          break;
        }

        case 3: { // shell indices
          getline (dat_file,line); // get shell size info
          stringstream temp_ss(line);

          int n_shells;
          vector< vector <vector<int> > > temp_shells;
          if (shell_term == 0) {
            n_shells = n_inter_shells;
          } else {
            n_shells = n_intra_shells;
          }

          temp_shells.resize(n_shells);

          int tot_shells = 0;

          for (int i = 0; i < n_shells; ++i){
            int sz;
            temp_ss >> sz;

            temp_shells[i].resize(sz);
            for (int j = 0; j < sz; ++j){
              temp_shells[i][j].resize(3);
            }
          }

          for (int i = 0; i < n_shells; ++i){
            int sz = temp_shells[i].size();
            for (int j = 0; j < sz; ++j){

              getline (dat_file,line); // get shell
              stringstream loop_ss(line);

              loop_ss >> temp_shells[i][j][0];
              loop_ss >> temp_shells[i][j][1];
              loop_ss >> temp_shells[i][j][2];
            }
          }

          if (shell_term == 0) {
            All_Eff_inter_shell_indices = temp_shells;
          } else {
            All_Eff_intra_shell_indices = temp_shells;
            file_section++;
          }

          shell_term++;
          break;
        }

        case 4:
          break;

      }
    }
    dat_file.close();
  }
  else cout << "Unable to open file: " << filename;

}

vector< vector< vector< complex<double> > > > Kp_tblg_construct::interpKP(double interp_scale, vector< vector< vector< complex<double> > > > high, vector< vector< vector< complex<double> > > > low){

  int s1 = high.size();
  int s2 = high[0].size();
  int s3 = high[0][0].size();

  vector< vector< vector< complex<double> > > > out;
  out.resize(s1);

  for (int i1 = 0; i1 < s1; ++i1){
    out[i1].resize(s2);
    for (int i2 = 0; i2 < s2; ++i2){
      out[i1][i2].resize(s3);
      for (int i3 = 0; i3 < s3; ++i3){
        out[i1][i2][i3] = high[i1][i2][i3]*(1.0 - interp_scale) + low[i1][i2][i3]*interp_scale;

      }
    }

  }

  return out;

}

// sets up the sparsity pattern for the Hamilotnian
void Kp_tblg_construct::prepare(){

  /* Some useful k-point info:
     our moire BZ gamma point is at [0 0]
     reciprocal vectors are given by hex_b1 and hex_b2 later in code
     kk1a/b/c give the K points, kk4 gives M point (depend on hex_b's)
     search for kk1a or kk4 in code to find their definition
  */

  // Currently only fully relaxed data available for this method

  // set a cut-off condition for the k-p model
  // that scales automatically with theta
  double hex_cut_fac0 = 5.0;
  double theta_cut_0 = 1.5;


  // keep all inter terms for now
  int max_inter_q = n_inter_couplings;
  int max_intra_q = n_intra_couplings;

  int max_inter_q_plusminus = 3; // only first k+- correction for interlayer
  int max_intra_q_plusminus = 0; // no k+- correction for intralayer

  // cut-off radius for our kp expansion
  double hex_cut_fac = max(hex_cut_fac0*sqrt(theta_cut_0/theta), hex_cut_fac0);

  // interpolate the kp-model parameters
  int min_theta_idx = -1;
  for (int i = 0; i < ntheta; ++i){
    if (theta - thetas_deg[i] < 0){
      min_theta_idx = i;
    }
  }

  if (min_theta_idx == -1 || min_theta_idx == ntheta-1){
    printf("ERROR! Selected theta outside of fitted data range: [%lf, %lf] \n",thetas_deg[0],thetas_deg[ntheta-1]);
    return;
  }

  double theta_high = thetas_deg[min_theta_idx];
  double theta_low = thetas_deg[min_theta_idx+1];
  double interp_scale = (theta_high - theta)/(theta_high - theta_low);

  //printf("interp: [%lf - %lf], %lf \n",theta_high, theta_low, interp_scale);

  vector< vector< vector<complex<double> > > > inter_k0  = interpKP(interp_scale, All_Eff_inter[min_theta_idx],         All_Eff_inter[min_theta_idx+1]);
  vector< vector< vector<complex<double> > > > inter_kp  = interpKP(interp_scale, All_Eff_inter_kplus[min_theta_idx],   All_Eff_inter_kplus[min_theta_idx+1]);
  vector< vector< vector<complex<double> > > > inter_km  = interpKP(interp_scale, All_Eff_inter_kminus[min_theta_idx],  All_Eff_inter_kminus[min_theta_idx+1]);
  vector< vector< vector<complex<double> > > > intra_bot = interpKP(interp_scale, All_Eff_intra_bot[min_theta_idx],     All_Eff_intra_bot[min_theta_idx+1]);
  vector< vector< vector<complex<double> > > > intra_top = interpKP(interp_scale, All_Eff_intra_top[min_theta_idx],     All_Eff_intra_top[min_theta_idx+1]);

  //printf("|w0| = %lf, |w1| = %lf \n", abs(inter_k0[0][0][0]), abs(inter_k0[0][0][1]));

  double rot_theta = theta*M_PI/180;
  double lattice_a = 1.42*sqrt(3);

  // number of DoF per k-point (for a single layer)
  unit_dim = 2;

  // sets the linear and quadratic term of the Dirac cones
  Dirac_v1 = 5.2268; // v_F
  Dirac_v2 = 2.2450*0; // no quadratic for now
  Dirac_diag = -1.1142*1.0; // linear diagonal term

  // some important length scales for the Brill. Zone
  double KD = 4*M_PI/3/lattice_a;
  double KTH = 2*KD*sin(rot_theta/2.0);

  double HEX_BLEN=KTH*sqrt(3);

  // basic pauli matrices
  // Not used??

  /*
  Eigen::Matrix2cd sigx;

  sigx(0,0) = 0;
  sigx(1,0) = 1.0;
  sigx(0,1) = 1.0;
  sigx(1,1) = 0;

  Eigen::Matrix2cd sigy;
  sigy(0,0) = 0;
  sigy(1,0) = 1i;
  sigy(0,1) = -1i;
  sigy(1,1) = 0;

  Eigen::Matrix2cd sigz;
  sigz(0,0) = 1.0;
  sigz(1,0) = 0;
  sigz(0,1) = 0;
  sigz(1,1) = -1.0;

  Eigen::Matrix2cd umat = expm(-1i*rot_theta*sigz/2);

  rsigx = umat.adjoint()*sigx*umat;
  rsigy = umat.adjoint()*sigy*umat;
  rsigz = umat.adjoint()*sigz*umat;
  */

  Vector2d hex_b1(HEX_BLEN*sqrt(3.0)/2.0, -HEX_BLEN*1.0/2.0);
  Vector2d hex_b2(HEX_BLEN*0.0, HEX_BLEN*1.0);

  Vector2d hex_shift = (-hex_b1+hex_b2)/3.0;

  // create hex table of possible intralayer interactions

  vector< vector<int> > intra_qs;
  vector< vector<int> > inter_qs;
  vector< Vector2d > intra_qs_vec;
  vector< Vector2d > interall_given_qs;

  // vars for intra_qs momenta calc
  Vector2d miller_b1 = hex_b1;
  Vector2d miller_b2 = hex_b1 + hex_b2;

  // vars for inter_qs momenta calc
  double a_rot = 120.0*M_PI/180.0;

  // !! might need to fix
  // r_mat = [cosd(a_rot) sind(a_rot); -sind(a_rot) cosd(a_rot)];
  Matrix2d r_mat;
  r_mat(0,0) = cos(a_rot);
  r_mat(0,1) = sin(a_rot);
  r_mat(1,0) = -sin(a_rot);
  r_mat(1,1) = cos(a_rot);
  Vector2d a1 = -2.0*hex_shift+hex_b2;
  Vector2d a2 = r_mat*a1;
  Vector2d a3 = r_mat*a2;

  for (int shell_idx = 0; shell_idx < All_Eff_intra_shell_indices.size(); ++shell_idx){
      vector< vector<int> > shell_here = All_Eff_intra_shell_indices[shell_idx];
      int n_momenta = shell_here.size();
      for (int n = 0; n < n_momenta; ++n){
          intra_qs.push_back(shell_here[n]);
          Vector2d q_here = ((double)shell_here[n][0])*miller_b1 + ((double)shell_here[n][1])*miller_b2;
          intra_qs_vec.push_back(q_here);
      }

  }

  // create hex table of possible interlayer interactions
  //int num_inter_qs = inter_k0.size();


  for (int shell_idx = 0; shell_idx < All_Eff_inter_shell_indices.size(); ++shell_idx){
      vector< vector<int> > shell_here = All_Eff_inter_shell_indices[shell_idx];
      int n_momenta = shell_here.size();
      for (int n = 0; n < n_momenta; ++n){
          inter_qs.push_back(shell_here[n]);
          Vector2d q_here = shell_here[n][0]*a1 + shell_here[n][1]*a2 + shell_here[n][2]*a3;
          interall_given_qs.push_back(q_here);

      }

  }

  vector<Vector2d> L12_qvecs = interall_given_qs;
  // loop over a large grid of hex points to find all valid couplings
  int hex_M = 100;

  int hex_index = 0;
  vector< Vector2d > hex_coor;
  double hex_cut= hex_cut_fac*HEX_BLEN; // 3.21*HEX_BLEN;

  G_to_index = MatrixXi::Constant(2*hex_M+1, 2*hex_M+1,-1);

  int ind = 0;
  for (int ind1 = -hex_M; ind1 < hex_M+1; ++ind1) {
      for (int ind2 = -hex_M; ind2 < hex_M+1; ++ind2) {

          Vector2d vec = hex_b1*ind1 + hex_b2*ind2 + hex_shift;

          if (sqrt(vec.dot(vec)) < hex_cut ){
              hex_coor.push_back(vec);

              G_to_index(ind1+hex_M,ind2+hex_M) = ind;
              ind++;
          }

      }
  }

  num_hex = hex_coor.size();

  index_to_G = MatrixXi::Constant(num_hex, 2,-1);

  for (int ind1 = -hex_M; ind1 < hex_M+1; ++ind1) {
      for (int ind2 = -hex_M; ind2 < hex_M+1; ++ind2) {
           int idx_here = G_to_index(ind1+hex_M,ind2+hex_M);
           if (idx_here != -1){
             index_to_G(idx_here,0) = ind1;
             index_to_G(idx_here,1) = ind2;

           }
       }
  }

  Vector3d moire_k_vec1(hex_b1[0], hex_b1[1], 0.0);
  Vector3d moire_k_vec2(hex_b2[0], hex_b2[1], 0.0);
  Vector3d moire_k_vec3(0.0,0.0,1.0);

  double vvv = abs(moire_k_vec1.dot( moire_k_vec2.cross(moire_k_vec3) ));
  Vector3d moire_L_x1 = 2.0*M_PI*moire_k_vec2.cross(moire_k_vec3)/vvv;
  Vector3d moire_L_x2 = 2.0*M_PI*moire_k_vec3.cross(moire_k_vec1)/vvv;
  Vector3d moire_L_x3 = 2.0*M_PI*moire_k_vec1.cross(moire_k_vec2)/vvv;

  // total number of DoFs for the KP hamiltonian
  // number of kpoints * number of DoF per k point * number of layers
  tot_dim = num_hex*unit_dim*2;


  hex_all_L1 = hex_coor;
  for (int j = 0; j < hex_coor.size(); ++j){
    hex_all_L2.push_back(-hex_coor[j]);
  }

  double vecthres = 1E-5;
  // create connection matrices

  int num_intra_qs = max_intra_q;
  int num_inter_qs = max_inter_q;

  int tot_num_G12 = num_inter_qs;

  // easier to make these complex<double> than deal with casting during H construction
  vector< vector< vector< complex<double> > > > connect_Mat_L1, connect_Mat_L2;
  vector< vector< vector< complex<double> > > > connect_Mat_L12;

  connect_Mat_L1.resize(num_hex);
  connect_Mat_L2.resize(num_hex);
  connect_Mat_L12.resize(num_hex);

  for (int i1 = 0; i1 < num_hex; ++i1){
    connect_Mat_L1[i1].resize(num_hex);
    connect_Mat_L2[i1].resize(num_hex);
    connect_Mat_L12[i1].resize(num_hex);
    for (int i2 = 0; i2 < num_hex; ++i2){

      connect_Mat_L1[i1][i2].resize(num_intra_qs);
      connect_Mat_L2[i1][i2].resize(num_intra_qs);

      for (int mono_i3 = 0; mono_i3 < num_intra_qs; ++mono_i3){

        connect_Mat_L1[i1][i2][mono_i3] = 0.0;
        connect_Mat_L2[i1][i2][mono_i3] = 0.0;
      }

      connect_Mat_L12[i1][i2].resize(num_inter_qs);

      for (int bi_i3 = 0; bi_i3 < num_inter_qs; ++bi_i3){
        connect_Mat_L12[i1][i2][bi_i3] = 0.0;
      }

    }
  }

  for (int ind1 = 0; ind1 < num_hex; ++ind1) {
    Vector2d pvec1_L1 = hex_all_L1[ind1];
    Vector2d pvec1_L2 = hex_all_L2[ind1];

      for (int ind2 = 0; ind2 < num_hex; ++ind2){

        Vector2d pvec2_L1 = hex_all_L1[ind2];
        Vector2d pvec2_L2 = hex_all_L2[ind2];

        //for (int indty = 0; indty < 3; ++indty){
        for (int indty = 0; indty < num_inter_qs; ++indty){
          Vector2d tmp_qvec = L12_qvecs[indty];
          Vector2d dvec_L12 = pvec2_L2 - pvec1_L1 - tmp_qvec;

          if (sqrt(dvec_L12.dot(dvec_L12)) < vecthres){
              connect_Mat_L12[ind2][ind1][indty] = 1.0;
              // if (ind2 == 0)
              // cout << "L_12 coupling: [" << ind1 << ", " << ind2 << ", " << indty << "] \n";
              // cout << "L1 G: \n" << index_to_G(ind1,0) << " " << index_to_G(ind1,1) << "\n";
              // cout << "L2 G: \n" << index_to_G(ind2,0) << " " << index_to_G(ind2,1) << "\n";
              // cout << "tmp_qvec: \n" << tmp_qvec << "\n";

            }
        }


        for (int indty = 0; indty < num_intra_qs; ++indty){

            Vector2d tmp_qvec = intra_qs_vec[indty];

            Vector2d dvec_L1 = pvec2_L1 - pvec1_L1 - tmp_qvec;
            Vector2d dvec_L2 = pvec2_L2 - pvec1_L2 - tmp_qvec;
            // dvec=pvec2-pvec1-(tmp_set(1)*sc_b1(1:2)+tmp_set(2)*sc_b2(1:2));

            if (sqrt(dvec_L1.dot(dvec_L1)) < vecthres)
                connect_Mat_L1[ind1][ind2][indty] = 1.0;

            if (sqrt(dvec_L2.dot(dvec_L2)) < vecthres)
                connect_Mat_L2[ind1][ind2][indty] = 1.0;

        }


    }
  }

  // for now these are just = 1, adds a local phase rotaion to each layer but doesn't seem needed
  expfac1 = exp(complex<double>(0.0, rot_theta/2.0));
  expfac2 = exp(complex<double>(0.0,-rot_theta/2.0));

  all_index_L1.resize(unit_dim*num_hex);
  all_index_L2.resize(unit_dim*num_hex);
  all_index.resize(2*unit_dim*num_hex);

  for (int i = 0; i < unit_dim*num_hex; ++i){
    all_index_L1[i] = i;
    all_index_L2[i] = i + unit_dim*num_hex;
  }

  for (int i = 0; i < 2*unit_dim*num_hex; ++i){
    all_index[i] = i;
  }

  /*
  all_index_L1 = reshape(1:(num_hex*unit_dim),unit_dim,num_hex);
  all_index_L2 = reshape((num_hex*unit_dim+1):(2*num_hex*unit_dim),unit_dim,num_hex);
  all_index = reshape(1:tot_dim,unit_dim,tot_dim/unit_dim);
  */

  Hmat_strain_L1 = MatrixXcd::Zero(tot_dim,tot_dim);
  Hmat_strain_L2 = MatrixXcd::Zero(tot_dim,tot_dim);

  /*
  MatrixXcd Hmat_strain_L1_kplus = MatrixXcd::Zero(tot_dim,tot_dim);
  MatrixXcd Hmat_strain_L2_kplus = MatrixXcd::Zero(tot_dim,tot_dim);
  MatrixXcd Hmat_strain_L1_kminus = MatrixXcd::Zero(tot_dim,tot_dim);
  MatrixXcd Hmat_strain_L2_kminus = MatrixXcd::Zero(tot_dim,tot_dim);
  */

  for (int indq = 0; indq < num_intra_qs; ++indq){

      vector< vector< complex<double> > > T_mat_L1 = intra_bot[indq];
      vector< vector< complex<double> > > T_mat_L2 = intra_top[indq];
      for (int r = 0; r < num_hex; ++r){
        for (int c = 0; c < num_hex; ++c){
          for (int o1 = 0; o1 < unit_dim; ++o1){
            for (int o2 = 0; o2 < unit_dim; ++o2){
              Hmat_strain_L1(all_index_L1[unit_dim*r+o1],all_index_L1[unit_dim*c+o2]) += connect_Mat_L1[r][c][indq]*T_mat_L1[o1][o2];
            }
          }

        }
      }
      //Hmat_strain_L1(all_index_L1(:),all_index_L1(:))=Hmat_strain_L1(all_index_L1(:),all_index_L1(:))+kron(squeeze(connect_Mat_L1(:,:,indq)),T_mat_L1);
      for (int r = 0; r < num_hex; ++r){
        for (int c = 0; c < num_hex; ++c){
          for (int o1 = 0; o1 < unit_dim; ++o1){
            for (int o2 = 0; o2 < unit_dim; ++o2){
              Hmat_strain_L2(all_index_L2[unit_dim*r+o1],all_index_L2[unit_dim*c+o2]) += connect_Mat_L2[r][c][indq]*T_mat_L2[o1][o2];
            }
          }

        }
      }
      //Hmat_strain_L2(all_index_L2(:),all_index_L2(:))=Hmat_strain_L2(all_index_L2(:),all_index_L2(:))+kron(squeeze(connect_Mat_L2(:,:,indq)),T_mat_L2);

      /*
      T_mat_L1_kplus=squeeze(All_Eff_intra_bot_kplus(:,:,indq));
      T_mat_L2_kplus=squeeze(All_Eff_intra_top_kplus(:,:,indq));
      Hmat_strain_L1_kplus(all_index_L1(:),all_index_L1(:))=Hmat_strain_L1_kplus(all_index_L1(:),all_index_L1(:))+kron(squeeze(connect_Mat_L1(:,:,indq)),T_mat_L1_kplus);
      Hmat_strain_L2_kplus(all_index_L2(:),all_index_L2(:))=Hmat_strain_L2_kplus(all_index_L2(:),all_index_L2(:))+kron(squeeze(connect_Mat_L2(:,:,indq)),T_mat_L2_kplus);

      T_mat_L1_kminus=squeeze(All_Eff_intra_bot_kminus(:,:,indq));
      T_mat_L2_kminus=squeeze(All_Eff_intra_top_kminus(:,:,indq));
      Hmat_strain_L1_kminus(all_index_L1(:),all_index_L1(:))=Hmat_strain_L1_kminus(all_index_L1(:),all_index_L1(:))+kron(squeeze(connect_Mat_L1(:,:,indq)),T_mat_L1_kminus);
      Hmat_strain_L2_kminus(all_index_L2(:),all_index_L2(:))=Hmat_strain_L2_kminus(all_index_L2(:),all_index_L2(:))+kron(squeeze(connect_Mat_L2(:,:,indq)),T_mat_L2_kminus);
      */
  }

  // make sure it is Hermitian (should not be necessary though!)
  Hmat_strain_L1 = (Hmat_strain_L1 + Hmat_strain_L1.adjoint().eval())/2.0;
  Hmat_strain_L2 = (Hmat_strain_L2 + Hmat_strain_L2.adjoint().eval())/2.0;

  MatrixXcd Hmat_inter = MatrixXcd::Zero(tot_dim,tot_dim);
  Hmat_inter_kplus0 = MatrixXcd::Zero(tot_dim/2,tot_dim/2);
  Hmat_inter_kminus0 = MatrixXcd::Zero(tot_dim/2,tot_dim/2);

  for (int indq = 0; indq < max_inter_q; ++indq){

      vector< vector< complex<double> > > T_tmp = inter_k0[indq];
      //cout << "inter_k0[" << indq << "] = [" << T_tmp[0][0] << ", " << T_tmp[0][1] << " | " << T_tmp[1][0] << ", " << T_tmp[1][1] << "] \n";

      for (int r = 0; r < num_hex; ++r){
        for (int c = 0; c < num_hex; ++c){
          for (int o1 = 0; o1 < unit_dim; ++o1){
            for (int o2 = 0; o2 < unit_dim; ++o2){
              Hmat_inter(all_index_L2[unit_dim*r+o1],all_index_L1[unit_dim*c+o2]) += connect_Mat_L12[r][c][indq]*T_tmp[o1][o2];
              //if (abs(connect_Mat_L12[r][c][indq]) > 0.0 && r == 0)
                //cout << connect_Mat_L12[r][c][indq]*T_tmp[o1][o2] << "added to Hmat_inter for " << r << " " << num_hex - c << " \n";
            }
          }

        }
      }
      //Hmat_inter(all_index_L2(:),all_index_L1(:))=Hmat_inter(all_index_L2(:),all_index_L1(:))+kron(squeeze(connect_Mat_L12(:,:,indq)),T_tmp);

      if (indq < max_inter_q_plusminus) { // only keep k+- for nearest inter terms
        vector< vector< complex<double> > > T_tmp_kplus = inter_kp[indq];
        for (int r = 0; r < num_hex; ++r){
          for (int c = 0; c < num_hex; ++c){
            for (int o1 = 0; o1 < unit_dim; ++o1){
              for (int o2 = 0; o2 < unit_dim; ++o2){
                Hmat_inter_kplus0(unit_dim*r+o1,unit_dim*c+o2) += connect_Mat_L12[r][c][indq]*T_tmp_kplus[o1][o2]*lattice_a;
              }
            }

          }
        }
        //Hmat_inter_kplus0(:,:)=Hmat_inter_kplus0(:,:)+kron(squeeze(connect_Mat_L12(:,:,indq)),T_tmp_kplus)*lattice_a;

        vector< vector< complex<double> > > T_tmp_kminus = inter_km[indq];
        for (int r = 0; r < num_hex; ++r){
          for (int c = 0; c < num_hex; ++c){
            for (int o1 = 0; o1 < unit_dim; ++o1){
              for (int o2 = 0; o2 < unit_dim; ++o2){
                Hmat_inter_kminus0(unit_dim*r+o1,unit_dim*c+o2) += connect_Mat_L12[r][c][indq]*T_tmp_kminus[o1][o2]*lattice_a;
              }
            }
          }
        }
        //Hmat_inter_kminus0(:,:)=Hmat_inter_kminus0(:,:)+kron(squeeze(connect_Mat_L12(:,:,indq)),T_tmp_kminus)*lattice_a;
      }

  }

  MatrixXcd Hmat_inter_adj = Hmat_inter.adjoint();

  //Hmat_inter = Hmat_inter ++ Hmat_inter.adjoint().eval();
  Hmat_inter = Hmat_inter + Hmat_inter_adj;// + Hmat_inter.adjoint().eval();

  MatrixXcd Hmat_inter_qplus0 = MatrixXcd::Zero(tot_dim,tot_dim);
  MatrixXcd Hmat_inter_qminus0 = MatrixXcd::Zero(tot_dim,tot_dim);

  for (int indq = 0; indq < max_inter_q_plusminus; ++indq){

      Vector2d given_q = interall_given_qs[indq];
      complex<double> given_qplus = complex<double>(given_q[0],given_q[1]);
      complex<double> given_qminus = complex<double>(given_q[0],-1.0*given_q[1]);

      vector< vector< complex<double> > > T_tmp_kplus = inter_kp[indq];

      for (int r = 0; r < num_hex; ++r){
        for (int c = 0; c < num_hex; ++c){
          for (int o1 = 0; o1 < unit_dim; ++o1){
            for (int o2 = 0; o2 < unit_dim; ++o2){
              Hmat_inter_qplus0(all_index_L2[unit_dim*r+o1],all_index_L1[unit_dim*c+o2]) += connect_Mat_L12[r][c][indq]*T_tmp_kplus[o1][o2]*lattice_a*given_qplus/2.0;
            }
          }

        }
      }
      // Hmat_inter_qplus0(all_index_L2(:),all_index_L1(:)) = Hmat_inter_qplus0(all_index_L2(:),all_index_L1(:))+kron(squeeze(connect_Mat_L12(:,:,indq)),T_tmp_kplus)*lattice_a*given_qplus/2;


      vector< vector< complex<double> > > T_tmp_kminus= inter_km[indq];

      for (int r = 0; r < num_hex; ++r){
        for (int c = 0; c < num_hex; ++c){
          for (int o1 = 0; o1 < unit_dim; ++o1){
            for (int o2 = 0; o2 < unit_dim; ++o2){
              Hmat_inter_qminus0(all_index_L2[unit_dim*r+o1],all_index_L1[unit_dim*c+o2]) += connect_Mat_L12[r][c][indq]*T_tmp_kminus[o1][o2]*lattice_a*given_qminus/2.0;
            }
          }

        }
      }
      //Hmat_inter_qminus0(all_index_L2(:),all_index_L1(:)) = Hmat_inter_qminus0(all_index_L2(:),all_index_L1(:))+kron(squeeze(connect_Mat_L12(:,:,indq)),T_tmp_kminus)*given_qminus/2;

  }

  double q_dep_fac = 1.0;
  Hmat_inter_qdep = Hmat_inter + q_dep_fac*(Hmat_inter_qplus0 + Hmat_inter_qminus0 + Hmat_inter_qplus0.adjoint().eval() + Hmat_inter_qminus0.adjoint().eval());


}

Vector2d Kp_tblg_construct::getK(){

  double rot_theta = theta*M_PI/180.0;
  double lattice_a = 1.42*sqrt(3.0);

  double KD = 4*M_PI/(3.0*lattice_a);
  double KTH = 2*KD*sin(rot_theta/2.0);
  double HEX_BLEN = KTH*sqrt(3.0);

  Vector2d hex_b1(HEX_BLEN*sqrt(3.0)/2.0,HEX_BLEN*-1.0/2.0);
  Vector2d hex_b2(0.0,HEX_BLEN*1.0);

  Vector2d K_pt = (hex_b1 - hex_b2)/3.0;
  return K_pt;

}

Vector2d Kp_tblg_construct::getM(){

  double rot_theta = theta*M_PI/180.0;
  double lattice_a = 1.42*sqrt(3.0);

  double KD = 4*M_PI/(3.0*lattice_a);
  double KTH = 2*KD*sin(rot_theta/2.0);
  double HEX_BLEN = KTH*sqrt(3.0);

  Vector2d hex_b1(HEX_BLEN*sqrt(3.0)/2.0,HEX_BLEN*-1.0/2.0);
  Vector2d hex_b2(0.0,HEX_BLEN*1.0);

  Vector2d M_pt = hex_b1/2.0;
  return M_pt;

}

Vector2d Kp_tblg_construct::getGamma(){

  Vector2d Gamma_pt(0.0,0.0);
  return Gamma_pt;

}

Matrix2cd Kp_tblg_construct::layer1Ham(Vector2d kknow){

  // the terms related to expfac2 need to be checked if quadratic Dirac terms nearest turned on
  // (might be misusing expfac1 vs expfac2)
  Matrix2cd out = Matrix2cd::Zero(2,2);

  // Dirac Cone approximation
  if (full_mono_ham == 0){
    out(0,0) = Dirac_diag*kknow.dot(kknow);
    out(0,1) = expfac1*Dirac_v1*complex<double>(kknow(1),-1.0*kknow(0))+pow(expfac1,2.0)*Dirac_v2*pow(complex<double>(kknow(0),-1.0*kknow(1)),2.0);
    out(1,0) = conj(out(0,1));
    out(1,1) = out(0,0);
    return out;
  }
  /* =
  (1,1) : [Dirac_diag*dot(kknow,kknow),
  (1,2) : expfac1*Dirac_v1*(-i*kknow(1)+kknow(2))+(expfac2^2)*Dirac_v2*(kknow(1)-i*kknow(2))^2;
  (2,1) : conj(expfac1)*Dirac_v1*(i*kknow(1)+kknow(2))+(expfac1^2)*Dirac_v2*(kknow(1)+i*kknow(2))^2,
  (2,2) : Dirac_diag*dot(kknow,kknow)];
  */


  // Otherwise do Full Hamiltonian

  double lattice_a = 1.42*sqrt(3.0);

  vector< vector<double> > a_base; // unrotated monolayer unit-cell
  a_base.resize(2);
  a_base[0].resize(2);
  a_base[0][0] =  lattice_a*sqrt(3.0)/2.0;
  a_base[0][1] = -lattice_a*1.0/2.0;
  a_base[1].resize(2);
  a_base[1][0] =  lattice_a*sqrt(3.0)/2.0;
  a_base[1][1] =  lattice_a*1.0/2.0;

  vector< vector<double> > a; // rotated unit-cell
  double theta_rad = -(theta/2.0)*(M_PI/180.0);
  a.resize(2);
  a[0].resize(2);
  a[0][0] = cos(theta_rad)*a_base[0][0] - sin(theta_rad)*a_base[0][1];
  a[0][1] = sin(theta_rad)*a_base[0][0] + cos(theta_rad)*a_base[0][1];
  a[1].resize(2);
  a[1][0] = cos(theta_rad)*a_base[1][0] - sin(theta_rad)*a_base[1][1];
  a[1][1] = sin(theta_rad)*a_base[1][0] + cos(theta_rad)*a_base[1][1];


  vector< vector<double> > b1 = getReciprocal(a); // reciprocal of monolayer unit-cell

	//printf("b1 = [%lf %lf; %lf %lf] \n",b1[0][0],b1[0][1],b1[1][0],b1[1][1]);

  // High Symm points of the monolayer BZ
	double gamma[2];
	double m[2];
	double k[2];

	gamma[0] = 0.0;
	gamma[1] = 0.0;
	m[0] = b1[1][0]/2.0;
	m[1] = b1[1][1]/2.0;

  double phi = M_PI/6;
  k[0] = (cos(phi)*m[0] - sin(phi)*m[1])/cos(phi);
  k[1] = (sin(phi)*m[0] + cos(phi)*m[1])/cos(phi);

  //k[0] = (1.0/(2.0*cos(M_PI/6)))*(cos(M_PI/2)*b1[1][0] + sin(M_PI/2)*b1[1][1]);
  //k[1] = (1.0/(2.0*cos(M_PI/6)))*(-1.0*sin(M_PI/2)*b1[1][0] + cos(M_PI/2)*b1[1][1]);

  // shift to unit-cell cone
  double q1 = kknow(0) + k[0];
  double q2 = kknow(1) + k[1];

  vector< vector<double> > orb_pos;
  orb_pos.resize(2);
  orb_pos[0].resize(2);
  orb_pos[0][0] = 0.0;
  orb_pos[0][1] = 0.0;
  orb_pos[1].resize(2);
  orb_pos[1][0] = (1.0/3.0)*(a[0][0] + a[1][0]);
  orb_pos[1][1] = (1.0/3.0)*(a[0][1] + a[1][1]);

  int L = 5;

	for (int i = -L; i < L+1; ++i){
		for (int j = -L; j < L+1; ++j){


      /*
			// unrotated variables
			double temp_x = i*a[0][0] + j*a[1][0];
			double temp_y = i*a[0][1] + j*a[1][1];
			double temp_z = 0;

			// rotated variables
			double x = temp_x*cos(theta) - temp_y*sin(theta);
			double y = temp_x*sin(theta) + temp_y*cos(theta);
			double z = temp_z;
      */

      double x = i*a[0][0] + j*a[1][0];
			double y = i*a[0][1] + j*a[1][1];

			array<int,2> grid_disp = {{ i, j }};

			vector< vector<double> > temp_o1_bloch_t_array;
			for (int o1 = 0; o1 < 2; ++o1){
				vector<double> temp_o2_bloch_t_array;
				for (int o2 = 0; o2 < 2; ++o2){

          double x_h = x + orb_pos[o2][0] - orb_pos[o1][0];
          double y_h = y + orb_pos[o2][1] - orb_pos[o1][1];

					// compute R(i,j) based on <a> and <theta>
					// find monolayer t_ij(R(i,j),o1,o2)
					// if t_ij != 0, add R and t_ij to respective vectors, and ++N_R;
					double t = grapheneIntralayerTerm(o1, o2, grid_disp);
          double phase = (q1*x_h + q2*y_h);
          out(o1,o2) = out(o1,o2) + polar(t,phase);

				} // end of o2 loop
			} // end of o1 loop
		} // end of j loop
	} // end of i loop

  return out;
}

Matrix2cd Kp_tblg_construct::layer2Ham(Vector2d kknow){

  // this line returns same effective monolayer hamiltonian as layer 1
  //return layer1Ham(kknow);

  // the terms related to expfac2 need to be checked if quadratic Dirac terms nearest turned on
  // (might be missusing expfac1 vs expfac2)
  Matrix2cd out = Matrix2cd::Zero(2,2);

  // Dirac Cone approximation
  if (full_mono_ham == 0){
    out(0,0) = Dirac_diag*kknow.dot(kknow);
    out(0,1) = expfac2*Dirac_v1*complex<double>(kknow(1),-1.0*kknow(0))+pow(expfac2,2.0)*Dirac_v2*pow(complex<double>(kknow(0),-1.0*kknow(1)),2.0);
    out(1,0) = conj(out(0,1));
    out(1,1) = out(0,0);
    return out;
  }
  /* =
  (1,1) : [Dirac_diag*dot(kknow,kknow),
  (1,2) : expfac2*Dirac_v1*(-i*kknow(1)+kknow(2))+(expfac2^2)*Dirac_v2*(kknow(1)-i*kknow(2))^2;
  (2,1) : conj(expfac2)*Dirac_v1*(i*kknow(1)+kknow(2))+(expfac1^2)*Dirac_v2*(kknow(1)+i*kknow(2))^2,
  (2,2) : Dirac_diag*dot(kknow,kknow)];
  */

  // Otherwise do full Hamiltonian

  double lattice_a = 1.42*sqrt(3.0);

  vector< vector<double> > a_base; // unrotated monolayer unit-cell
  a_base.resize(2);
  a_base[0].resize(2);
  a_base[0][0] =  lattice_a*sqrt(3.0)/2.0;
  a_base[0][1] = -lattice_a*1.0/2.0;
  a_base[1].resize(2);
  a_base[1][0] =  lattice_a*sqrt(3.0)/2.0;
  a_base[1][1] =  lattice_a*1.0/2.0;

  vector< vector<double> > a; // rotated unit-cell
  double theta_rad = (theta/2.0)*(M_PI/180.0);
  a.resize(2);
  a[0].resize(2);
  a[0][0] = cos(theta_rad)*a_base[0][0] - sin(theta_rad)*a_base[0][1];
  a[0][1] = sin(theta_rad)*a_base[0][0] + cos(theta_rad)*a_base[0][1];
  a[1].resize(2);
  a[1][0] = cos(theta_rad)*a_base[1][0] - sin(theta_rad)*a_base[1][1];
  a[1][1] = sin(theta_rad)*a_base[1][0] + cos(theta_rad)*a_base[1][1];

  vector< vector<double> > b1 = getReciprocal(a); // reciprocal of monolayer unit-cell

	//printf("b1 = [%lf %lf; %lf %lf] \n",b1[0][0],b1[0][1],b1[1][0],b1[1][1]);

  // High Symm points of the monolayer BZ
	double gamma[2];
	double m[2];
	double k[2];

	gamma[0] = 0.0;
	gamma[1] = 0.0;
	m[0] = b1[1][0]/2.0;
	m[1] = b1[1][1]/2.0;

  double phi = M_PI/6;
  k[0] = (cos(phi)*m[0] - sin(phi)*m[1])/cos(phi);
  k[1] = (sin(phi)*m[0] + cos(phi)*m[1])/cos(phi);

  //k[0] = (1.0/(2.0*cos(M_PI/6)))*(cos(M_PI/2)*b1[1][0] + sin(M_PI/2)*b1[1][1]);
  //k[1] = (1.0/(2.0*cos(M_PI/6)))*(-1.0*sin(M_PI/2)*b1[1][0] + cos(M_PI/2)*b1[1][1]);

  // shift to unit-cell cone
  double q1 = kknow(0) + k[0];
  double q2 = kknow(1) + k[1];

  vector< vector<double> > orb_pos;
  orb_pos.resize(2);
  orb_pos[0].resize(2);
  orb_pos[0][0] = 0.0;
  orb_pos[0][1] = 0.0;
  orb_pos[1].resize(2);
  orb_pos[1][0] = (1.0/3.0)*(a[0][0] + a[1][0]);
  orb_pos[1][1] = (1.0/3.0)*(a[0][1] + a[1][1]);

  int L = 5;

	for (int i = -L; i < L+1; ++i){
		for (int j = -L; j < L+1; ++j){


      /*
			// unrotated variables
			double temp_x = i*a[0][0] + j*a[1][0];
			double temp_y = i*a[0][1] + j*a[1][1];
			double temp_z = 0;

			// rotated variables
			double x = temp_x*cos(theta) - temp_y*sin(theta);
			double y = temp_x*sin(theta) + temp_y*cos(theta);
			double z = temp_z;
      */

      double x = i*a[0][0] + j*a[1][0];
			double y = i*a[0][1] + j*a[1][1];

			array<int,2> grid_disp = {{ i, j }};

			vector< vector<double> > temp_o1_bloch_t_array;
			for (int o1 = 0; o1 < 2; ++o1){
				vector<double> temp_o2_bloch_t_array;
				for (int o2 = 0; o2 < 2; ++o2){

          double x_h = x + orb_pos[o2][0] - orb_pos[o1][0];
          double y_h = y + orb_pos[o2][1] - orb_pos[o1][1];

					// compute R(i,j) based on <a> and <theta>
					// find monolayer t_ij(R(i,j),o1,o2)
					// if t_ij != 0, add R and t_ij to respective vectors, and ++N_R;
					double t = grapheneIntralayerTerm(o1, o2, grid_disp);
          double phase = (q1*x_h + q2*y_h);
          out(o1,o2) = out(o1,o2) + polar(t,phase);

				} // end of o2 loop
			} // end of o1 loop
		} // end of j loop
	} // end of i loop

  return out;

}

Matrix2cd Kp_tblg_construct::layer1GradHam(Vector2d kknow, int dim){

  if (dim != 0 && dim != 1){
    throw "ERROR: Invalid dimension for layer1GradHam!";
  }

  if (full_mono_ham != 0) {
    throw "ERROR: layer1GradHam only implemented for Dirac cone approximation!";
  }

  Matrix2cd out = Matrix2cd::Zero(2,2);

  /*
  out(0,0) = Dirac_diag*kknow.dot(kknow);
  out(0,1) = expfac1*Dirac_v1*complex<double>(kknow(1),-1.0*kknow(0))+pow(expfac2,2.0)*Dirac_v2*pow(complex<double>(kknow(0),-1.0*kknow(1)),2.0);
  out(1,0) = conj(out(0,1));
  out(1,1) = out(0,0);
  */

  // Dirac Cone approximation
  // if (full_mono_ham == 0){
    out(0,0) = Dirac_diag*2.0*kknow[dim];
    if (dim == 0){
      out(0,1) = expfac1*Dirac_v1*complex<double>(0.0,-1.0) + pow(expfac1,2.0)*Dirac_v2*complex<double>(2.0, 0.0)*complex<double>(kknow(0),-1.0*kknow(1));
    } else if (dim == 1){
      out(0,1) = expfac1*Dirac_v1*complex<double>(1.0, 0.0) + pow(expfac1,2.0)*Dirac_v2*complex<double>(0.0,-2.0)*complex<double>(kknow(0),-1.0*kknow(1));
    }
    out(1,0) = conj(out(0,1));
    out(1,1) = out(0,0);
    return out;
  // }

}

Matrix2cd Kp_tblg_construct::layer2GradHam(Vector2d kknow, int dim){

  if (dim != 0 && dim != 1){
    throw "ERROR: Invalid dimension for layer1GradHam!";
  }

  if (full_mono_ham != 0) {
    throw "ERROR: layer1GradHam only implemented for Dirac cone approximation!";
  }

  Matrix2cd out = Matrix2cd::Zero(2,2);


  /*  out(0,0) = Dirac_diag*kknow.dot(kknow);
      out(0,1) = expfac2*Dirac_v1*complex<double>(kknow(1),-1.0*kknow(0))+pow(expfac2,2.0)*Dirac_v2*pow(complex<double>(kknow(0),-1.0*kknow(1)),2.0);
      out(1,0) = conj(out(0,1));
      out(1,1) = out(0,0);
  */

  // Dirac Cone approximation
  // if (full_mono_ham == 0){
    out(0,0) = Dirac_diag*2.0*kknow[dim];
    if (dim == 0){
      out(0,1) = expfac2*Dirac_v1*complex<double>(0.0,-1.0) + pow(expfac2,2.0)*Dirac_v2*complex<double>(2.0, 0.0)*complex<double>(kknow(0),-1.0*kknow(1));
    } else if (dim == 1){
      out(0,1) = expfac2*Dirac_v1*complex<double>(1.0, 0.0) + pow(expfac2,2.0)*Dirac_v2*complex<double>(0.0,-2.0)*complex<double>(kknow(0),-1.0*kknow(1));
    }
    out(1,0) = conj(out(0,1));
    out(1,1) = out(0,0);
    return out;
  // }

}


double Kp_tblg_construct::grapheneIntralayerTerm(int orbit_row, int orbit_col, array<int,2> vector){
  // convert displacement vector to Miller indices
  std::array<int, 2> hom_vec {{   2 * vector[0] +     vector[1],
                                   -1 * vector[0] +     vector[1] }};
       // redundant 3rd coordinate: -1 * vector[0] - 2 * vector[1]

    /* Shift the arrow vector by the orbital coordinates */

    if (orbit_col == 0 && orbit_row == 1)
    {
        hom_vec[0] -=  1;
        // hom_vec[2] -= -1;
    }
    else if (orbit_col == 1 && orbit_row == 0)
    {
        hom_vec[0] +=  1;
        // hom_vec[2] += -1;
    }

    /**
     * Compute the distance to the origin:
     * we use the identity r = x^2 + y^2 + z^2 = x^2 + y^2 + (x+y)^2 = 2 * (x * (x+y) + y^2)
     */
    int r = hom_vec[0] * (hom_vec[0] + hom_vec[1]) + hom_vec[1]*hom_vec[1];

    const double onsite = 0.3504;
	  //const double t_arr[9] = {0, -2.8922, 0, 0, 0, 0, 0, 0, 0}; // only nearest neighbors (debug/testing)!
    //const double t_arr[9] = {onsite, -2.8922, 0.2425, -0.2656, 0.0235, 0.0524, -0.0209, -0.0148, -0.0211};

    // coefficients from strain modeling only up to 3rd nearest-neighbor (used in generating k-dot-p terms)
    const double t_arr[9] = { 0.254*3.0, -2.822, 0.254, -0.180, 0.0, 0.0, 0.0, 0.0, 0.0};
    //const double t_arr[9] = {0.0, -2.822, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    switch (r)
    {
        case 0:
            return t_arr[0];
        case 1:
            return t_arr[1];
        case 3:
            return t_arr[2];
        case 4:
            return t_arr[3];
        case 7:
            return t_arr[4];
        case 9:
            return t_arr[5];
        case 12:
            return t_arr[6];
        case 13:
            return t_arr[7];
        case 16:
            return t_arr[8];
        default:
            return 0.0;
    }

}

// returns reciprocal of a given pair of 2d lattice vectors, a_in
vector< vector<double> > Kp_tblg_construct::getReciprocal(vector< vector<double> > a_in){

	// Here we assume a is a 2x2 matrix
	vector< vector<double> > a;
	a.resize(3);
	a[0].resize(3);
	a[0][0] = a_in[0][0];
	a[0][1] = a_in[0][1];
	a[0][2] = 0.0;
	a[1].resize(3);
	a[1][0] = a_in[1][0];
	a[1][1] = a_in[1][1];
	a[1][2] = 0.0;
	a[2].resize(3);
	a[2][0] = 0.0;
	a[2][1] = 0.0;
	a[2][2] = 1.0;

	vector< vector<double> > b;
	b.resize(2);
	for (int i = 0; i < 2; ++i){
		b[i].resize(2);
	}

	double denom1 = 0.0;
	double denom2 = 0.0;
	//double denom3 = 0.0;

	for (int j = 0; j < 3; ++j) {
		denom1 += a[0][j]*crossProd(a[1],a[2],j);
		denom2 += a[1][j]*crossProd(a[2],a[0],j);
		//denom3 += a[2][j]*crossProd(a[0],a[1],j);
	}

	for (int k = 0; k < 2; ++k){
		b[0][k] = 2*M_PI*crossProd(a[1],a[2],k)/denom1;
		b[1][k] = 2*M_PI*crossProd(a[2],a[0],k)/denom2;
		//b[2][k] = 2*M_PI*crossProd(a[0],a[1],k)/denom3;
	}

	/*
	for (int l = 0; l < 3; ++l){
		printf("b[%d] = [%lf, %lf, %lf] \n", l, b[l][0], b[l][1], b[l][2]);
	}
	*/

	return b;
}

double Kp_tblg_construct::crossProd(vector<double> x, vector<double> y, int dim){

	if (dim == 0){
		return ( x[1]*y[2] - x[2]*y[1] );
	} else if (dim == 1) {
		return ( x[2]*y[0] - x[0]*y[2] );
	} else if (dim == 2) {
		return ( x[0]*y[1] - x[1]*y[0] );
	} else {
		return 0;
	}

}


// returns complex matrix at given supercell momentum
MatrixXcd Kp_tblg_construct::getH(Vector2d k){

  // only need to construct the inter kplus/kminus terms here
  // and the in-plane (monolayer) terms.
  // The other precomputed terms are not k-dependent!

  vector< Vector2d > shift_klist_L1;
  vector< Vector2d > shift_klist_L2;

  shift_klist_L1.resize(num_hex);
  shift_klist_L2.resize(num_hex);

  for (int i = 0; i < num_hex; ++i){
    shift_klist_L1[i] = hex_all_L1[i] + k;
    shift_klist_L2[i] = hex_all_L2[i] + k;
  }

  /*
  MatrixXcd kplus_L1 = MatrixXd::Zero(2*num_hex,2*num_hex);
  MatrixXcd kminus_L1 = MatrixXd::Zero(2*num_hex,2*num_hex);
  MatrixXcd kplus_L2 = MatrixXd::Zero(2*num_hex,2*num_hex);
  MatrixXcd kminus_L2 = MatrixXd::Zero(2*num_hex,2*num_hex);

  for (int i = 0; i < num_hex; ++i){
    complex<double> L1_khere = complex<double>(shift_klist_L1[i][0],shift_klist_L1[i][1]);
    complex<double> L2_khere = complex<double>(shift_klist_L2[i][0],shift_klist_L2[i][1]);

    kplus_L1(2*i, 2*i)      = L1_khere;
    kplus_L1(2*i+1, 2*i+1)  = L1_khere;
    kminus_L1(2*i, 2*i)     = conj(L1_khere);
    kminus_L1(2*i+1, 2*i+1) = conj(L1_khere);

    kplus_L2(2*i, 2*i)      = L2_khere;
    kplus_L2(2*i+1, 2*i+1)  = L2_khere;
    kminus_L2(2*i, 2*i)     = conj(L2_khere);
    kminus_L2(2*i+1, 2*i+1) = conj(L2_khere);

  }

  // only gives H_12, need to to adjoint later for H_21

  // following line takes up majority of run time for this method (Eigen matrix-matrix is slow, single threaded)
  // so instead lets try to calculate this in the loop for kplus_L# instead
  MatrixXcd Hmat_half_inter_kk = Hmat_inter_kplus0*kplus_L1+Hmat_inter_kminus0*kminus_L1;
  */

  /*
  MatrixXcd kplus_L1 = kron(diag(shift_klist_L1(:,1)+i*shift_klist_L1(:,2)),eye(2));
  MatrixXcd kminus_L1 = kron(diag(shift_klist_L1(:,1)-i*shift_klist_L1(:,2)),eye(2));
  MatrixXcd kplus_L2 = kron(diag(shift_klist_L2(:,1)+i*shift_klist_L2(:,2)),eye(2));
  MatrixXcd kminus_L2 = kron(diag(shift_klist_L2(:,1)-i*shift_klist_L2(:,2)),eye(2));
  */

vector< complex<double> > kplus_L1(tot_dim/2);
vector< complex<double> > kminus_L1(tot_dim/2);
// unused (get them by doing adjoint later)
//vector< complex<double> > kplus_L2(tot_dim/2);
//vector< complex<double> > kminus_L2(tot_dim/2);


  for (int i = 0; i < num_hex; ++i){
    complex<double> L1_khere = complex<double>(shift_klist_L1[i][0],shift_klist_L1[i][1]);
    complex<double> L2_khere = complex<double>(shift_klist_L2[i][0],shift_klist_L2[i][1]);

    kplus_L1[2*i]     = L1_khere;
    kplus_L1[2*i+1]   = L1_khere;
    kminus_L1[2*i]    = conj(L1_khere);
    kminus_L1[2*i+1]  = conj(L1_khere);

    /*
    kplus_L2[2*i]     = L2_khere;
    kplus_L2[2*i+1]   = L2_khere;
    kminus_L2[2*i]    = conj(L2_khere);
    kminus_L2[2*i+1]  = conj(L2_khere);
    */

  }

  MatrixXcd Hmat_inter_kk = MatrixXcd::Zero(tot_dim,tot_dim);

  int rspan = tot_dim/2;
  int cspan = tot_dim/2;
  for (int r = 0; r < rspan; ++r){
    for (int c = 0; c < cspan; ++c){
      Hmat_inter_kk(all_index_L2[r],all_index_L1[c]) += Hmat_inter_kplus0(r,c)*kplus_L1[c] +  Hmat_inter_kminus0(r,c)*kminus_L1[c];
    }
  }


  // MATLAB vers for H_12:
  /*
  Hmat_inter_kk(all_index_L2(:),all_index_L1(:)) = ...
      Hmat_inter_kk(all_index_L2(:),all_index_L1(:)) +  ...
      (Hmat_inter_kplus0*kplus_L1+Hmat_inter_kminus0*kminus_L1)*1;
  */

  // H_inter = H_12 + H_21
  Hmat_inter_kk = Hmat_inter_kk + Hmat_inter_kk.adjoint().eval();


  MatrixXcd Hmat;
  // add in all precomputed terms
  Hmat = strain_fac*(Hmat_strain_L1 + Hmat_strain_L2) + inter_fac*(Hmat_inter_qdep + Hmat_inter_kk);

  // add the intralayer part
  // MATLAB vers:
  // Hmat(all_index_L1(:,indh),all_index_L1(:,indh))=Hmat(all_index_L1(:,indh),all_index_L1(:,indh))+Layer1_ham(shift_klist_L1(indh,:));
  // Hmat(all_index_L2(:,indh),all_index_L2(:,indh))=Hmat(all_index_L2(:,indh),all_index_L2(:,indh))+Layer2_ham(shift_klist_L2(indh,:));

  for (int i = 0; i < num_hex; ++i){

    Matrix2cd H_1 = layer1Ham(shift_klist_L1[i]);
    Matrix2cd H_2 = layer2Ham(shift_klist_L2[i]);

    int start_idx_L1 = all_index_L1[2*i];
    int start_idx_L2 = all_index_L2[2*i];

    for (int o1 = 0; o1 < unit_dim; ++o1) {
      for (int o2 = 0; o2 < unit_dim; ++o2){

        Hmat(start_idx_L1+o1, start_idx_L1+o2) += H_1(o1, o2);
        Hmat(start_idx_L2+o1, start_idx_L2+o2) += H_2(o1, o2);

      }
    }
  }

  return Hmat;


}

MatrixXi Kp_tblg_construct::getGToIndex(){

  return G_to_index;

}

MatrixXi Kp_tblg_construct::getIndexToG(){

  return index_to_G;

}


int Kp_tblg_construct::getSize(){
  return tot_dim;
}

MatrixXcd Kp_tblg_construct::getGradH(Vector2d k, int dim){


    // only need to construct the inter kplus/kminus terms here
    // and the in-plane (monolayer) terms.
    // The other precomputed terms are not k-dependent!

    vector< Vector2d > shift_klist_L1;
    vector< Vector2d > shift_klist_L2;

    shift_klist_L1.resize(num_hex);
    shift_klist_L2.resize(num_hex);

    for (int i = 0; i < num_hex; ++i){
      shift_klist_L1[i] = hex_all_L1[i] + k;
      shift_klist_L2[i] = hex_all_L2[i] + k;
    }

    MatrixXcd Hmat_inter_kk = MatrixXcd::Zero(tot_dim,tot_dim);

    int rspan = tot_dim/2;
    int cspan = tot_dim/2;
    for (int r = 0; r < rspan; ++r){
      for (int c = 0; c < cspan; ++c){

        //Hmat_inter_kk(all_index_L2[r],all_index_L1[c]) += Hmat_inter_kplus0(r,c)*kplus_L1[c] +  Hmat_inter_kminus0(r,c)*kminus_L1[c];

        if (dim == 0){
          Hmat_inter_kk(all_index_L2[r],all_index_L1[c]) += Hmat_inter_kplus0(r,c) + Hmat_inter_kminus0(r,c);
        } else if (dim == 1){
          Hmat_inter_kk(all_index_L2[r],all_index_L1[c]) += complex<double>(0.0,1.0)*(Hmat_inter_kplus0(r,c) - Hmat_inter_kminus0(r,c));
        }

      }
    }


    // MATLAB vers for H_12:
    /*
    Hmat_inter_kk(all_index_L2(:),all_index_L1(:)) = ...
        Hmat_inter_kk(all_index_L2(:),all_index_L1(:)) +  ...
        (Hmat_inter_kplus0*kplus_L1+Hmat_inter_kminus0*kminus_L1)*1;
    */

    // H_inter = H_12 + H_21
    Hmat_inter_kk = Hmat_inter_kk + Hmat_inter_kk.adjoint().eval();


    MatrixXcd GradHmat;
    // add in all precomputed terms
    GradHmat = inter_fac*Hmat_inter_kk;

    // add the intralayer part
    // MATLAB vers:
    // Hmat(all_index_L1(:,indh),all_index_L1(:,indh))=Hmat(all_index_L1(:,indh),all_index_L1(:,indh))+Layer1_ham(shift_klist_L1(indh,:));
    // Hmat(all_index_L2(:,indh),all_index_L2(:,indh))=Hmat(all_index_L2(:,indh),all_index_L2(:,indh))+Layer2_ham(shift_klist_L2(indh,:));

    for (int i = 0; i < num_hex; ++i){

      Matrix2cd H_1 = layer1GradHam(shift_klist_L1[i], dim);
      Matrix2cd H_2 = layer2GradHam(shift_klist_L2[i], dim);

      int start_idx_L1 = all_index_L1[2*i];
      int start_idx_L2 = all_index_L2[2*i];

      for (int o1 = 0; o1 < unit_dim; ++o1) {
        for (int o2 = 0; o2 < unit_dim; ++o2){

          GradHmat(start_idx_L1+o1, start_idx_L1+o2) += H_1(o1, o2);
          GradHmat(start_idx_L2+o1, start_idx_L2+o2) += H_2(o1, o2);

        }
      }
    }

    return GradHmat;

}
