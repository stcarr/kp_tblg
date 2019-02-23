/*
 * File:   main.cpp
 * Author: Stephen Carr
 *
 * Created on January 18, 2018, 4:28 PM
 */

#include "kp_tblg_construct.h"
#include <fstream>
#include <math.h>
#include <iostream>

using namespace std;

int main(int argc, char** argv) {

  double theta = 1.0;//1.00;

  // empty constructor
  Kp_tblg_construct kp_tool;


  // loads relaxation data (for all thetas, ~5 degrees to ~0.5 degrees)
  string filename = "full_relax_kp_01-06-2019.dat";
  kp_tool.loadFiles(filename);
  printf("Done with file loads \n");

  // precompute most parts of the k-p Hamiltonian
  kp_tool.setTwist(theta);
  kp_tool.prepare();
  printf("Done with H setup \n");

  // Save gamma point H to file
  /*

  Vector2d k_gamma(0.0,0.0);
  MatrixXcd H_gamma = kp_tool.getH(k_gamma);

  ofstream file;
  file.open("H_gamma.txt");
  int ncols = H_gamma.cols();
  int nrows = H_gamma.rows();

  if (file.is_open())
  {
    for (int r = 0; r < nrows; ++r){
      for (int c = 0; c < ncols; ++c){

        if ( H_gamma(r,c).imag() >= 0)
          file << H_gamma(r,c).real() << "+" << H_gamma(r,c).imag() << "i" << "    ";
        else
        file << H_gamma(r,c).real() << H_gamma(r,c).imag() << "i" << "    ";

      }
      file << '\n';

    }
  }
  file.close();

  */

  //Vector2d zero_k(0.0,0.0);
  //Matrix2cd test_mat = kp_tool.layer1Ham(zero_k);
  //cout << test_mat;
  // Computes Eigenvalues:
  // /*

  double rot_theta = theta*M_PI/180;
  double lattice_a = 1.42*sqrt(3);

  // some important length scales for the Brill. Zone
  double KD = 4*M_PI/3/lattice_a;
  double KTH = 2*KD*sin(rot_theta/2.0);
  double HEX_BLEN=KTH*sqrt(3);

  Vector2d hex_b1(HEX_BLEN*sqrt(3.0)/2.0, -HEX_BLEN*1.0/2.0);
  Vector2d hex_b2(HEX_BLEN*0.0, HEX_BLEN*1.0);
  Vector2d hex_shift = (-hex_b1+hex_b2)/3.0;


  // Default K point sampling, goes along:
  // K G M K

   // K
   Vector2d kk1a = (hex_b1-hex_b2)/3.0;
   Vector2d kk1b = (hex_b1+2.0*hex_b2)/3.0;
   Vector2d kk1c = (-2.0*hex_b1-hex_b2)/3.0;

   // Ks'
   Vector2d kk2a = -(hex_b1-hex_b2)/3.0;
   Vector2d kk2b = -(hex_b1+2.0*hex_b2)/3.0;
   Vector2d kk2c = -(-2.0*hex_b1-hex_b2)/3.0;

   Vector2d kk3 = 0.0*hex_b1;
   Vector2d kk4 = hex_b1/2.0;

   int num_sections = 3;
   int nk = 20;
   vector < Vector2d > k_scan;

   vector <Vector2d> k_endpts;
   k_endpts.resize(2*(num_sections + 1));

   k_endpts[0] = kk1a;
   k_endpts[1] = kk3;
   k_endpts[2] = kk4;
   k_endpts[3] = kk1a;

   k_endpts[4] = -kk1a;
   k_endpts[5] = -kk3;
   k_endpts[6] = -kk4;
   k_endpts[7] = -kk1a;

   Matrix2cd test_mat = kp_tool.layer1Ham(0.5*kk1a+hex_shift);
   cout << test_mat << "\n";
   Vector2d kk1a_mirror;
   kk1a_mirror(0) = cos(-M_PI/3.0)*kk1a(0) - sin(-M_PI/3.0)*kk1a(1);
   kk1a_mirror(1) = sin(-M_PI/3.0)*kk1a(0) + cos(-M_PI/3.0)*kk1a(1);

   printf("kk1a   = [%lf, %lf] \n",kk1a(0),kk1a(1));
   printf("kk1a_m = [%lf, %lf] \n",kk1a_mirror(0),kk1a_mirror(1));

   test_mat = kp_tool.layer2Ham(0.5*kk1a_mirror-hex_shift);
   cout << test_mat << "\n";

   for (int valley = 0; valley < 2; ++valley){
     for (int n_sec = 0; n_sec < num_sections; ++n_sec){
       Vector2d vStart = k_endpts[valley*4 + n_sec];
       Vector2d vEnd = k_endpts[valley*4 + n_sec +1];

       for (int k = 0; k < nk; ++k){

         double dk = ((double) k)/((double) nk);
         Vector2d k_here = (1.0-dk)*vStart + dk*vEnd;
         k_scan.push_back(k_here);

       }
     }
   }

   int tot_k = k_scan.size();

   int tot_dim = kp_tool.getSize();

   SelfAdjointEigenSolver<MatrixXcd> es(tot_dim);
   ofstream file;
   file.open("eigenvals.txt");

   for (int k = 0; k < tot_k; ++k){

     MatrixXcd H = kp_tool.getH(k_scan[k]);
     printf("Done with H construct %d/%d \n",k+1,tot_k);

     es.compute(H,EigenvaluesOnly);
     file << es.eigenvalues().transpose() << endl;

     printf("Done with k-point %d/%d \n",k+1,tot_k);

   }
   file.close();
  // End of Eigenvalue calc

  return 0;
}
