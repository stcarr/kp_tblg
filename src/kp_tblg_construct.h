/*
 * File:   kp_tblg_construct.h
 * Author: Stephen Carr
 *
 * Created on January 13, 2019, 3:41 PM
 */

#ifndef KP_TBLG_CONSTRUCT_H
#define KP_TBLG_CONSTRUCT_H

// Dense matrices only, can't eigensolve sparse ones!
#include <Eigen/Dense>

#include <vector>
#include <complex>
#include <array>

using namespace std;
using namespace Eigen;

class Kp_tblg_construct {

  private:

    double theta; // twisting angle
    double inter_fac; // scales interlayer couplings
    double strain_fac; // scales strain (intralayer) couplings
    int full_mono_ham; // 0: Use dirac-cone approx, 1: Use full monolayer Hamiltonians


    // data storage sizes
    int ntheta, n_inter_couplings, n_intra_couplings,n_inter_shells,n_intra_shells;

    vector<double> thetas_deg;

		vector< vector < vector < vector< complex<double> > > > >  All_Eff_inter;
    vector< vector < vector < vector< complex<double> > > > >  All_Eff_inter_kplus;
    vector< vector < vector < vector< complex<double> > > > >  All_Eff_inter_kminus;
    vector< vector < vector < vector< complex<double> > > > >  All_Eff_intra_bot;
    vector< vector < vector < vector< complex<double> > > > >  All_Eff_intra_top;

    vector< vector <vector<int> > >     All_Eff_inter_shell_indices;
    vector< vector <vector<int> > >     All_Eff_intra_shell_indices;

    // terms for Monolayer Ham construction
    complex<double> expfac1;
    complex<double> expfac2;
    double Dirac_v1;
    double Dirac_v2; // no quadratic for now
    double Dirac_diag;

    // variables for H construction

    int num_hex;
    int unit_dim;
    int tot_dim;
    vector<Vector2d> hex_all_L1;
    vector<Vector2d> hex_all_L2;

    vector<int> all_index_L1;
    vector<int> all_index_L2;
    vector<int> all_index;

    MatrixXcd Hmat_inter_qdep;
    MatrixXcd Hmat_inter_kplus0;
    MatrixXcd Hmat_inter_kminus0;

    MatrixXcd Hmat_strain_L1;
    MatrixXcd Hmat_strain_L2;

    MatrixXi index_to_G;
    MatrixXi G_to_index;


	public:

    // Empty constructor
    Kp_tblg_construct();

    // Copy constructor
    Kp_tblg_construct(const Kp_tblg_construct& orig);

    // Destructor
    ~Kp_tblg_construct();

    // constructor with angle supplied
    Kp_tblg_construct(double theta_in);

    // define twisting angle, can be changed before call to prepare
    void setTwist(double theta_in);

    // set the strength of Interlayer and Strain terms in the Hamiltonian
    void setInterFac(double interfac_in);
    void setStrainFac(double strainfac_in);

    // 0: use Dirac cone approximation, 1: Use full monolayer Hamiltonian
    void setFullMonoHam(int fullmonoham_in);

    // load kp data from file
    void loadFiles(string filename);

    // for custom 1D interpolation
    vector< vector< vector< complex<double> > > > interpKP(double interp_scale, vector< vector< vector< complex<double> > > > high, vector< vector< vector< complex<double> > > > low);

    // sets up the sparsity pattern for the Hamilotnian
    void prepare();

    // returns high symmetry pts in supercell BZ (for the currently set theta)
    Vector2d getK();
    Vector2d getM();
    Vector2d getGamma();


    // simple functions to return 2x2 monolayer Hamlitonian at a given k
    Matrix2cd layer1Ham(Vector2d kknow);
    Matrix2cd layer2Ham(Vector2d kknow);

    Matrix2cd layer1GradHam(Vector2d kknow, int dim);
    Matrix2cd layer2GradHam(Vector2d kknow, int dim);

    // returns monolayer couplings for accurate onsite Hamiltonian construction
    double grapheneIntralayerTerm(int orbit_row, int orbit_col, array<int,2> vector);

    // for computing reciprocal lattice
    vector< vector<double> > getReciprocal(vector< vector<double> > a_in);
    double crossProd(vector<double> x, vector<double> y, int dim);

    // returns complex matrix at given supercell momentum
    MatrixXcd getH(Vector2d k);

    MatrixXi getGToIndex();
    MatrixXi getIndexToG();


    // returns the k-basis for the most recent run
    MatrixXcd getKbasis();



    // returns the size of matrix returned by getH with current settings
    int getSize();

    // returns k-derivative of H at given supercell momentum
    MatrixXcd getGradH(Vector2d k, int dim);


};

#endif /* KP_TBLG_CONSTRUCT_H */
