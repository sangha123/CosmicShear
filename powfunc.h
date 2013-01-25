#include <vector>
#include <iostream>
#include <string>
#include <fstream>
using namespace std;

class Cosmology{
 public:
  Cosmology(vector<double> z, vector<double> k, vector<double> l,vector<double> zbin,vector<double> theta,double Omega_M=0.272,double Omega_L=0.728,double Omega_k=1.e-10,double h=0.703,double Omega_R=8.4e-5,double keq=0.073,double ns=1,double z_CMB=1100.0,double scale_sm=8.07,double H0=100.,double w=-1,double sigmaz=0.001,double zbias=0.0, double alpha=1.197, double beta=1.193, double z00=0.555,double na=40.0);
  double Friedmann(double a0);
  double dEada(double a0);
  double Comoving_int(double a0);
  double Comoving_dist(double z1, double z2);
  double RedshiftfromDist(double D);
 
  vector<double> num_dens();
  vector<vector <double> > nbinz();

  //na: number of galaxies per steradian
  vector<double> W_lens(double z0);
  vector<vector< vector<double> > >Cl_shear(vector<vector<double> > P,double kmax);
  vector<double> Growth();
  //vector< vector<double> > Zi(vector< vector<double> > P);
  //vector< vector <vector<double> > > Covariance(double fsky, int ngal, double sigma_e,vector< vector<double> > P);
 protected:
  double Omega_M, Omega_L,Omega_k,h,Omega_R,keq,ns,z_CMB,scale_sm, H0,w,coverh0;
  int n_z,n_l,n_k,num_bin,n_theta;
  double sigmaz,zbias,na,z00,alpha,beta;
  int func(double a0, const double y[], double f[]);
  int jac (double a0, const double y[], double *dfdy,double dfdt[]);
  vector<double> z;
  vector<double> l;
  vector<double> a;
  vector<double> k;
  vector<double> zbin;
  vector<double> theta;
};

template <class C, double (C::*func)(double)>
double call_func(double x, void *cxt) {
  C *c = static_cast<C*>(cxt);
  return (c->*func)(x);
}; 

template <class C, int (C::*func)(double, const double[], double[])>
int call_func2(double x1, const double y[], double f[],void *cxt) {
  C *c = static_cast<C*>(cxt);
  return (c->*func)(x1, y, f);
};

template <class C, int (C::*func)(double, const double[], double*, double[])>
  int call_func3(double x2, const double y[], double *dfdy,double dfdt[],void *cxt) {
  C *c = static_cast<C*>(cxt);
  return (c->*func)(x2, y, dfdy, dfdt);
};

