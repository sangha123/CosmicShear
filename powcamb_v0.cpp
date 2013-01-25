//This code takes the power spectrum from the camb and returns the shear power spectrum.
// call [paramfile]
// 
// 
// om h^2 (This includes dark matter and baryonic matter)
// omb h^2
// ns
// sigma_8
// w
// h
// name of the file containing the matter power spectrum (Example: /users/astro/sdeb/camb/NLResults/concat.dat)
// name of the outputfile for Cl's
// name of the file containing redshifts
// alpha (Photoz parameters Ma Hu & Huterer)
// beta
// na
// sigmaz
// zbias
// kmax
// id  (id=0 represents camb, id=1 represents Emulator. This is important since emulator does not require multiplication by h.)
// zbin (This is an array of the form 0.6 1.2 2.0 3.0)
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include <algorithm>
#include "powfunc.h"
#include <string>
#include <sstream>
#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <ctime>
#include <sstream>
using namespace std;

int main(int argc,char *argv[]){

  time_t start;
  double duration;
  start = time(0);

  int n=10000,n_theta=100;
  double alpha,beta,kmax;
  double zpivot,na,sigmaz,zbias; // parameters for the number density of 
  vector<double> z;
  vector<double>k,k0, P0,eta;
  vector<double> l;
  vector<double> zbin;
  vector<double> theta;
  double val,om,omb,ns,w,sigma_8,number,h;
  int num_bin;
  int id;
  /////////////////////////////////////////////////
  // These are some default values.
  //
  double Omega_M_new=0.1375;
  double ns_new=0.961;
  double scale_sm_new=0.7558;
  double w_new=-0.998;
  double Omega_L_new=1-Omega_M_new;
  double z_CMB_new=1100;
  double Omega_k_new=0;
  double h_new=0.69;
  double Omega_R_new=0;
  double keq_new=0;
  double H0_new=100;
  double z0_new=0.373;
  //
  //////////////////////////////////////////////////

  if (argc < 2) {
    cerr << "Usage: " << argv[0] << " [input-file-name3: paramfile] "<<"\n";
    return 1;
  }

  fstream configfile;
  configfile.open(argv[1],ios::in);
  if (! configfile.is_open()) {
    cerr<<" File does not exist"<<"\n";
    return(0);
  }

  char zfile[500],outfile[500],pkfile_camb[500],Clfile_camb[500],Clfile_emu[500],zbinfilecall[500];
  
  string line;

  configfile>>om;
  configfile>>omb;
  configfile>>ns;
  configfile>>sigma_8;
  configfile>>w;
  
  configfile>>h;
  configfile>>pkfile_camb;
  configfile>>Clfile_camb;
  configfile>>zfile;

  configfile>>alpha;

  //getline(configfile,line);
  //cerr<<"Line= "<<line<<"\n";
  configfile>>beta;
  //getline(configfile,line);
  configfile>>na;
  //getline(configfile,line);
  configfile>>sigmaz;
  //getline(configfile,line);
  configfile>>zbias;
  configfile>>kmax;
  configfile>>id;
  cerr<<"Type=  "<<id<<"\n";
  double d00;
  if (id==0) kmax=kmax/h;
  while(getline(configfile,line,',')){
    stringstream ss(line);
    ss>>d00;
    zbin.push_back(d00);
    
  }
  cerr<<"Size zbin= "<<zbin.size()<<"\n";


  configfile.close();
  cerr<<"#/str  "<<na<<"\n";

  cerr<<"zfile ="<<zfile<<"\n";
  cerr<<"printing input params  "<<Clfile_camb<<" "<<om<<"\n";
  cerr<<"beta=  "<<beta<<"  "<<na<<"  "<<sigmaz<<"  "<<zbias<<"\n" ;
  
  Omega_M_new=(om+omb)/(h*h);
  ns_new=ns;
  scale_sm_new=sigma_8;
  w_new=w;
  Omega_L_new=1-Omega_M_new;
  z_CMB_new=1100;
  Omega_k_new=0;
  h_new=h;
  Omega_R_new=0;
  keq_new=0;
  H0_new=100;
  
  int c=0;
  fstream  redshift;
  redshift.open(zfile,ios::in);
  double d1,d2;
  while(!(redshift.eof() || redshift.peek()== EOF )){
    redshift >>d1;
    z.push_back(d1);
    getline(redshift,line);
  }
  redshift.close();

  int n_z=z.size();
  fstream pkfile;
  
  pkfile.open(pkfile_camb,ios::in);

  while(!(pkfile.eof() || pkfile.peek()== EOF )){
    pkfile >>d1>>d2;
    k0.push_back(d1);
    P0.push_back(d2);
    getline(pkfile,line);
  }
  pkfile.close();

  num_bin=zbin.size();
  cerr<<"numbin  "<<num_bin<<"\n";

  cerr<<"k0_size= "<<k0.size()<<"\n";
  cerr<<"z_size = "<<n_z<<"\n";
  int n_k=(k0.size()+1)/n_z;
  cerr<<"size_k1  "<<n_k<<"\n";
  for(int i=0;i<n_k;i++){
    if (id==0) k.push_back(k0[i]*h);
    if (id==1) k.push_back(k0[i]);
  }
  vector< vector<double> > Pk;
 
  cerr<<"size_k  "<<k.size()<<"\n";
  int n_pk=P0.size();
  cerr<<"P_size= "<<n_pk<<"\n";
  cerr<<"zbin_size= "<<zbin.size()<<"\n";
  /*for(int i=0;i<num_bin;i++)
    {
      val=i;
      zbin.push_back(val*1.2+0.001);
      cerr<<zbin[i]<<"\n";
    }
  */
  cerr<<zbin[0]<<"\n";
  for(int i=0;i<n;i++){
    val=i+2;
    l.push_back(val);
  }

  for(int i=0;i<n_theta;i++){
      val=(i+0.5)*60./206265.;
      theta.push_back(val);
  }
  cerr<<"now cosmology "<<"\n";
  
  Cosmology C(z,l,k,zbin,theta,Omega_M_new=Omega_M_new,Omega_L_new=Omega_L_new,Omega_k_new=Omega_k_new,h_new=h_new, Omega_R_new=Omega_R_new,keq_new=keq_new,ns_new=ns_new,z_CMB_new=z_CMB_new,scale_sm_new=scale_sm_new,H0_new=H0_new,w_new=w_new,sigmaz=sigmaz,zbias=zbias, alpha=alpha,beta=beta,z0_new=z0_new,na=na);
  
  vector<double> W;
  for(int i=0;i<n_z;i++){
    W.push_back(C.W_lens(z[i])[0]);
    eta.push_back(C.Comoving_dist(0,z[i]));
  }
  vector<double> nz=C.num_dens();
  vector<vector<double> > nbinz=C.nbinz();
  cerr<<"now growth"<<"\n";
  vector<double> G=C.Growth();
  ofstream gw;
  gw.open("nz.dat");
  for(int i=0;i<n_z;i++){
    gw << z[i] <<"  "<<nz[i]<<"  "<<nbinz[0][i]<<"  "<<W[i]<<"  "<<eta[i]<<"\n";
  }
  gw.close();
  

  //Calculating nbar for everybin
  vector<double> nbar;
  double sum,dz;
  dz=z[10]-z[9];
  //vector<vector<double> > nbinz=C.nbinz(); 
  for(int i=0;i<num_bin-1;i++){
    sum=0;
    for(int j=0;j<n_z;j++){
      sum+= nbinz[i][j]*dz;
    }
    nbar.push_back(sum);
  }
  
  ofstream nbarfile;
  nbarfile.open("nbar.dat");
  for(int i=0;i<num_bin-1;i++){
    nbarfile <<nbar[i]<<"\n";
  } 
  nbarfile.close();
  vector<double> dum_pk;

  for(int i=0;i<n_z;i++){
    dum_pk.push_back(0);
  }
  for (int j=0;j<n_k;j++){
    Pk.push_back(dum_pk);
  }
 
  for(int j=0;j<n_z;j++){
    for(int i=0;i<n_k;i++){
      if (id==0) Pk[i][j]=P0[i+j*n_k]/(pow(h_new,3));
      if (id==1) Pk[i][j]=P0[i+j*n_k];
    }
  }

  ofstream cldata, zidata, covdata;
  cldata.open(Clfile_camb);
  cerr<<"Now Pk  "<<"\n";
  
  vector<vector<vector<double> > > Cl=C.Cl_shear(Pk,kmax);
  cerr <<"Size_Cl=  " << Cl.size() <<"  "<< Cl[1].size() <<"  "<< Cl[0].size() <<"\n";
  for (int i=0;i<Cl.size();i++){                                                                                                                                           for (int p=0;p<Cl[1].size();p++){
      for (int j=0;j<Cl[0][0].size();j++){                                                                                    
	cldata<< i <<"  "<< p <<"  "<< j <<"  "<<l[j]<<"  "<<Cl[i][p][j]<<"\n";                                                                                                            
      }                                                                                                                                                                }
  }                                                                                                                            
  cldata.close();   
  
  /*
  vector< vector<double> > Cl_kbin1,Cl_kbin2;
  vector<double> Cl0,Cl1;
  for(int i=0;i<n;i++) {
    Cl0.push_back(0);
    Cl1.push_back(0);
  }
  int nkbin=10;
  for(int i=0;i<n;i++) {
    Cl0[i]=Cl[0][i];
    Cl1[i]=Cl[1][i];
  }

  Cl_kbin1.push_back(Cl0);
  Cl_kbin2.push_back(Cl1);
  
  for (int i=0;i<nkbin;i++){
    kmax=2.0*i+0.1;
    vector<vector<double> > Clk=C.Cl_shear(Pk,kmax);
    for (int j=0;j<Clk[0].size();j++){ 
      Cl0[j]=(Clk[0][j]);
      Cl1[j]=Clk[1][j];
    }
    Cl_kbin1.push_back(Cl0);
    Cl_kbin2.push_back(Cl1);
  }

  for (int i=0;i<Cl_kbin1.size();i++){                                                                                                                  
    for (int j=0;j<Cl_kbin1[0].size();j++){                                                                                                             
      cldata<< i <<"  "<< j <<"  "<<l[j]<<"  "<<Cl_kbin1[i][j]<<"  "<<Cl_kbin2[i][j]<<"\n";
    }                                                                                                                                                 
  }                                                                                                                                                     
  cldata.close(); 
  */
  duration = ( time(0) - start );
  cerr<<"Total Time=  "<<duration<<"\n";

}
    
