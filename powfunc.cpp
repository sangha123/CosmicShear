#include "powfunc.h"
#include "math.h"
#include <cmath>
#include <algorithm>
#include "gsl/gsl_integration.h"
#include "gsl/gsl_sf_erf.h"
#include "gsl/gsl_errno.h"
#include "gsl/gsl_spline.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_gamma.h>
#include <omp.h>
          
#define pi 3.14
using namespace std;

Cosmology::Cosmology(vector<double> z_new,vector<double> l_new,vector<double> k_new,vector<double> zbin_new,vector<double> theta_new,double Omega_M_new,double Omega_L_new,double Omega_k_new,double h_new,double Omega_R_new,double keq_new,double ns_new,double z_CMB_new,double scale_sm_new,double H0_new,double w_new,double sigmaz_new,double zbias_new, double alpha_new, double beta_new, double z0_new,double na_new)
{
  Omega_L = Omega_L_new;
  Omega_M = Omega_M_new;
  //<< "Omega_L= " <<Omega_L << "\n";
  Omega_k = Omega_k_new;
  Omega_R = Omega_R_new;
  h = h_new;
  keq= keq_new*Omega_M*h;
  ns=ns_new;
  z_CMB=z_CMB_new;
  scale_sm=scale_sm_new; // 8 h^-1 Mpc
  H0=H0_new;
  coverh0=3000.0/h;
  
  w=w_new;
  sigmaz=sigmaz_new;
  zbias=zbias_new;
  alpha=alpha_new;
  beta=beta_new;
  z00=z0_new;					\
  na=na_new;
  // This determines the size of z,k,l
  n_z=z_new.size();
  n_k=k_new.size();
  n_l=l_new.size();
  n_theta=theta_new.size();
  num_bin=zbin_new.size();
  //cerr<<"num_bin=  "<<num_bin<<"\n";
  for(int i=0;i<n_theta;i++){
    theta.push_back(theta_new[i]);
  }

  for(int i=0;i<num_bin;i++){ zbin.push_back(zbin_new[i]);}
  for (int i=0;i<n_z;i++){ z.push_back(z_new[i]);
    a.push_back(1./(1+z_new[i]));}
  for (int i=0;i<n_l;i++){ l.push_back(l_new[i]);}
  for (int i=0;i<n_k;i++){ k.push_back(k_new[i]);}
  cerr<< "Printing parameters "<<Omega_L<<"  "<<Omega_M<<" "<<coverh0<<"  "<<ns<<" "<<w<<"\n";
}

double Cosmology::Friedmann(double a0)
{
  double val=0;
  
  val=(Omega_M/pow(a0,3) + Omega_L/pow(a0,3*(1+w)));
  return val;
}

double Cosmology::dEada(double a0)
{
  // This is derivative of Friedmann. If you update Friedmann ... please change this.

  return (-3*Omega_M/pow(a0,4));
}

double Cosmology::Comoving_int(double a0)
{
  return (coverh0/(a0*a0*pow(Friedmann(a0),0.5)));
}

double Cosmology::Comoving_dist(double z1, double z2){
  
  if (z2<z1) return(0);
  else{
    size_t neval;
    double a1,a2;
    double result, error;
    gsl_function F;
    
    F.function = call_func<Cosmology, &Cosmology::Comoving_int>;
    a1=1./(1+z1);
    a2=1./(1+z2);
    
    F.params = this;
    gsl_integration_qng(&F,a2,a1,1e-7, 1e-7,
			&result, &error,&neval);
    
    return(result);}
}

double Cosmology::RedshiftfromDist(double D)
{
  double z,eps,del_D,delz,a0;
  size_t neval;
  
  double D_new, error;
  gsl_function F;

  F.function = call_func<Cosmology, &Cosmology::Comoving_int>;
  F.params = this;
  
  z=D/coverh0;
  eps=0.01;
  del_D=1236.0;
  while (abs(del_D)>eps)
    {
      a0=1./(1.+z);
      gsl_integration_qng(&F,a0,1,1e-7, 1e-7,&D_new, &error,&neval);
      del_D=D-D_new;
      delz=0.01*del_D/D;
      z=z+delz;
    }
  return(z);
}

vector<double> Cosmology::num_dens()
{
  //cerr<<"alpha="<<alpha<<"\n";
  double B=0;
  double dz=z[10]-z[9];
  vector<double> nz;
  for(int i=0;i<n_z;i++){
    //n(z) =  n_gal /(2 z0^3 ) z^2 exp(-z/z0)
    //cerr<< alpha << "  " << beta << "  " <<z00 << "\n";
    nz.push_back(na*0.5*(pow(z[i],alpha)/pow(z00,3))* exp(-pow((z[i]/z00),beta)));
    //nz.push_back((beta/(z00*gsl_sf_gamma((1+alpha)/beta)))*pow((z[i]/z00),alpha)*exp(-pow((z[i]/z00),beta)));
    B+=nz[i]*dz;
  }
  //cerr<<"norm=  "<<B<<"\n";
  //for(int i=0;i<n_z;i++){
  //  nz[i]=nz[i]/B;}  
  
  return (nz);

}

vector<vector<double> > Cosmology::nbinz()
{
  vector<double> nz=num_dens();
  
  double xi1,xi,xi1_errf,xi_errf;
  vector<vector<double> > ni;
  vector<double> dum;
  double dz=z[10]-z[9];
  for(int i=0;i<n_z;i++)
    {
      dum.push_back(0);
    }
  for(int j=0;j<num_bin-1;j++)
    {
      ni.push_back(dum);  
    }
  for(int i=0;i<num_bin-1;i++){
    for(int j=0;j<n_z;j++){
      xi1=(zbin[i+1]-z[j]+zbias)/(pow(2,0.5)*sigmaz);
      xi=(zbin[i]-z[j]+zbias)/(pow(2,0.5)*sigmaz);
      xi1_errf=gsl_sf_erf(xi1);
      xi_errf=gsl_sf_erf(xi);
      //cout<<xi1_errf<<"  "<<xi_errf<<"  "<<xi1<<"  "<<xi<<"\n";
      ni[i][j]=0.5*nz[j]*(xi1_errf-xi_errf);
      //ni[i][j]=nz[j];
      //cout<<nz[j]<<"  "<<zbin[i]<<"\n";
    }
  }
  /*double sum;
  for (int i=0;i<num_bin-1;i++){
    sum=0;
    for(int j=0;j<n_z;j++){
      sum+=ni[i][j]*dz;
    }
    //cerr<<sum<<"\n";
    for(int j=0;j<n_z;j++){
      ni[i][j]=ni[i][j]/sum;
    }
    }*/
  return(ni);
  
}

vector<double> Cosmology::W_lens(double z0)
{
  vector<double> eta_z0z,eta_z,int_result,result; 
  double dz=z[10]-z[9];
  double sum,eta_z0,val;
  int i,j;
  vector<vector<double> > ni=nbinz();
 
  for (int i=0;i<num_bin-1;i++){                                                                                                           sum=0;                                                                                                                                for(int j=0;j<n_z;j++){                                                                                                                 sum+=ni[i][j]*dz;                                                                                                                   }                                                                                                                                                                                                                                                                       for(int j=0;j<n_z;j++){                                                                                                                 ni[i][j]=ni[i][j]/sum;                                                                                                           
    }                                                                                                                                    }
  eta_z0=Comoving_dist(0,z0);

  for(i=0;i<n_z;i++)
    {
      eta_z0z.push_back(Comoving_dist(z0,z[i]));
      eta_z.push_back(Comoving_dist(0,z[i]));
    }
  #pragma omp parallel private(i,j) reduction(+:sum)  
  for(i=0;i<num_bin-1;i++)
    {
      sum=0;
      
      for(j=0;j<n_z;j++){
	//cerr<<i<<"  "<<ni[i][j]<<"  "<<eta_z0z[j]<<"  "<<z0<<"  "<<z[j]<<"  "<<1.5*Omega_M*(1+z0)*eta_z0*(dz*ni[i][j]*eta_z0z[j]/eta_z[j])/pow(coverh0,2)<<"\n";
	val=1.5*Omega_M*(1+z0)*eta_z0*(dz*ni[i][j]*eta_z0z[j]/eta_z[j])/pow(coverh0,2);
	if (j == 0) val=0;
	sum+=val;
	//cerr<<z[i]<<"  "<<val<<" "<<sum<<"\n";
      }
      //cerr<<sum<<"\n";
      #pragma omp critical
      {result.push_back(sum);}
    }
  return(result);
}

vector<vector<vector <double> > > Cosmology::Cl_shear(vector<vector <double> > P,double kmax)
{
  double dz=z[10]-z[9];
  size_t neval=n_k;
  vector<vector<double> > W,k_new,P_new,dumdum;
  vector<vector<vector<double> > >C_ll;
  int i,j,m,p;
  vector<double> a,da,deta,dum,eta,dum_k,WW,dum_cl,k0;
  //cerr<<"num_bin= "<<num_bin<<"\n";
  for(i=0;i<n_z;i++)
    {
      dum.push_back(0);
      eta.push_back(0);
      a.push_back(1./(1.+z[i]));
      da.push_back(dz/pow((1+z[i]),2) );
      deta.push_back(da[i]*Comoving_int(a[i]));
    }
  
  for(j=0;j<n_l;j++)
    {
      k_new.push_back(dum);
      P_new.push_back(dum);
      dum_cl.push_back(0);
    }

  
  for(j=0;j<num_bin-1;j++)
    {
      W.push_back(dum);
      dumdum.push_back(dum_cl);
      
    }

  for(j=0;j<num_bin-1;j++)
    {
      C_ll.push_back(dumdum);
    }

  //cerr<<"Size Cl  "<<C_ll.size()<<"  "<<C_ll[0].size()<<"  "<<C_ll[0][0].size()<<"\n";

  for(i=0;i<n_k;i++){
    dum_k.push_back(0);
    k0.push_back(0);
  }
  //ofstream pspl,dumk;
  //pspl.open("pspl.dat");
  //dumk.open("dumk.dat");
#pragma omp parallel for private(i,j,m,WW) firstprivate(k0,dum_k)       
  for(i=0;i<n_z;i++)
    {
      gsl_interp_accel *acc = gsl_interp_accel_alloc ();
      gsl_spline *spline = gsl_spline_alloc (gsl_interp_cspline, neval);
      eta[i]=Comoving_dist(0,z[i]);
      
      WW = W_lens(z[i]);

      for(m=0;m<n_k;m++)
	{
	  k0[m]=k[m];
	  dum_k[m]=(P[m][i]);
	  //dumk<<i<<"  "<<m<<"  "<<k0[m]<<"  "<<dum_k[m]<<"\n";
	}
      
      gsl_spline_init(spline, &k0[0], &dum_k[0], neval);

	for(j=0;j<n_l;j++){
	  k_new[j][i]=l[j]/(eta[i]);
	  if (k_new[j][i]<*min_element(k.begin(),k.end())) {P_new[j][i]=(P[0][i]/pow(k[0],ns))*pow(k_new[j][i],ns);
	  }
	  
	  else if (k_new[j][i]>= kmax) {
	    P_new[j][i]=0.0;
	  }
	  else{
	    P_new[j][i]=gsl_spline_eval(spline, k_new[j][i], acc);
	  }
	  //pspl<<j<<"  "<<i<<"  "<<z[i]<<"  "<<k_new[j][i]<<"  "<<P_new[j][i]<<"\n";
	}
      for(j=0;j<num_bin-1;j++){
	W[j][i]=WW[j];
	  }
      
  gsl_spline_free (spline);
  gsl_interp_accel_free (acc);
    }
  //pspl.close();
  //dumk.close();
  
  double sum;
  #pragma omp parallel for private(i,j,m,p) reduction(+:sum)
    
  for(p=0;p<num_bin-1;p++){
    for(i=0;i<num_bin-1;i++){
      for(j=0;j<n_l;j++) {
	sum=0;
	for(m=0;m<n_z;m++){
	  //cerr<<W[p][m]*W[i][m]*P_new[j][m]*deta[m]/pow(eta[m],2)<<"\n";
	  if (m>0) sum+=W[p][m]*W[i][m]*P_new[j][m]*deta[m]/pow(eta[m],2);
	}
	C_ll[p][i][j]=sum;
	//cerr<<C_ll[p][i][j]<<"\n";
      }}}
  
  return C_ll;

}
/*
vector< vector<double> > Cosmology::Zi(vector< vector<double> > P)
// This is the real space correlation function.

{
  vector< vector<double> > zi;
  vector<double> dum_zi;
  double sum,j0;
  double kmax=20;
  double dl=l[1]-l[0];
  vector< vector< vector<double> > > Cl=Cl_shear(P,kmax);
  for (int i=0;i<n_theta;i++){
    dum_zi.push_back(0.0);
  }

  for (int i=0;i<num_bin-1;i++){
    zi.push_back(dum_zi);
  }

  //cerr<<n_theta<<"  "<<num_bin<<"  zi=  "<<zi.size()<<"  "<<zi[0].size()<<"\n";
  for (int i=0;i<n_theta;i++)
    {
      for (int m=0;m<num_bin-1;m++){
	sum=0;
	for(int j=0;j<n_l;j++){
	  j0=gsl_sf_bessel_J0(l[j]*theta[i]);
	  //cerr<<"TEST  "<<l[j]<<"  "<<theta[i]<<"  "<<"  "<<j0<<"\n";
	  sum+=(dl*l[j]/2*pi)*Cl[m][j]*j0;
	    }
      
	zi[m][i]=sum;
	  }
    }
  return(zi);

}

*/
/*
vector< vector<double> > Cosmology::Map(vector< vector<double> > P)
//function calculates aperturemass
{



}
*/

//func, jac are functions required by gsl to solve the ODE for the growth function. They are private.

int Cosmology::func(double a0, const double y[], double f[])
{
  f[0] = y[1];
  f[1] = -(3./a0+0.5*dEada(a0)/Friedmann(a0))*y[1]+(3.*Omega_M/(2*pow(a0,5)*Friedmann(a0)))*y[0];
  
  return GSL_SUCCESS;
}

int Cosmology::jac (double a0, const double y[], double *dfdy, 
	 double dfdt[])
{
  
  gsl_matrix_view dfdy_mat 
    = gsl_matrix_view_array (dfdy, 2, 2);
       gsl_matrix * m = &dfdy_mat.matrix; 
       gsl_matrix_set (m, 0, 0, 0.0);
       gsl_matrix_set (m, 0, 1, 1.0);
       gsl_matrix_set (m, 1, 0, 3.0*Omega_M/(2*pow(a0,5)*Friedmann(a0)));
       gsl_matrix_set (m, 1, 1, -(3./a0+0.5*dEada(a0)/Friedmann(a0)));
       dfdt[0] = 0.0;
       dfdt[1] = 0.0;
       return GSL_SUCCESS;
}

vector<double> Cosmology::Growth()
{
  vector<double> G,G_new;
  vector<double> z_new,a_new;
  for(int i=0;i<n_z;i++){
    z_new.push_back(z[i]);}
  
  int zstart=floor(*max_element(z.begin(),z.end())+0.4);
  //cerr<<"  "<<zstart<<"\n";
  for(int i=zstart;i<(z_CMB);i++)
    {
      z_new.push_back(i);
    }
  

  gsl_odeiv2_system sys = {call_func2<Cosmology, &Cosmology::func>, call_func3<Cosmology, &Cosmology::jac>, 2,this};
    gsl_odeiv2_driver * d = 
      gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk8pd,
				     1e-6, 1e-6, 0.0);
    
    int nz_new=z_new.size();
    for(int i=0;i<nz_new;i++)
      {
	a_new.push_back(1./(1+z_new[i]));
	G_new.push_back(0.0);
      }
    double t = a_new[nz_new-1];
    double y[2] = {a_new[nz_new-1],1.0};  
    for (int i = 1; i < nz_new; i++)
      {
	int status = gsl_odeiv2_driver_apply (d, &t, a_new[nz_new-i], y);
     	if (status != GSL_SUCCESS)
	  {
	    printf ("error, return value=%d\n", status);
	    break;
	  }
	G_new[nz_new-i]=y[0];
	  //printf ("%.5e %.5e %.5e\n", t, y[0], y[1]);
      }
    
    for (int i=0;i<n_z;i++){
      G.push_back(G_new[n_z-i]);
    }
    
    for(int i=0;i<n_z;i++){
      G[i]=G[i]/G[0];
    }

    gsl_odeiv2_driver_free (d);
    //cerr<<"TEST2"<<"\n";
    return(G);

}

/*
vector< vector <vector<double> > > Cosmology::Covariance(double fsky, int ngal, double sigma_e,vector< vector<double> > P)
//fsky is the fraction of sky covered by the survey
//ngal : number of galaxies per square arc-minute.
{

  vector< vector <vector <double> > > Cov;
  vector< vector <double> >dum_cc;

  vector<double> dum_cov;
  double dl=l[1]-l[0];
  double j01,j02,sum;
  double kmax=20.0;
  vector< vector<double> > Cl=Cl_shear(P,kmax);

  for(int i=0;i<n_theta;i++){
      dum_cov.push_back(0);
    }

  for(int i=0;i<n_theta;i++){
    dum_cc.push_back(dum_cov);
  }

  for(int i=0;i<num_bin-1;i++){
    Cov.push_back(dum_cc);
  }

  for(int i=0;i<num_bin-1;i++){
    for(int j=0;j<n_theta;j++){
      for(int m=0;m<n_theta;m++){
	sum=0;
	for (int p=0;p<n_l;p++){
	  j01=gsl_sf_bessel_J0(l[p]*theta[j]);
	  j02=gsl_sf_bessel_J0(l[p]*theta[m]);
	  //cerr<<j01<<"  "<<j02<<"\n";
	  sum+=dl*l[p]*j01*j02*pow((Cl[i][p]+0.5*pow(sigma_e,2)/ngal),2);
	}
	Cov[i][j][m]=sum;
      }
    }
  }

  return(Cov);
}

*/
