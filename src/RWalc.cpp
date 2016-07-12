/* Rwalc.cpp

This code is heavily based on and is a simmplification of the TMB code
from the argosTrack package (https://github.com/calbertsen/argosTrack).

*/

#include <TMB.hpp>

using namespace density;


/** \brief Multivariate t distribution with user supplied scale matrix

Class to evaluate the negative log density of a multivariate t distributed variable with general scale matrix Sigma and location vector 0 and df degrees of freedom.
*/
template <class Type>
class MVT_t: public MVNORM_t<Type>
{
  Type df;
  bool useNorm;

public:
  MVT_t()
    : MVNORM_t<Type>()
  {
    useNorm = false;
  };
  MVT_t(Type df_)
    : MVNORM_t<Type>()
  {
    df = df_;
    useNorm = false;
  }
  MVT_t(matrix<Type> Sigma_, Type df_)
    : MVNORM_t<Type>(Sigma_)
  {
    df = df_;
    useNorm = false;
  }
  MVT_t(matrix<Type> Sigma_, Type df_, bool useNorm_)
    : MVNORM_t<Type>(Sigma_)
  {
    df = df_;
    useNorm = useNorm_;
  }

  void setdf(Type df_){
    df = df_;
  }

  /** \brief Evaluate the negative log density */
  Type operator()(vector<Type> x){
    Type p = x.size();
    //Lange et al. 1989 http://www.jstor.org/stable/2290063
    Type tdens = -lgamma(Type(0.5)*(df+p))+lgamma(Type(0.5)*df)+p*Type(0.5)*log(df)+p*lgamma(Type(0.5))-Type(0.5)*this->logdetQ + Type(0.5)*(df+p)*log(Type(1.0)+this->Quadform(x)/df);
    Type ndens = -Type(.5)*this->logdetQ + Type(.5)*this->Quadform(x) + p*Type(log(sqrt(2.0*M_PI)));

    if(useNorm) return ndens; else return tdens;
  }
};




template<class Type>
Type objective_function<Type>::operator() ()
{


  DATA_MATRIX(y);               // (lon, lat) observations
  DATA_MATRIX(w);               // Error weights
  DATA_VECTOR(dt);              // Time steps
  DATA_IVECTOR(obs);            // Indices of the observed states
  DATA_IVECTOR(seg);            // Indices of the observed states
  DATA_SCALAR(tdf);             // Degrees of freedom of the t distribution
  DATA_VECTOR(bshrink);         // Shrinkage penalty for beta
  PARAMETER_VECTOR(logBeta);
  PARAMETER_VECTOR(logSigma);
  PARAMETER_VECTOR(logTau);
  PARAMETER_MATRIX(mu);
  PARAMETER_MATRIX(nu);


  vector<Type> beta = exp(logBeta);
  vector<Type> sigma2 = exp(Type(2.0)*logSigma);

  // Backtransform parameters
  vector<Type> sigma = exp(logSigma);
  vector<Type> tau = exp(logTau);

  // Response distribution
  matrix<Type> H(2,2);
  H.setZero();
  H(0,0) = exp(2.0*logTau(0));
  H(1,1) = exp(2.0*logTau(1));

  MVT_t<Type> nll_dist = MVT_t<Type>(H,Type(tdf),tdf <= 0);
  //MVNORM_t<Type> nll_dist = MVNORM<Type>(H);


  Type nll = 0.0;
  vector<Type> eta(4);
  matrix<Type> Q(4,4);
  Q.setZero();

  for(int i=0; i<dt.size(); ++i) {
    if(seg(i+1)==seg(i)) {
      eta(0) = mu(i+1,0)-(mu(i,0)+nu(i,0)*(1.0-exp(-beta(0)*dt(i)))/beta(0));
      eta(1) = nu(i+1,0)-nu(i,0)*exp(-beta(0)*dt(i));

      eta(2) = mu(i+1,1)-(mu(i,1)+nu(i,1)*(1.0-exp(-beta(1)*dt(i)))/beta(1));
      eta(3) = nu(i+1,1)-nu(i,1)*exp(-beta(1)*dt(i));

      Q(0,0) = sigma2(0)*(dt(i)-2.0*(1.0-exp(-beta(0)*dt(i)))/beta(0)+(1.0-exp(-2.0*beta(0)*dt(i)))/(2.0*beta(0)))/pow(beta(0),2.0);
      Q(1,1) = sigma2(0)*(1.0-exp(-2.0*beta(0)*dt(i)))/(2.0*beta(0));
      Q(1,0) = sigma2(0)*(1.0-2.0*exp(-beta(0)*dt(i))+exp(-2.0*beta(0)*dt(i)))/(2.0*pow(beta(0),2.0));
      Q(0,1) = Q(1,0);

      Q(2,2) = sigma2(1)*(dt(i)-2.0*(1.0-exp(-beta(1)*dt(i)))/beta(1)+(1.0-exp(-2.0*beta(1)*dt(i)))/(2.0*beta(1)))/pow(beta(1),2.0);
      Q(3,3) = sigma2(1)*(1.0-exp(-2.0*beta(1)*dt(i)))/(2.0*beta(1));
      Q(2,3) = sigma2(1)*(1.0-2.0*exp(-beta(1)*dt(i))+exp(-2.0*beta(1)*dt(i)))/(2.0*pow(beta(1),2.0));
      Q(3,2) = Q(2,3);

      nll += MVNORM<Type>(Q)(eta);
    }
  }

  vector<Type> epsilon(2);
  for(int i=0; i<obs.size(); ++i) {
    int k = obs(i)-1;
    epsilon(0) = (y(i,0)-mu(k,0))/w(i,0);
    epsilon(1) = (y(i,1)-mu(k,1))/w(i,1);
    nll += nll_dist(epsilon);
  }

  nll += bshrink(0)*beta(0)*beta(0)+bshrink(1)*beta(1)*beta(1);

  ADREPORT(beta)
  ADREPORT(sigma)
  ADREPORT(tau)
  return nll;
}


