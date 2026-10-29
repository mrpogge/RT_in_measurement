#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List elo_cpp(
    IntegerVector position,// vector specifying where the data from wach person starts and ends
    IntegerVector item_ids,// item IDs for all responses
    IntegerVector X, // response accuracy data
    NumericVector S, // item scores
    NumericVector delta,// item parameter estimates
    int N, // number of persons
    double K, // K-value
    NumericVector z,// vector of item scores [for PCM and NRM]
    double r, // penalty P in for the CRTM
    int model){ // model indicator [0 - RM, 1 - SRT, 2 - DSRT, 3 - CRT, 4 - PCM, 5 - NRM]
  
  double MSE=0.0;//
  double theta=0;
  double MSEp=0.0;
  
  for (int p = 0; p < N; ++p) {
    MSEp=0.0;// computing person specific MSE
    theta=0;
    int start=position[p];
    int end=position[p+1];
    for(int i=start;i<end;++i){// loop over responses of person p
      int j = item_ids[i]; 
      double s=S[i];
      int x=X[i];
      double Ex=0;
      double Es=0;

// for each model compute the expected value for the accuracy and the expected score
      switch(model){
      case 0:{
        double D=exp(theta+delta[j]);
        Ex=D/(1+D);
        Es=Ex;
        break;
      }
      
      case 1:{
        double D=(theta+delta[j]);
        double eD=exp(D);
        double eD2=eD*eD;
        Ex=eD/(1+eD);
        Es=(eD2+1)/(eD2-1)-1/D;
        break;
      }
        
      case 2:{
        double D = theta + delta[j];
        double D2 = exp(D);           // exp(D)
        double D1 = sqrt(D2);         // exp(0.5 * D)
        double D3 = D2 * D1;    // exp(1.5 * D)
        Ex=D2/(1+D2);
        Es=(0.5*D1+D2+1.5*D3)/(1+D1+D2+D3);
        break;
      }  
        
      case 3:{
        double D=theta+delta[j];
        double eD=exp(D);
        double DerD=D*exp(-r*D);
        Ex=(eD-1)/(eD-1+DerD);
        Es=(-r*D*DerD+(D-1)*eD+1)/(D*(eD-1+DerD));
        break;
      }
        
      case 4:{
        double eT = exp(theta);
        double sqrt_eT = sqrt(eT);  // equivalent to exp(0.5 * D)
        double eD1 = exp(delta[j*3]);
        double eD2 = exp(delta[j*3+1]);
        double eD3 = exp(delta[j*3+2]);
        double p1 = sqrt_eT*eD1;
        double p2 = eT * eD2;
        double p3 = eT * sqrt_eT * eD3;
        double denom = 1.0 + p1+p2+p3;
        Ex = (p2+p3) / denom;
        Es = (0.5 * p1 + p2 + 1.5 * p3) / denom;
        break;
      }  
      
      case 5:{
        double e1=exp(z[0]*theta+delta[j*2]);
        double e2=exp(z[1]*theta+delta[j*2+1]);
        double denom=1+e1+e2;
        Ex=(e1+e2)/denom;
        Es=(z[0]*e1+z[1]*e2)/denom;
        break;
      }
     
    }
      double change = K * (s-Es); // elo update
      theta += change;
      double error=(x-Ex)*(x-Ex); // squared error
      MSEp+=error;// adding up all errors for the specific person
      S[i]=theta;// we re-use the vector with the item scores to save the Elo ratings
    }
    MSEp/=(end-start); // divide the MSEp by the number of responses the person had
    MSE+=MSEp;// adding up all person-specific MSEs
  }
  MSE/=N; // divide the total MSE by the number of persons
  return List::create(
    Named("MSE") = MSE,
    Named("MSEp") = MSEp,
    Named("Theta")=S
  );
}


