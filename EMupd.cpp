#include <Rcpp.h>

using namespace Rcpp;


// [[Rcpp::export]]

// Function for the EM algorithm of the models with only 1 item parameter
List EMupd_cpp(
    IntegerVector Position,    // for each person the position at which their data start (length N+1, the last element is equal to n_resp-1)
    IntegerVector I,           // item index
    NumericVector S,           // item scores
    NumericVector Weight,      // weights of the persons
    NumericVector th,          // nodes (K)
    NumericVector x,           // starting values
    NumericVector w,           // log-weights for nodes
    int N,                     // number of persons
    int n,                     // number of items
    int K,                     // number of nodes
    int model,                 // model 0 - Rasch, 2 - DSRT
    double r,                  // penalty for incorrect responses
    int max_iterEM,             // maximum number of iterations for the full EM   
    int max_iterNR,             // maximum number of iterations with NR in the Mstep
    double tol) {              // maximum difference from previous iteration
  
  NumericMatrix DIF(K,n);
  NumericMatrix lDIF(K,n);
  // Initialize minus gradient
  NumericVector G(n + 1);
  // Allocate Jacobian
  NumericVector Ji(n);
  NumericVector Jis(n);
  NumericVector L(K);
  NumericMatrix W(K, n);
  NumericMatrix R(K, n);
  NumericVector delta(n + 1);
  NumericVector eap(N);
  NumericVector post_var(N);
  
  int iter=0;
  for (; iter < max_iterEM; ++iter) {
    
    // --- E-step ---
    NumericVector x_old=x;
    
    // Step 0: Pre-computations for each combination of node and item of the log-likelihood
    if(model==0){
      for(int k=0;k<K;++k){
        double th_k=th[k];
        double slope=x[n];
        for(int i=0;i<n;++i){
          double dif=slope*th_k+x[i];
          DIF(k,i)=dif;
          lDIF(k,i)=log1p(exp(dif));
        }
      }
    }
    if(model==1){
      for(int k=0;k<K;++k){
        double th_k=th[k];
        double slope=x[n];
        for(int i=0;i<n;++i){
          double dif=slope*th_k+x[i];
          DIF(k,i)=dif;
          lDIF(k,i)=log((exp(dif)-exp(-dif))/(dif));
        }
      }
    }
    if(model==2){
      for(int k=0;k<K;++k){
        double th_k=th[k];
        double slope=x[n];
        for(int i=0;i<n;++i){
          double dif=slope*th_k+x[i];
          DIF(k,i)=dif;
          lDIF(k,i)=log1p(exp(0.5*dif)+exp(dif)+exp(1.5*dif));
        }
      }
    }
    if(model==3){
      for(int k=0;k<K;++k){
        double th_k=th[k];
        double slope=x[n];
        for(int i=0;i<n;++i){
          double dif=slope*th_k+x[i];
          DIF(k,i)=dif;
          lDIF(k,i)=log((exp(dif)-1+dif*exp(-r*dif))/dif);
        }
      }
    }
    
    int start=0;
    int end=0;
    W.fill(0.0); // sets all entries of W to zero
    R.fill(0.0); // sets all entries of R to zero
    
// Compute the posterior distribution for each person
    for(int j=0;j<N;++j){
      double weight = Weight[j];
      start=Position[j];
      end=Position[j+1];
      // step 1: Take the log-prior
      L=clone(w);

      // step 2: Add log-likelihoods of each response
      for(int m=start;m<end;++m){
        int i=I[m];
        double s=S[m];
        for(int k=0;k<K;++k){
          double dif = DIF(k,i);  
          double ldif=lDIF(k,i);	
          L[k] += s * dif - ldif;
        }
      }
      
      // step 3: Normalize the posterior with softmax
      // find the maximum value
      double max_L = L[0];
      for (int k = 1; k < K; ++k) {
        if (L[k] > max_L) max_L = L[k];
      }
      // Compute exponentials and sum
      double sum_exp = 0.0;
      for (int k = 0; k < K; ++k) {
        double val = exp(L[k] - max_L);
        L[k] = val;
        sum_exp += val;
      }
      // Normalize
      for (int k = 0; k < K; ++k) {
        L[k] /= sum_exp;
      }
      
      //Step 4: accumulate W and R
      for (int m = start; m < end; ++m) {
        int i = I[m];
        double sw=S[m]*weight;// weighted scores of person j
        for (int k = 0; k < K; ++k) {
          double post = L[k];
          W(k, i) += post*weight; 
          R(k, i) += post*sw;  
        }
      }
      
      // Compute the posterior mean and variance
      double av=0;
      double av2=0;
      for(int k=0;k<K;++k){
        double t=th[k];
        double WW=L[k];
        av+=t*WW;
        av2+=t*t*WW;
      }
      eap[j]=av;
      post_var[j]=av2-av*av;
    }
    
    // --- M-step (Newton-Raphson) ---
    for(int iterNR=0;iterNR<max_iterNR;++iterNR){
      
      double slope = x[n];  // s
      
      // Compute the gradient and the Jacobian
      double gr_s = 0.0;
      double Js=0.0;
      if (model == 0) {
        // Rasch model: logistic response function
        for (int i = 0; i < n; ++i) {
          double gr_i=0.0;
          double jac_i=0.0;
          double jac_is=0.0;
          for (int k = 0; k < K; ++k) {
            double t=th[k];
            double D=slope * t + x[i];
            double expected=1.0 / (1.0 + std::exp(-D));
            double var = expected * (1.0 - expected);  // variance of Bernoulli
            double wp=R(k, i) - W(k, i) * expected;
            gr_i -=wp;
            gr_s-=t*wp;
            Js-= t * t * W(k, i) * var;
            jac_is -= t * W(k, i) * var;
            jac_i -= W(k, i) * var;
          }
          G[i]=gr_i;
          Jis[i]=jac_is;
          Ji[i]=jac_i;
        }
      }
      if (model == 1) {
        for (int i = 0; i < n; ++i) {
          double gr_i=0.0;
          double jac_i=0.0;
          double jac_is=0.0;
          for (int k = 0; k < K; ++k) {
            double t=th[k];
            double D=slope * t + x[i];
            double eD1=exp(D);
            double eD2=exp(-D);
            double expected=((D-1)*eD1+(D+1)*eD2)/(D*(eD1-eD2));
            double D2=D*D;
            double se2=(eD1*(1-2/D+2/D2)-eD2*(1+2/D+2/D2))/(eD1-eD2);
            double var = se2-expected*expected;  // variance is E(S^2|D)-(E(S|D))^2;
            double wp=R(k, i) - W(k, i) * expected;
            gr_i -=wp;
            gr_s-=t*wp;
            Js-= t * t * W(k, i) * var;
            jac_is -= t * W(k, i) * var;
            jac_i -= W(k, i) * var;
          }
          G[i]=gr_i;
          Jis[i]=jac_is;
          Ji[i]=jac_i;
        }
      }
      if (model == 2) {
        // DSRT model
        std::vector<double> z = {0.0, 0.5, 1.0, 1.5};
        int M = z.size();
        for (int i = 0; i < n; ++i) {
          double gr_i=0.0;
          double jac_i=0.0;
          double jac_is=0.0;
          for (int k = 0; k < K; ++k) {
            double t=th[k];
            double D=slope * t + x[i];
            double expected=0.0;
            double var=0.0;
            std::vector<double> probs(M);
            double denom = 0.0;
            for (int m = 0; m < M; ++m) {
              probs[m] = std::exp(z[m] * D);
              denom += probs[m];
            }
            for (int m = 0; m < M; ++m) {
              probs[m] /= denom;
              expected += z[m] * probs[m];
            }
            for (int m = 0; m < M; ++m) {
              var += probs[m] * std::pow(z[m] - expected, 2);
            }
            double wp=R(k, i) - W(k, i) * expected;
            gr_i-=wp;
            gr_s-=t*wp;
            Js-= t * t * W(k, i) * var;
            jac_is -= t * W(k, i) * var;
            jac_i -= W(k, i) * var;
          }
          G[i]=gr_i;
          Jis[i]=jac_is;
          Ji[i]=jac_i;
        }
      }
      if (model == 3) {
        for (int i = 0; i < n; ++i) {
          double xx=x[i];
          double gr_i=0.0;
          double jac_i=0.0;
          double jac_is=0.0;
          for (int k = 0; k < K; ++k) {
            double t=th[k];
            double D=slope * t + xx;
            double eD=exp(D);
            double DerD=D*exp(-r*D);
            double D2=D*D;
            double expected=((D-1)*eD+1-r*D*DerD)/(D*(eD-1+DerD));
            double se2=(((D2-2*D+2)*eD-2)/D2+r*r*DerD)/(eD-1+DerD);
            double var=se2-expected*expected;
            double wp=R(k, i) - W(k, i) * expected;
            gr_i -=wp;
            gr_s-=t*wp;
            Js-= t * t * W(k, i) * var;
            jac_is -= t * W(k, i) * var;
            jac_i -= W(k, i) * var;
          }
          G[i]=gr_i;
          Jis[i]=jac_is;
          Ji[i]=jac_i;
        }
      }
      
      G[n]=gr_s;
      
      // Solve delta=(J)^(-1)*G using the fact that the Jacobian is bordered diagonal
      
      // Step 1: Solve D * x1 = g1
      Rcpp::NumericVector x1(n);
      for (int i = 0; i < n; ++i) {
        x1[i] = G[i] / Ji[i];
      }
      
      // Step 2: Compute inner sums
      double sum_vx1 = 0.0;
      double sum_v2_over_D = 0.0;
      for (int i = 0; i < n; ++i) {
        sum_vx1 += Jis[i] * x1[i];
        sum_v2_over_D += (Jis[i] * Jis[i]) / Ji[i];
      }
      
      // Step 3: Solve for delta[n]
      delta[n] = (G[n] - sum_vx1) / (Js - sum_v2_over_D);
      
      // Step 4: Update delta[0 to n-1]
      for (int i = 0; i < n; ++i) {
        delta[i] = x1[i] - (Jis[i] * delta[n]) / Ji[i];
      }
      
      
      NumericVector x_new = x + delta;
      
      double max_deltaNR = max(abs(delta));
      if (max_deltaNR < tol) {
        break;
      }
      
      x = x_new;
      
      
      
    }
    
    NumericVector delta_outer=x_old-x;
    
    double max_delta=max(abs(delta_outer));
    
    if (max_delta < tol) {
      break;
    }
    
  }  
  
  return Rcpp::List::create(
    Rcpp::Named("par") = x,
    Rcpp::Named("iter") = iter,
    Rcpp::Named("EAP")=eap,
    Rcpp::Named("Post_var")=post_var
  );
  
  
  
}