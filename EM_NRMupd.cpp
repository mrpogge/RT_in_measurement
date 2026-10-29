#include <Rcpp.h>
using namespace Rcpp;

// helper function to compute a determinant of a 3x3 matrix
double det3x3(NumericMatrix B) {
  return B(0,0)*(B(1,1)*B(2,2) - B(1,2)*B(2,1))
  - B(0,1)*(B(1,0)*B(2,2) - B(1,2)*B(2,0))
  + B(0,2)*(B(1,0)*B(2,1) - B(1,1)*B(2,0));
}

// helper function to invert a 3x3 matrix
NumericMatrix invert3x3(NumericMatrix B) {
  NumericMatrix C(3,3);
  double d = det3x3(B);
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      NumericMatrix sub(2,2);
      int r = 0;
      for (int a = 0; a < 3; ++a) {
        if (a == i) continue;
        int c = 0;
        for (int b = 0; b < 3; ++b) {
          if (b == j) continue;
          sub(r,c++) = B(a,b);
        }
        r++;
      }
      double minor = sub(0,0)*sub(1,1) - sub(0,1)*sub(1,0);
      C(i,j) = std::pow(-1, i + j) * minor / d;
    }
  }
  NumericMatrix result(3,3);
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      result(i,j) = C(j,i);
  return result;
}

// helper function to compute a determinant of a 2x2 matrix
double det2x2(NumericMatrix B){
  return( B(0,0)*B(1,1)-B(0,1)*B(1,0));
}


// helper function to invert a 2x2 matrix
NumericMatrix invert2x2(NumericMatrix B){
  NumericMatrix C(2,2);
  double d=det2x2(B);
  C(0,0)=B(1,1)/d;
  C(1,1)=B(0,0)/d;
  C(0,1)=-B(0,1)/d;
  C(1,0)=C(0,1);
  return(C);
}

// helper function for the Newton-Raphton algorithm
NumericVector computeDelta(List Blocks, NumericVector Jis, double Js, NumericVector G, int n_par) {
  int n = Blocks.size();
  int dim = n_par * n;
  
  // Step 0: Precompute all inverted blocks
  std::vector<NumericMatrix> invBlocks(n);
  for (int k = 0; k < n; ++k) {
    invBlocks[k] = (n_par == 3) ? invert3x3(Blocks[k]) : invert2x2(Blocks[k]);
  }
  
  // Step 1: Compute u = invD %*% Jis using invBlocks
  NumericVector u(dim);
  for (int k = 0; k < n; ++k) {
    NumericMatrix inv = invBlocks[k];
    for (int i = 0; i < n_par; ++i)
      for (int j = 0; j < n_par; ++j)
        u[n_par * k + i] += inv(i, j) * Jis[n_par * k + j];
  }
  
  // Step 2: Sc = Js - t(Jis) %*% u
  double sum = 0.0;
  for (int i = 0; i < dim; ++i) {
    sum += Jis[i] * u[i];
  }
  double Sc = Js - sum;
  
  // Step 3: Compute delta = -G %*% D without building D
  NumericVector delta(dim + 1);
  
  for (int j = 0; j < dim; ++j) {
    double temp = 0.0;
    
    // invD contribution
    int block_id = j / n_par;
    int col = j % n_par;
    NumericMatrix inv = invBlocks[block_id];
    
    for (int row = 0; row < n_par; ++row) {
      int i = block_id * n_par + row;
      temp += G[i] * inv(row, col);
    }
    
    // rank-1 update contribution: (u %*% t(u)) / Sc
    for (int i = 0; i < dim; ++i) {
      temp += G[i] * (u[i] * u[j]) / Sc;
    }
    
    delta[j] = -temp;
  }
  
  // Last element of delta (scalar parameter)
  double last = 0.0;
  for (int i = 0; i < dim; ++i) {
    last += G[i] * (-u[i] / Sc);
  }
  last += G[dim] * (1.0 / Sc);
  delta[dim] = -last;
  
  return delta;
}

// [[Rcpp::export]]


// EM algorithm for the models with 2 or 3 item parameters
List EMnrm_npar_cpp(
    IntegerVector Position,    // for each person the first position for their responses
    IntegerVector I,           // item index
    NumericVector S,           // item scores
    NumericVector Weight,      // weights of the persons
    NumericVector th,          // nodes (K)
    NumericVector x,           // starting values
    NumericVector w,           // log-weights for nodes
    int N,                     // number of persons
    int n,                     // number of items
    int K,                     // number of nodes
    int max_iterEM,             // maximum number of iterations for the full EM   
    int max_iterNR, // maximum number of iterations with NR in the Mstep
    NumericVector scores, // scores corresponding to the different response categories
    NumericVector RR, // length n_par*n sufficient statistics for the intercepts
    double tol,               // maximum difference from previous iteration  
    int n_par) {              // number of item parameters
  
  NumericMatrix DIF1(K,n);
  NumericVector delta(n*n_par + 1);
  NumericMatrix DIF2(K,n);
  NumericMatrix DIF3(K,n);
  NumericMatrix lDIF(K,n);
  // Initialize minus gradient
  NumericVector G(n_par*n + 1);
  // Allocate Jacobian
  NumericVector Jis(n_par*n);
  NumericVector L(K);
  NumericMatrix W(K, n);
  NumericMatrix R(K, n);
  double Js=0.0;
  List Blocks(n);
  NumericVector gr_i(n_par);
  NumericVector jac_i(n_par);
  NumericVector jac_ij = (n_par == 3) ? NumericVector(3) : NumericVector(1);
  NumericVector jac_is(n_par);
  NumericVector eap(N);
  NumericVector post_var(N);
  
  int iter = 0;
  for (; iter < max_iterEM; ++iter) {
    
    // --- E-step ---
    NumericVector x_old=x;
    
    // Step 0: Pre-computations for each combination of node and item
    for(int k=0;k<K;++k){
      double th_k=th[k];
      double slope=x[n_par*n];
      for(int i=0;i<n;++i){
        double dif1=scores[0]*slope*th_k+x[i*n_par];
        DIF1(k,i)=dif1;
        double dif2=scores[1]*slope*th_k+x[1+i*n_par];
        DIF2(k,i)=dif2;
        if(n_par==3){
          double dif3=scores[2]*slope*th_k+x[2+i*n_par];
          DIF3(k,i)=dif3;
          lDIF(k,i)=log1p(exp(dif1)+exp(dif2)+exp(dif3));
        }
        if(n_par==2){
          lDIF(k,i)=log1p(exp(dif1)+exp(dif2));
        }
      }
    }
    
    
    int start=0;
    int end=0;
    W.fill(0.0); // sets all entries of W to zero
    R.fill(0.0); // sets all entries of R to zero
    

// Compute the posterior distributions
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
          double dif1 = DIF1(k,i);  
          double dif2 = DIF2(k,i);  
          double dif3 = DIF3(k,i);// if n_par==2 then it will be 0  
          double ldif=lDIF(k,i);	
          L[k] += dif1*(s==scores[0])+dif2*(s==scores[1])-ldif;
          if(n_par==3){L[k]+=dif3*(s==scores[2]);}
          
        }
      }
      
      // step 3: Normalize the posterior with softmax
      // find the maximum value
      double max_L = L[0];
      for (int k = 1; k < K; ++k) {
        if (L[k] > max_L) max_L = L[k];
      }
      // Compute exponentials and the normalising constant
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
    
    // W and R are input for the M-step
    
    // --- M-step (Newton-Raphson) ---
    for(int iterNR=0;iterNR<max_iterNR;++iterNR){
      
      double slope = x[n_par*n];  // s
      
      // Compute the gradient and the Jacobian
      
      double gr_s = 0.0;
      Js = 0.0;
      
      for (int i = 0; i < n; ++i) {
        gr_i.fill(0.0);
        jac_i.fill(0.0);
        jac_ij.fill(0.0);
        jac_is.fill(0.0);
        
        for (int j = 0; j < n_par; ++j) {
          gr_i[j] = RR[i + j * n];
        }
        
        for (int k = 0; k < K; ++k) {
          double t = th[k];
          NumericVector eD(n_par), e(n_par);
          
          double Z=1;
          for (int j = 0; j < n_par; ++j) {
            eD[j] = exp(scores[j] * slope * t + x[j + i * n_par]);
            Z+=eD[j];
          }
          
          for (int j = 0; j < n_par; ++j) {
            e[j] = eD[j] / Z;
          }
          
          double expected = 0.0, expected2 = 0.0;
          for (int j = 0; j < n_par; ++j) {
            expected += scores[j] * e[j];
            expected2 += scores[j] * scores[j] * e[j];
          }
          
          double ww = W(k, i);
          for (int j = 0; j < n_par; ++j) {
            gr_i[j] -= ww * e[j];
            jac_i[j] -= ww * e[j] * (1 - e[j]);
            jac_is[j] += t * ww * e[j] * (expected - scores[j]);
          }
          
          gr_s += t * (R(k, i) - ww * expected);
          Js -= t * t * ww * (expected2 - expected * expected);
          
          jac_ij[0] += ww * e[0] * e[1]; // 12
          if(n_par==3){
            jac_ij[1] += ww * e[0] * e[2]; // 13
            jac_ij[2] += ww * e[1] * e[2]; // 23
          }
          
        }
        
        for (int j = 0; j < n_par; ++j) {
          G[j + i * n_par] = gr_i[j];
          Jis[j + i * n_par] = jac_is[j];
        }
        
        
        NumericMatrix block(n_par, n_par);
        for(int j=0; j<n_par; ++j){
          block(j,j)=jac_i[j];
        }
        block(0, 1) = jac_ij[0];
        block(1, 0) = jac_ij[0];
        if(n_par==3){
          block(0, 2) = jac_ij[1];
          block(1, 2) = jac_ij[2];
          block(2, 0) = jac_ij[1]; 
          block(2, 1) = jac_ij[2]; 
        }
        Blocks[i] = block;
      }
      
      G[n_par * n] = gr_s;
      
      //Solve delta=-(J)^(-1)*G using the fact that the Jacobian is blocked bordered diagonal
      
      delta=computeDelta(Blocks,Jis,Js,G,n_par);
      
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