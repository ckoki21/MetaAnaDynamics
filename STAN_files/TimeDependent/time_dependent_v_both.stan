// see https://khakieconomics.github.io/2018/02/24/Regime-switching-models.html
//metaphase model. based on armond et al 2015 plos comp biol

//This version includes updated prior for L based on nocodazole data
// Also includes angle to the metaphase plate (input as data measured via KiT) 
//This adjusts the spring based forces which may not act perpendicular to metaphase plate

functions {
 matrix constructTransitionMatrix(real p_icoh, real p_coh){
    matrix[4,4] P;
    real q_coh;
    real q_icoh;
    q_icoh = 1-p_icoh;
    q_coh = 1-p_coh;
    P[1,1] = p_icoh*p_icoh;
    P[2,1] = p_coh*q_coh;
    P[3,1] = p_coh*q_coh;
    P[4,1] = q_icoh*q_icoh;
    P[1,2] = p_icoh*q_icoh;
    P[2,2] = p_coh*p_coh;
    P[3,2] = q_coh*q_coh;
    P[4,2] = p_icoh*q_icoh;
    P[1,3] = p_icoh*q_icoh;
    P[2,3] = q_coh*q_coh;
    P[3,3] = p_coh*p_coh;
    P[4,3] = p_icoh*q_icoh;
    P[1,4] = q_icoh*q_icoh;
    P[2,4] = p_coh*q_coh;
    P[3,4] = p_coh*q_coh;
    P[4,4] = p_icoh*p_icoh;
    return P;
  }
  matrix odeUpdateMatrix(real[] th, int t, real dt){
    matrix[2,6] M;
    M[1,1] = -th[3]- th[2];
    M[1,2] = th[3];
    M[1,3] = -th[5]*exp(th[10]*(t-1)*dt);
    M[1,4] = -th[5]*exp(th[10]*(t-1)*dt);
    M[1,5] = -th[4]*exp(th[9]*(t-1)*dt);
    M[1,6] = -th[4]*exp(th[9]*(t-1)*dt);
    M[2,1] = th[3];
    M[2,2] = -th[3] - th[2];
    M[2,3] = th[5]*exp(th[10]*(t-1)*dt);
    M[2,4] = th[4]*exp(th[9]*(t-1)*dt);
    M[2,5] = th[5]*exp(th[10]*(t-1)*dt);
    M[2,6] = th[4]*exp(th[9]*(t-1)*dt);
    return M;
  }
  vector odeUpdateVector(real[] th, real angleTheta){
    vector[2] mu;
    mu[1] = th[3]*th[8]*angleTheta;
    mu[2] = -th[3]*th[8]*angleTheta;
    return mu;
  }
  vector construct_transition_column(real p_icoh, real p_coh,int sigma){
    vector[4] P_col;
    real q_icoh = 1-p_icoh;
    real q_coh = 1-p_coh;
    if (sigma==1){
    P_col = [p_icoh*p_icoh,
             p_coh*q_coh,
             p_coh*q_coh,
             q_icoh*q_icoh]';
    } else if(sigma==2){
    P_col = [p_icoh*q_icoh,
              p_coh*p_coh,
              q_coh*q_coh,
              p_icoh*q_icoh]';
    } else if (sigma==3){
    P_col = [p_icoh*q_icoh,
              q_coh*q_coh,
              p_coh*p_coh,
              p_icoh*q_icoh]';
    } else if (sigma==4){
    P_col = [q_icoh*q_icoh,
              p_coh*q_coh,
              p_coh*q_coh,
              p_icoh*p_icoh]';
      } else {
      P_col=rep_vector(0,4);
    }
    return(P_col);
    }
}
data {
  int<lower=1> Frames;
  int<lower=0> T0;
  int<upper=Frames> T1;
  matrix[Frames,2] y; //trajectory data from each sister
  real dt;    
  int nStates; 
  vector[nStates] sigma0;
  vector[Frames] cos_theta; //the angle to metaphase plate is phi
}
transformed data {
    vector[2] x0;
  //  real b; //sharpness of sigmoid switch to anaphase. Fixed not inferred. 
    x0 = y[T0+1, ]';
  //  b = dt/2;
}
parameters {
  real<lower=0> tau;
  real<lower=0> alpha;
  real<lower=0> kappa;
  real<upper=0> v_minus_0;
  real<lower=0> v_plus_0;
  real<lower=0,upper=1> p_icoh;
  real<lower=0,upper=1> p_coh;
  real<lower=0> L;
  real v_minus_1;
  real v_plus_1;
}
transformed parameters{
  real theta[10]; //for metaphase ODE model
  matrix[Frames, nStates] eta;
  matrix[Frames, nStates] xi;
  vector[Frames] f;
  matrix[nStates,nStates] P;
  vector[nStates+2] auxStates;
  
  theta[1]=tau;
  theta[2]=alpha;
  theta[3]=kappa;
  theta[4]=v_minus_0;
  theta[5]=v_plus_0;
  theta[6]=p_icoh;
  theta[7]=p_coh;
  theta[8]=L;
  theta[9]=v_minus_1;
  theta[10]=v_plus_1;
  
for (t in 1:T0) {
  for (j in 1:nStates) {
    eta[t,j] = 0;
    xi[t,j] = 0;
  }
  f[t] = 1;
}
for (t in (T1+1):Frames) {
  for (j in 1:nStates) {
    eta[t,j] = 0;
    xi[t,j] = 0;
  }
  f[t] = 1;
}  
  // fill in etas
  for (t in (T0+1):T1) {
    if(t==(T0+1)) {
      for (j in 1:nStates){
        eta[t,j] = exp(normal_lpdf(y[t]| x0, sqrt(dt/tau)));
      }
    } else {
        if(sum(y[t-1,])>10^8 || sum(y[t,])>10^8){
           for (j in 1:nStates){
           eta[t,j] = 1;
            }
          }else{
            for (j in 1:4){
            //construct vector from appending y[t-1] and hidden state eg [0 1 0 0]
               for (i in 1:(nStates+2)){
        if (i<=2) {
          auxStates[i] = y[t-1,i];
        } else if (i==(j+2)) {
          auxStates[i] = 1.0;
        } else {
          auxStates[i] = 0.0;
        }
      }
        eta[t,j] = exp(normal_lpdf(y[t]' | y[t-1]' + dt*odeUpdateMatrix(theta, t-1, dt)*auxStates + dt*odeUpdateVector(theta, cos_theta[t]), sqrt(dt/theta[1])));
      }
    }
  }
}
  
  P = constructTransitionMatrix(p_icoh,p_coh);
  // work out likelihood contributions
  for(t in (T0+1):T1) {
  // for the first observation
  if(t==(T0+1)) {
    //replacing xi[t-1] by sigma0
    //For computational reasons I replaced the f[t] and xi with the following.
//    f[t] = dot_product((sigma0')*P, eta[t]');
     f[t]=(sum(exp(log(sigma0'*P)+log(eta[t]))));
    //xi[t] = (((sigma0')*P) .* eta[t]) ./ f[t];
    xi[t]=(exp((log(sigma0'*P)+log(eta[t]))-log(f[t])));
    } else {
    // and for the rest
    
    //f[t] = dot_product(xi[t-1]*P, eta[t]');
    f[t]=(sum(exp(log(xi[t-1]*P)+log(eta[t]))));
    //xi[t] = ((xi[t-1]*P) .* eta[t]) ./ f[t];
    xi[t]=(exp((log(xi[t-1]*P)+log(eta[t]))-log(f[t])));
    }
  }
}
model {
//   priors ...
//    tau ~ gamma(0.5,1.0/1000);
  target+= gamma_lpdf(tau|0.5,10^(-3));
//   alpha ~ normal(0.01,100) T[0,];
  target+= normal_lpdf(alpha|0.01,0.1)-normal_lccdf(0|0.01,0.1); //T[0,] DO I NEED THIS?
//   kappa ~ normal(0.05,100) T[0,];
//   target+= normal_lpdf(kappa|0.05,0.1)-normal_lccdf(0|0.05,0.1);
  target+= normal_lpdf(kappa|0.05,0.1)-normal_lccdf(0|0.05,0.1); //is this better? Variances are markedtly different than the ones at Ar15
//   v_minus ~ normal(-0.03,0.1) T[,0];
  target+= normal_lpdf(v_minus_0|-0.03,0.1)-normal_lcdf(0|-0.03,0.1);
//   v_plus ~ normal(0.03,0.1) Tr[0,];
  target+=normal_lpdf(v_plus_0|0.03,0.1)-normal_lccdf(0|0.03,0.1);  
//  p_icoh ~ beta(2,1);
  target+=beta_lpdf(p_icoh|2,1);
//  p_coh ~ beta(2.5,1);
  target+=beta_lpdf(p_coh|2.5,1);
//  L ~ normal(0.790,0.119) T[0,];
  target+= normal_lpdf(L|0.790,0.119)-normal_lccdf(0|0.790,0.119);
  target+= normal_lpdf(v_minus_1|0, log(10)/600);
  target+= normal_lpdf(v_plus_1|0, log(10)/600);
  target += sum(log(f));
}

 generated quantities {
   vector[Frames] log_lik;
//sampling from hidden states could go here
   int sigma_sim[Frames];
   vector[nStates] state_probs = xi[T1]';
   vector[nStates] conditional_state_probs; //on log scale
   vector[nStates] P_col;
   int frame;
    
  log_lik = rep_vector(0,Frames); // repeat the size T vector(column) consisting of copies of 0
   for (t in 1:Frames){
    log_lik[t] = log(f[t]);
   }
  
   for (t in (T1+1):Frames) {
     sigma_sim[t] = 0; //for time points with missing data
  }
   for (t in 1:T0) {
    sigma_sim[t] = 0; //for time points with missing data
   }


   //start with the final time pt
  sigma_sim[T1] = categorical_rng(state_probs);
  for (t in 1:(T1-T0-1)){
     frame = T1+1-t; //from T1 to T1+1-(T1-T0-1)=T0+2
    P_col = construct_transition_column(p_icoh,p_coh,sigma_sim[frame]);
     conditional_state_probs = log(P_col) + log(xi[frame-1]');
    sigma_sim[frame-1] =categorical_logit_rng(conditional_state_probs);
  }
    for (t in 1:Frames) {
    log_lik[t] = log(f[t]);
    
  }
}

