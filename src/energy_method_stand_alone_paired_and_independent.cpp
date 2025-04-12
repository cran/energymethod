#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;


double MIXED_STAT(int K, int T, arma::vec data){
  
  arma::mat reshaped_data=arma::reshape(data, T, K);
  arma::rowvec col_norms_squared = arma::sum(arma::square(reshaped_data), 0);
  return(arma::max(col_norms_squared)/((double )T));
  
}
/*

double MIXED_STAT(int K, int T, const arma::vec& data) {
  arma::mat reshaped_data(data.memptr(), T, K, false, true);
  return arma::max(arma::sum(arma::square(reshaped_data), 0)) / static_cast<double>(T);
}
 */

double compute_p_val(double test_stat, const arma::vec& empirical_cdf){
  double denom = (double)empirical_cdf.n_elem;
  //double rank=(double)(arma::index_min(arma::abs(empirical_cdf - test_stat)));
  
  arma::vec::const_iterator it = std::lower_bound(empirical_cdf.begin(), empirical_cdf.end(), test_stat);
  return (std::distance(it, empirical_cdf.end())/denom);
}


// [[Rcpp::export]]
List energy_method(arma::cube sample_1, arma::cube sample_2, int num_bootstrap_reps,int seed, std::string type) {
  
    Rcpp::RNGScope scope;
    Environment base_env("package:base");
    Function set_seed_r = base_env["set.seed"];
    set_seed_r(seed);

  
  if(sample_1.n_rows!=sample_2.n_rows){
    stop("The two samples must have the same number of channels");
  }
  if(sample_1.n_cols!=sample_2.n_cols){
    stop("The two samples must have the same number of recordings per channel");
  }
  
  if(type=="paired"){
    if(sample_1.n_slices!=sample_2.n_slices){
      stop("The two samples must each have the same sample size for a paired test");
    }
  }
  

  int K = sample_2.n_rows;
  int T=sample_2.n_cols;
//stores the bootstrap replicates for estimating c.d.f
  arma::vec bootstrap_replicates(num_bootstrap_reps);
  //stores the global test_stat
  
  double global_test_stat=0;
  //stores test stat for each of the K channels
  arma::vec channel_wise_test_stat = arma::vec(K, arma::fill::zeros);  
  //stores the p-value for each of the K channels
  arma::vec channel_wise_p_val = arma::vec(K, arma::fill::zeros); 
  //stores the global p_val
  double global_p_val=1;
  
  //the sample sizes of the two samples; if the test is paired these are required to be the same, and 
  //in the code for the paired test below, a common sample_size int is used.
  int sample_size_1=sample_1.n_slices;
  int sample_size_2=sample_2.n_slices;
 
  if(type=="paired"){
    //confirmed with stop condition that the number of slices is the same across both samples, so 
    //the choice to set sample_size according to sample_2 is arbitrary; same with K and T.
    int sample_size=sample_size_2;
    //int K = sample_2.n_rows;
    //int T=sample_2.n_cols;
    //construct differente of two samples. 
  arma::cube differenced_sample=sample_1-sample_2;
  //clear memory as we no longer need original sample_1 and sample_2
  sample_1.reset();
  sample_2.reset();
//take difference of two samples for used in paired test
  arma::mat sample(sample_size, K*T);
  for (int i = 0; i < sample_size; i++) { 
    sample.row(i) = arma::vectorise(differenced_sample.slice(i).t()).t(); // Transpose to make it a row
  }
  //now free-up memory from structured differenced_sample as it is no longer needed
  differenced_sample.reset();
  
  //take column mean of the difference vector; test is based on this mean
  arma::vec col_means(K*T);
  col_means=mean(sample, 0).t();
  
 
  arma::mat sample_de_meaned(sample_size, K*T);
  sample_de_meaned= sample.each_row()- col_means.t();
  
  
  arma::mat sample_de_meaned_multipliers(sample_size, K*T);
  arma::vec multipliers(sample_size);
  
  arma::vec col_means_bootstrap(K*T);
 


  
  
  
  
  
  
  
  for(int b=0; b<num_bootstrap_reps; b++){
    multipliers= arma::randn(sample_size);
    sample_de_meaned_multipliers=sample_de_meaned.each_col() % multipliers;
    col_means_bootstrap=mean(sample_de_meaned.each_col() % multipliers,0).t();
    
    bootstrap_replicates(b)=((double)sample_size)*MIXED_STAT(K, T, col_means_bootstrap);
  
  }
 bootstrap_replicates=sort(bootstrap_replicates);

  
  global_test_stat=((double)sample_size)*MIXED_STAT(K,T,col_means);
 
  
  global_p_val =compute_p_val(global_test_stat,bootstrap_replicates);
 
  for(int k=0; k <K; k++){
    channel_wise_test_stat(k)=((double)sample_size)*MIXED_STAT(1,T,col_means.subvec(k*T, (k+1)*T-1));
  
    
    channel_wise_p_val(k) = compute_p_val(channel_wise_test_stat(k),bootstrap_replicates);
  
  }

  }
  
  if(type=="independent"){
    int sample_size_1=sample_1.n_slices;
    int sample_size_2=sample_2.n_slices;
    //stop condition requires that K and T are the same for both samples, so chose sample_2 arbitrarily:
    int K = sample_2.n_rows;
    int T=sample_2.n_cols;
    
    double fac = ((double)sample_size_1)*((double)sample_size_2)/(((double)sample_size_1)+((double)sample_size_2));
    //double sqrt_ss_1=std::sqrt(sample_size_1);
    //double sqrt_ss_2=std::sqrt(sample_size_2);
    
    arma::mat sample_1_f(sample_size_1, K*T);
    arma::mat sample_2_f(sample_size_2, K*T);
  

  
    
    
    for (int i = 0; i < sample_size_1; i++) { 
      sample_1_f.row(i) = arma::vectorise(sample_1.slice(i).t()).t(); // Transpose to make it a row
    }
    for (int i = 0; i < sample_size_2; i++) { 
      sample_2_f.row(i) = arma::vectorise(sample_2.slice(i).t()).t(); // Transpose to make it a row
    }
    
    sample_1.reset();
    sample_2.reset();
    
    arma::vec col_means_1(K*T);
    arma::vec col_means_2(K*T);
    
    arma::vec col_means_differenced(K*T);
    
    col_means_1=mean(sample_1_f, 0).t();
    col_means_2=mean(sample_2_f, 0).t();
    
    col_means_differenced = col_means_1-col_means_2;
    
    arma::mat sample_1_f_de_meaned(sample_size_1, K*T);
    arma::mat sample_2_f_de_meaned(sample_size_2, K*T);
    sample_1_f_de_meaned= sample_1_f.each_row()-col_means_1.t();
    sample_2_f_de_meaned= sample_2_f.each_row()-col_means_2.t();
    
    arma::mat sample_1_f_de_meaned_multipliers(sample_size_1, K*T);
    arma::mat sample_2_f_de_meaned_multipliers(sample_size_2, K*T);
    arma::vec multipliers_1(sample_size_1);
    arma::vec multipliers_2(sample_size_2);
    
    arma::vec col_means_1_bootstrap(K*T);
    arma::vec col_means_2_bootstrap(K*T);

    
    for(int b=0; b<num_bootstrap_reps; b++){
      multipliers_1= arma::randn(sample_size_1);
      multipliers_2= arma::randn(sample_size_2);
      sample_1_f_de_meaned_multipliers=sample_1_f_de_meaned.each_col() % multipliers_1;
      sample_2_f_de_meaned_multipliers=sample_2_f_de_meaned.each_col() % multipliers_2;
      
      col_means_1_bootstrap=mean( sample_1_f_de_meaned_multipliers,0).t();
      col_means_2_bootstrap=mean( sample_2_f_de_meaned_multipliers,0).t();
      
      bootstrap_replicates(b)=fac*MIXED_STAT(K, T, col_means_1_bootstrap-col_means_2_bootstrap);
      
    }
    global_test_stat=fac*MIXED_STAT(K,T,col_means_differenced);
    
    
    global_p_val =compute_p_val(global_test_stat,bootstrap_replicates);
    
    for(int k=0; k <K; k++){
      channel_wise_test_stat(k)=fac*MIXED_STAT(1,T,col_means_differenced.subvec(k*T, (k+1)*T-1));
      
      
      channel_wise_p_val(k) = compute_p_val(channel_wise_test_stat(k),bootstrap_replicates);
      
    }
    
    
  }
  
  
  
  
  List out;
  
  out["sample_size_1"]=sample_size_1;
  out["sample_size_2"]=sample_size_2;
  out["number_of_channels"] = K;
  out["number_of_recordings_per_channel"]=T;
  out["global_p_val"]=global_p_val;
  out["channel_wise_p_val"]=channel_wise_p_val;

  return(out);
  
  
}
