

clear all; close all; clc;
%%% NUCLEOTIDE ORIENATION:     %A = 1, C = 2, G = 3, T = 4
   
   

%%%%%%%%%%%% READ IN DATA %%%%%%%%%%%%

   Data_leaves = fastaread('data_ward.txt');
   Data_Matrix_Int = fn_Seq_to_Matrix(Data_leaves);
   
   
   
   subst_model = 3;
                       %%%% JC69 = 1
                       %%%% Kimura 2 parameters = 2
                       %%%% F84 = 3
                       %%%% HKY85 = 4

   
   
%%% UPGMA tree.... this is to get the first tree
%     UPGMAtree = fn_UPGMA_tree(Data_leaves,1);
%     h = plot(UPGMAtree,'orient','top');
%     title('UPGMA Distance Tree of Primates using Jukes-Cantor model');
%     ylabel('Evolutionary distance')

 %%% 2... Result from UPGMA    
%         UPGMA_arrangee = [6 2 3;...
%                           7 5 6;...
%                           8 1 7;...
%                           9 4 8];
                          
 


            
            
       
%%%%%%%%%%%% WATTERSON ESTIMATES %%%%%%%%%%%%
         %[a, b, c, d] = fn_Watterson(Data_Matrix_Int)
         
         

  
  
  
  
  
  
  
  
%%%%%%%%%%%% GENERAL VARIABLES %%%%%%%%%%%%


no_leaves = size(Data_Matrix_Int,1); 
n = no_leaves;
Data = Data_Matrix_Int;



N = 500; %100 %number of samples we supposed to have 
N_mcmc = 50; %10

%%% EPSILON
            be_ep = 0;  % beginning of epsilon
            en_ep = 1;  % end of epsilon
            int_ep = 0.01; % interval of epsilon

            Epsil = be_ep:int_ep:en_ep;
            T = length(Epsil);

%%% WEIGHTS
            wt = zeros(T,N);
            wt(1,:) = 1/N;   % assign equal weight to all 
            
            
                
%%%% TREE MATRIX HEAD
            abl = [1:n]';
            roco = zeros(n,6);
            roco(:,1) = abl;
            TREE_HEAD = roco; % designed so that we will save memory




            
            
switch subst_model
    
    
case 1    %%% JC69
        
               
%%%%%%%%%%%% CONSTRUCT THE FIRST ARRANGEMENT OF THE TREE %%%%%%%%%%%%
        
       

               
         
%%% Number of parameters to sample
           
       B = 2; 
           %%%% 1. tree
           %%%% 2. Theta
           
           
%%% STORAGE
            TREE_SAMPLES = zeros(n-1,6,N);
            THETA_SAMPLES = zeros(1,N);
           




%%% PRIORS DISTRIBUTION PARAMETERS
    %%% 1... THETA
        theta_max = 1.0;
        theta_min = 0.002;
        delta_theta = 0.5;  % less than or equal to theta_max 
        

    %%% 2... k_ratio
    
        k_ratio = 1;
       
    %%% 3... purine_pyrimidine ratio
        
        pu_py_ratio = 1;
       
     
        
    %%% 4... frequency get the empirical freq of nucleotides to fix the parameter
         
          nt_freq = [0.25 0.25 0.25 0.25];
          
          
           %distances = seqpdist(Data_leaves,'Method','Jukes-Cantor','Alpha','DNA');
        [mat,distances] = fn_seqdistmatrix(Data_leaves,subst_model,nt_freq);
         UPGMA_arrangee =  fn_UPGMA_tree_Structure(distances,no_leaves);
             
         
    
for t = 1:T 
   t 
   T
    
    
    if t == 1
                  
            for d = 1:N
                   %%% 1..... sample theta
                       theta_sample = unifrnd(theta_min,theta_max);
                       THETA_SAMPLES(d) = theta_sample;
                   %%% 2..... sample tree and change its arm to unit of number of substitutions
                       %TreeU = fn_simu_tree(n,theta_sample);
                       TreeU = fn_simu_tree_UPGMA(n,UPGMA_arrangee);
                       TREE_SAMPLES(:,:,d) = TreeU(n+1:end,:);
                  
                   %%% ..... sample substitution model parameters
            end

   
    
    
    else 
        
            ep_change = Epsil(t) - Epsil(t-1);
            ep_t = Epsil(t);

            %%% Computation of the weights for each sample

                wt_lana = wt(t-1,:);
                wt_unnorm = zeros(1,N);

                for a = 1:N
                    
                    %%% 1. Get theta
                        theta_sam = THETA_SAMPLES(a);
                    
                    %%% 2. Get the tree
                    tree_i = [TREE_HEAD; TREE_SAMPLES(:,:,a)];
                    
                   
                    
                    %%% 3. Calculate the Likelihood
                    [lik,log_10_lik] = fn_likelihood_general(subst_model,tree_i,theta_sam,Data,k_ratio,pu_py_ratio,nt_freq);
                
                    wt_lana_i = wt_lana(a);
                    
                    wt_unnorm(a) = (wt_lana_i)*(10)^(log_10_lik*ep_change);
                    %wt_unnorm(a) = wt_lana_i*(lik)^(ep_change);
                    
                    
                    
                    
                end

                
              %%% Normalization of weights 
                wt_loni = wt_unnorm/sum(wt_unnorm);
                wt(t,:) = wt_loni;

              %%% Resampling step

                  neff = 1/sum(wt_loni.^2);

                  if neff < N/10

                      rr = randsample([1:N],N,true,wt_loni);
                      
                     %%% 1 Resample the trees
                           TREE_SAMPLES = TREE_SAMPLES(:,:,rr);
                      
                     %%% 2 Resample theta
                      
                      wt(t,:) = 1/N;
                  end


                  
                  
        %%% MOVE THE SAMPLES

              for m = 1:N
                  
                  %%% Get the jth parameter
                  Tree_init = [TREE_HEAD; TREE_SAMPLES(:,:,m)];
                  theta_j = THETA_SAMPLES(m);

                  
                      %%% MCMC step
                        
                          ALL_TREE_MCMC = zeros(2*n-1,6,N_mcmc + 1);
                             ALL_TREE_MCMC(:,:,1) = Tree_init;
                             
                          ALL_THETA_MCMC = zeros(1,N_mcmc + 1);
                             ALL_THETA_MCMC(1) = theta_j;
                             
                          
                          for i = 2:N_mcmc + 1
                               
                              Tree_lana = ALL_TREE_MCMC(:,:,i-1);
                              theta_lana = ALL_THETA_MCMC(i-1);

                                 for b = 1:B
                                     
                                            %%% FOR THETHA
                                             if b == 1
                                                   
                                                   %%% proposal for theta
                                                    theta_prop = fn_theta_proposal(theta_lana,theta_max,theta_min,delta_theta);
                                                 
                                                   
                                                    
                                                     %%%%% prior ratio 
                                                       prior_ratio = fn_coalesce_prior_ratio(Tree_lana,Tree_lana,n,theta_prop,theta_lana);
                                                 
                                                     %%%%% likelihood ratio 
                                                       [lik_num, log_10_lik_num] = fn_likelihood_general(subst_model,Tree_lana,theta_prop,Data,k_ratio,pu_py_ratio,nt_freq);
                                                       [lik_den, log_10_lik_den] = fn_likelihood_general(subst_model,Tree_lana,theta_lana,Data,k_ratio,pu_py_ratio,nt_freq);
                                                       likelihood_ratio = 10^((log_10_lik_num - log_10_lik_den)*ep_t);
                                                           
                                                       
                                                       
                                                      
                                                  %%% get the acceptance ratio

                                                    A = likelihood_ratio*prior_ratio;

                                                    alpha = min(1,A);

                                                    u = unifrnd(0,1);

                                                    if u <= alpha  % consider the new tree
                                                       ALL_THETA_MCMC(i) = theta_prop; 
                                                    else  % take the old tree
                                                       ALL_THETA_MCMC(i) = theta_lana;   
                                                    end
                                                    
                                             else  %%% FOR the Tree
                                                 
                                                 theta_nibi = ALL_THETA_MCMC(i);
                                                 
                                                   %%% get a new tree from the proposal distribution
                                                   Tree_proposed = fn_BKYF_proposal(Tree_lana,n);
                                                    
                                                                                         
                                                   
                                                   
                                                   %%%%% prior ratio 
                                                       %prior_ratio = fn_coalesce_prior_ratio(Tree_proposed,Tree_lana,n,theta_nibi,theta_nibi);
                                                 
                                                     %%%%% likelihood ratio 
                                                       [lik_num, log_10_lik_num] = fn_likelihood_general(subst_model,Tree_proposed,theta_nibi,Data,k_ratio,pu_py_ratio,nt_freq);
                                                       [lik_den, log_10_lik_den] = fn_likelihood_general(subst_model,Tree_lana,theta_nibi,Data,k_ratio,pu_py_ratio,nt_freq);
                                                       likelihood_ratio = 10^((log_10_lik_num - log_10_lik_den)*ep_t);
                                                       
                                                       
                                                       
                                                      
                                                  %%% get the acceptance ratio
                                                    %A = likelihood_ratio*prior_ratio;
                                                    A = likelihood_ratio;

                                                    alpha = min(1,A);

                                                    u = unifrnd(0,1);

                                                    if u <= alpha  % consider the new tree
                                                       ALL_TREE_MCMC(:,:,i) = Tree_proposed; 
                                                    else  % take the old tree
                                                       ALL_TREE_MCMC(:,:,i) = Tree_lana;   
                                                    end
                                                   
                                                   
                                             end
                                 end  % end for b = 1:B
                             
                          end

                          
                  %%% moved particles
                  TREE_SAMPLES(:,:,m) = ALL_TREE_MCMC(n+1:end,:,end); 
                  THETA_SAMPLES(m) = ALL_THETA_MCMC(end);
              end 
      
      
      
      
      
   end
end


%%% FINAL WEIGHT
wt_final = wt(end,:);

%%% Estimates
    %%% 1.  THETA
      
            mean_theta = sum(wt_final.*THETA_SAMPLES)
            
            
    %%% 2.  Estimates relating to the TREE....e.g time to MRCA
         
         

     
    
    
    

    
case 2    %%% K80   
        
   
    
    
                    
%%%%%%%%%%%% CONSTRUCT THE FIRST ARRANGEMENT OF THE TREE %%%%%%%%%%%%

        
 

        
        %%%%%%%%%%%% MODEL SPECIFIC VARIABLES %%%%%%%%%%%%



%%% Number of parameters to sample
           
       B = 3; 
           %%%% 1. tree
           %%%% 2. Theta
           %%%% 3. k_ratio
          
           
           
%%% STORAGE
            TREE_SAMPLES = zeros(n-1,6,N);
            THETA_SAMPLES = zeros(1,N);
            K_RATIO_SAMPLES = zeros(1,N);
        




%%% PRIORS DISTRIBUTION PARAMETERS
    %%% 1... THETA
        theta_max = 1;
        theta_min = 0.002;
        delta_theta = 0.5;  % less than or equal to theta_max 
        

    %%% 2... k_ratio
    
        k_ratio_max = 5;
        k_ratio_min = 0;
        delta_k_ratio = 3;  % less than or equal to theta_max 
        
        
    %%% 3... purine_pyrimidine ratio
        
        pu_py_ratio = 1;
       
        
        
    %%% 4... frequency get the empirical freq of nucleotides to fix the parameter
         
          nt_freq = [0.25 0.25 0.25 0.25];
        

         %distances = seqpdist(Data_leaves,'Method','Jukes-Cantor','Alpha','DNA');
        [mat,distances] = fn_seqdistmatrix(Data_leaves,subst_model,nt_freq);
         UPGMA_arrangee =  fn_UPGMA_tree_Structure(distances,no_leaves);

                 


    
for t = 1:T %T
   t 
   T
    
    
    
    if t == 1
                  
            for d = 1:N
                 
                   %%% 1..... sample k_ratio 
                       k_ratio_sample = unifrnd(k_ratio_min,k_ratio_max);
                       K_RATIO_SAMPLES(d) = k_ratio_sample;
                   
                   %%% 2..... sample theta
                       theta_sample = unifrnd(theta_min,theta_max);
                       THETA_SAMPLES(d) = theta_sample;
                   %%% 3..... sample tree and change its arm to unit of number of substitutions
                       %TreeU = fn_simu_tree(n,theta_sample);
                       TreeU = fn_simu_tree_UPGMA(n,UPGMA_arrangee);
                       TREE_SAMPLES(:,:,d) = TreeU(n+1:end,:);
                  
                   %%% ..... sample substitution model parameters
            end

   
    
    
    else 
        
            ep_change = Epsil(t) - Epsil(t-1);
            ep_t = Epsil(t);

            %%% Computation of the weights for each sample

                wt_lana = wt(t-1,:);
                wt_unnorm = zeros(1,N);

                for a = 1:N
                    
                    %%% 1. theta
                    theta_sam = THETA_SAMPLES(a);
                    
                    %%% 2. Get the tree
                    tree_i = [TREE_HEAD; TREE_SAMPLES(:,:,a)];
                    
                    %%% 3. Get the subst model parameters
                        k_part = K_RATIO_SAMPLES(a);

                    
                    %%% 4. Calculate the Likelihood
                    [lik,log_10_lik] = fn_likelihood_general(subst_model,tree_i,theta_sam,Data,k_part,pu_py_ratio,nt_freq);
                
                    wt_lana_i = wt_lana(a);
                    
                    wt_unnorm(a) = (wt_lana_i)*(10)^(log_10_lik*ep_change);
                    %wt_unnorm(a) = wt_lana_i*(lik)^(ep_change);
                    
                    
                    
                end

                
              %%% Normalization of weights 
                wt_loni = wt_unnorm/sum(wt_unnorm);
                wt(t,:) = wt_loni;

              %%% Resampling step

                  neff = 1/sum(wt_loni.^2);

                  if neff < N/10

                      rr = randsample([1:N],N,true,wt_loni);
                      
                     %%% 1 Resample the trees
                           TREE_SAMPLES = TREE_SAMPLES(:,:,rr);
                      
                     %%% 2 Resample theta
                      
                      wt(t,:) = 1/N;
                  end


                  
                  
        %%% MOVE THE SAMPLES

              for m = 1:N
                  
                  %%% Get the jth parameter
                  Tree_init = [TREE_HEAD; TREE_SAMPLES(:,:,m)];
                  theta_j = THETA_SAMPLES(m);
                  k_value = K_RATIO_SAMPLES(m);
                  
                      %%% MCMC step
                        
                          ALL_TREE_MCMC = zeros(2*n-1,6,N_mcmc + 1);
                             ALL_TREE_MCMC(:,:,1) = Tree_init;
                             
                          ALL_THETA_MCMC = zeros(1,N_mcmc + 1);
                             ALL_THETA_MCMC(1) = theta_j;
                            
                          ALL_K_RATIO_MCMC = zeros(1,N_mcmc + 1);
                             ALL_K_RATIO_MCMC(1) = k_value;
                          
                          for i = 2:N_mcmc + 1
                               
                              Tree_lana = ALL_TREE_MCMC(:,:,i-1);
                              theta_lana = ALL_THETA_MCMC(i-1);
                              k_ratio_lana = ALL_K_RATIO_MCMC(i-1);
                              
                                 for b = 1:B
                                             if b == 1 %%% FOR k_ratio
                                                 
                                                 
                                                   %%% proposal for k_ratio
                                                    k_ratio_prop = fn_theta_proposal(k_ratio_lana,k_ratio_max,k_ratio_min,delta_k_ratio);
                                                 
                                                  
                                                   %%% get the likelihood ratio
                                                    [lik_num, log_10_lik_num] = fn_likelihood_general(subst_model,Tree_lana,theta_lana,Data,k_ratio_prop,pu_py_ratio,nt_freq);
                                                    [lik_den, log_10_lik_den]  = fn_likelihood_general(subst_model,Tree_lana,theta_lana,Data,k_ratio_lana,pu_py_ratio,nt_freq);
                                                       
                                                   likelihood_ratio = 10^((log_10_lik_num - log_10_lik_den)*ep_t);
                                                    
                                                        
                                                      
                                                  %%% get the acceptance ratio

                                                    A = likelihood_ratio;
                                                    alpha = min(1,A);

                                                    u = unifrnd(0,1);

                                                    if u <= alpha  % consider the new tree
                                                       ALL_K_RATIO_MCMC(i) = k_ratio_prop; 
                                                    else  % take the old tree
                                                       ALL_K_RATIO_MCMC(i) = k_ratio_lana;   
                                                    end
                                                 
                                                 
                                             elseif b == 2  %%% FOR THETHA
                                                 
                                                 k_ratio_nibi = ALL_K_RATIO_MCMC(i);
                                                   
                                                   %%% proposal for theta
                                                    theta_prop = fn_theta_proposal(theta_lana,theta_max,theta_min,delta_theta);
                                                 
                                                    
                                                   %%%%% prior ratio 
                                                    prior_ratio = fn_coalesce_prior_ratio(Tree_lana,Tree_lana,n,theta_prop,theta_lana);
                                                 
                                                  %%% get the likelihood ratio
                                                   [lik_num, log_10_lik_num] = fn_likelihood_general(subst_model,Tree_lana,theta_prop,Data,k_ratio_nibi,pu_py_ratio,nt_freq);
                                                   [lik_den, log_10_lik_den] = fn_likelihood_general(subst_model,Tree_lana,theta_lana,Data,k_ratio_nibi,pu_py_ratio,nt_freq);
                                                       
                                                   likelihood_ratio = 10^((log_10_lik_num - log_10_lik_den)*ep_t);
                                                   
                                                       
                                                       
                                                    
                                                    A = likelihood_ratio*prior_ratio;

                                                    alpha = min(1,A);

                                                    u = unifrnd(0,1);

                                                    if u <= alpha  % consider the new tree
                                                       ALL_THETA_MCMC(i) = theta_prop; 
                                                    else  % take the old tree
                                                       ALL_THETA_MCMC(i) = theta_lana;   
                                                    end
                                                    
                                             else  %%% FOR the Tree
                                                 
                                                 k_ratio_nibi = ALL_K_RATIO_MCMC(i);
                                                 theta_nibi = ALL_K_RATIO_MCMC(i);
                                                 
                                                   %%% get a new tree from the proposal distribution
                                                   Tree_proposed = fn_BKYF_proposal(Tree_lana,n);
                                                    
                                                   %%%%% prior ratio 
                                                   % prior_ratio = fn_coalesce_prior_ratio(Tree_proposed,Tree_lana,n,theta_nibi,theta_nibi);
                                                 
                             
                                                    %%% get the likelihood ratio
                                                   [lik_num, log_10_lik_num] = fn_likelihood_general(subst_model,Tree_proposed,theta_nibi,Data,k_ratio_nibi,pu_py_ratio,nt_freq);
                                                   [lik_den, log_10_lik_den] = fn_likelihood_general(subst_model,Tree_lana,theta_nibi,Data,k_ratio_nibi,pu_py_ratio,nt_freq);
                                                       
                                                   likelihood_ratio = 10^((log_10_lik_num - log_10_lik_den)*ep_t);
                                                   
                                                            
                                                        
                                                       
                                                       
                                                   %%% get the acceptance ratio
                                                   A = likelihood_ratio;
                                                   alpha = min(1,A);

                                                    u = unifrnd(0,1);

                                                    if u <= alpha  % consider the new tree
                                                       ALL_TREE_MCMC(:,:,i) = Tree_proposed; 
                                                    else  % take the old tree
                                                       ALL_TREE_MCMC(:,:,i) = Tree_lana;   
                                                    end
                                                   
                                                   
                                             end
                                 end  % end for b = 1:B
                             
                          end

                          
                  %%% moved particles
                  TREE_SAMPLES(:,:,m) = ALL_TREE_MCMC(n+1:end,:,end); 
                  THETA_SAMPLES(m) = ALL_THETA_MCMC(end);
              end 
      
      
      
      
      
   end
end


%%% FINAL WEIGHT
wt_final = wt(end,:);

%%% Estimates
    %%% 1.  THETA
      
            mean_theta = sum(wt_final.*THETA_SAMPLES)
            
            
    %%% 2.  K_ratio
         mean_k_ratio = sum(wt_final.*K_RATIO_SAMPLES)
         
     %%% 2.  Estimates relating to the TREE....e.g time to MRCA
    
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
case 3    %%% F84
        
        
    
                   
%%%%%%%%%%%% CONSTRUCT THE FIRST ARRANGEMENT OF THE TREE %%%%%%%%%%%%

        

        
        %%%%%%%%%%%% MODEL SPECIFIC VARIABLES %%%%%%%%%%%%



%%% Number of parameters to sample
           
       B = 4; 
           %%%% 1. tree
           %%%% 2. Theta
           %%%% 3. k_ratio
           %%%% 4. nucleotide freq 
           
           
%%% STORAGE
            TREE_SAMPLES = zeros(n-1,6,N);
            THETA_SAMPLES = zeros(1,N);
            K_RATIO_SAMPLES = zeros(1,N);
            FREQ_SAMPLES = zeros(4,N);




%%% PRIORS DISTRIBUTION PARAMETERS
    %%% 1... THETA
        theta_max = 1;
        theta_min = 0.04;
        delta_theta = 0.8;  % less than or equal to theta_max 
        

    %%% 2... k_ratio
    
        k_ratio_max = 500;
        k_ratio_min = 0;
        delta_k_ratio = 400;  % less than or equal to theta_max 
        
        
            
    %%% 3... purine_pyrimidine ratio
        
        pu_py_ratio = 1;
      
        
    %%% 4... frequency get the empirical freq of nucleotides to fix the parameter
         
             n_A = length(find(Data_Matrix_Int == 1));
             n_C = length(find(Data_Matrix_Int == 2));
             n_G = length(find(Data_Matrix_Int == 3));
             n_T = length(find(Data_Matrix_Int == 4));
             
             
             f_A = n_A/(n_A + n_C + n_G + n_T);
             f_C = n_C/(n_A + n_C + n_G + n_T);
             f_G = n_G/(n_A + n_C + n_G + n_T);
             f_T = n_T/(n_A + n_C + n_G + n_T);
             
          emp_freq = [f_A f_C f_G f_T];
          w_peak1 = 100;
          

       %distances = seqpdist(Data_leaves,'Method','Jukes-Cantor','Alpha','DNA');
        [mat,distances] = fn_seqdistmatrix(Data_leaves,subst_model,emp_freq);
         UPGMA_arrangee =  fn_UPGMA_tree_Structure(distances,no_leaves);
               
    



              
    
for t = 1:T
   t 
    
    
    
    if t == 1
                  
            for d = 1:N
                 
                   %%% 1..... sample k_ratio 
                       k_ratio_sample = unifrnd(k_ratio_min,k_ratio_max);
                       K_RATIO_SAMPLES(d) = k_ratio_sample;
                       
                   %%% 2..... sample theta 
                       freq_sample = drchrnd(w_peak1*emp_freq,1);
                       FREQ_SAMPLES(:,d) = freq_sample;
          
                   
                   %%% 3..... sample theta
                       theta_sample = unifrnd(theta_min,theta_max);
                       THETA_SAMPLES(d) = theta_sample;
                   %%% 4..... sample tree and change its arm to unit of number of substitutions
                       %TreeU = fn_simu_tree(n,theta_sample);
                       TreeU = fn_simu_tree_UPGMA(n,UPGMA_arrangee);
                       TREE_SAMPLES(:,:,d) = TreeU(n+1:end,:);
                  
                   %%% ..... sample substitution model parameters
            end

   
    
    
    else 
        
            ep_change = Epsil(t) - Epsil(t-1);
            ep_t = Epsil(t);

            %%% Computation of the weights for each sample

                wt_lana = wt(t-1,:);
                wt_unnorm = zeros(1,N);

                for a = 1:N
                    
                    %%% 1. Get the theta
                        theta_sam = THETA_SAMPLES(a);
                    
                    %%% 2. Get the subst model parameters
                        k_part = K_RATIO_SAMPLES(a);
                    
                    %%% 3. freq
                        freq_sam = FREQ_SAMPLES(:,a);
                    
                    %%% 4. Get the tree
                    tree_i = [TREE_HEAD; TREE_SAMPLES(:,:,a)];
                    
                    

                    %%% 5. Calculate the Likelihood
                    [lik,log_10_lik] = fn_likelihood_general(subst_model,tree_i,theta_sam,Data,k_part,pu_py_ratio,freq_sam);


                    wt_lana_i = wt_lana(a);
                    
                    wt_unnorm(a) = (wt_lana_i)*(10)^(log_10_lik*ep_change);
                    %wt_unnorm(a) = wt_lana_i*(lik)^(ep_change);
                end

                
              %%% Normalization of weights 
                wt_loni = wt_unnorm/sum(wt_unnorm);
                wt(t,:) = wt_loni;

              %%% Resampling step

                  neff = 1/sum(wt_loni.^2);

                  if neff < N/10

                      rr = randsample([1:N],N,true,wt_loni);
                      
                     %%% 1 Resample the trees
                           TREE_SAMPLES = TREE_SAMPLES(:,:,rr);
                      
                     %%% 2 Resample theta
                      
                      wt(t,:) = 1/N;
                  end


                  
                  
        %%% MOVE THE SAMPLES

              for m = 1:N
                  
                  %%% Get the jth parameter
                  Tree_init = [TREE_HEAD; TREE_SAMPLES(:,:,m)];
                  theta_j = THETA_SAMPLES(m);
                  k_value = K_RATIO_SAMPLES(m);
                  freq_val = FREQ_SAMPLES(:,m);
                  
                      %%% MCMC step
                        
                          ALL_TREE_MCMC = zeros(2*n-1,6,N_mcmc + 1);
                             ALL_TREE_MCMC(:,:,1) = Tree_init;
                             
                          ALL_THETA_MCMC = zeros(1,N_mcmc + 1);
                             ALL_THETA_MCMC(1) = theta_j;
                            
                          ALL_K_RATIO_MCMC = zeros(1,N_mcmc + 1);
                             ALL_K_RATIO_MCMC(1) = k_value;
                             
                          ALL_FREQ_MCMC = zeros(4,N_mcmc + 1);  
                             ALL_FREQ_MCMC(:,1) = freq_val; 
                          
                          for i = 2:N_mcmc + 1
                               
                              Tree_lana = ALL_TREE_MCMC(:,:,i-1);
                              theta_lana = ALL_THETA_MCMC(i-1);
                              k_ratio_lana = ALL_K_RATIO_MCMC(i-1);
                              freq_lana = ALL_FREQ_MCMC(:,i-1);
                                 for b = 1:B
                                             if b == 1 %%% FOR k_ratio
                                                 
                                                 
                                                   %%% proposal for k_ratio
                                                    k_ratio_prop = fn_theta_proposal(k_ratio_lana,k_ratio_max,k_ratio_min,delta_k_ratio);
                                                 
                                                   
                                                    %%% get the likelihood ratio
                                                   [lik_num, log_10_lik_num] = fn_likelihood_general(subst_model,Tree_lana,theta_lana,Data,k_ratio_prop,pu_py_ratio,freq_lana);
                                                   [lik_den, log_10_lik_den] = fn_likelihood_general(subst_model,Tree_lana,theta_lana,Data,k_ratio_lana,pu_py_ratio,freq_lana);
                                                      
                                                   %%% get the acceptance ratio

                                                    A = 10^((log_10_lik_num - log_10_lik_den)*ep_t);

                                                    alpha = min(1,A);

                                                    u = unifrnd(0,1);

                                                    if u <= alpha  % consider the new tree
                                                       ALL_K_RATIO_MCMC(i) = k_ratio_prop; 
                                                    else  % take the old tree
                                                       ALL_K_RATIO_MCMC(i) = k_ratio_lana;   
                                                    end
                                                 
                                               
                                             elseif b == 2  %%% FOR FREQUENCY
                                                   
                                                     %%% proposal for k_ratio
                                                    freq_prop = drchrnd(w_peak1*freq_lana',1);
                                                 
                                                   
                                                    k_ratio_nibi = ALL_K_RATIO_MCMC(i); 
                                                    
                                                    %%% prior ratio
                                                    prior_num = dirpdftunji(freq_prop,w_peak1*emp_freq);
                                                    prior_den = dirpdftunji(freq_lana',w_peak1*emp_freq);
                                                      prior_ratio = prior_num/prior_den;
                                                     
                                                    
                                                    %%% get the likelihood ratio
                                                   [lik_num, log_10_lik_num] = fn_likelihood_general(subst_model,Tree_lana,theta_lana,Data,k_ratio_nibi,pu_py_ratio,freq_prop);
                                                   [lik_den, log_10_lik_den] = fn_likelihood_general(subst_model,Tree_lana,theta_lana,Data,k_ratio_nibi,pu_py_ratio,freq_lana);
                                                       
                                                   likelihood_ratio = 10^((log_10_lik_num - log_10_lik_den)*ep_t);
                                                     
                                                     %%% proposal ratio
                                                      prop_num = dirpdftunji(freq_prop,w_peak1*freq_lana');
                                                      prop_den = dirpdftunji(freq_lana',w_peak1*freq_prop);
                                                         prop_ratio = prop_num/prop_den; 
                                                      
                                                  %%% get the acceptance ratio

                                                    A = likelihood_ratio*prior_ratio*prop_ratio;

                                                    alpha = min(1,A);

                                                    u = unifrnd(0,1);

                                                    if u <= alpha  % consider the new tree
                                                       ALL_FREQ_MCMC(:,i) = freq_prop'; 
                                                    else  % take the old tree
                                                       ALL_FREQ_MCMC(:,i) = freq_lana;   
                                                    end
                                                 
                                        
                                                 
                                                 
                                             elseif b == 3  %%% FOR THETHA
                                                 
                                                 k_ratio_nibi = ALL_K_RATIO_MCMC(i);
                                                 freq_nibi = ALL_FREQ_MCMC(:,i);
                                                   
                                                   %%% proposal for theta
                                                    theta_prop = fn_theta_proposal(theta_lana,theta_max,theta_min,delta_theta);
                                               
                                                  
                                                  %%% get the prior ratio 
                                                   prior_ratio = fn_coalesce_prior_ratio(Tree_lana,Tree_lana,n,theta_prop,theta_lana);
                                                 %%% get the likelihood ratio
                                                   [lik_num, log_10_lik_num] = fn_likelihood_general(subst_model,Tree_lana,theta_prop,Data,k_ratio_nibi,pu_py_ratio,freq_nibi);
                                                  [lik_den, log_10_lik_den] = fn_likelihood_general(subst_model,Tree_lana,theta_lana,Data,k_ratio_nibi,pu_py_ratio,freq_nibi);
                                                     lik_ratio = 10^((log_10_lik_num - log_10_lik_den)*ep_t);   
                                                   
                                                     
                                                    
                                                   
                                                    A = lik_ratio*prior_ratio; 

                                                    alpha = min(1,A);

                                                    u = unifrnd(0,1);

                                                    if u <= alpha  % consider the new tree
                                                       ALL_THETA_MCMC(i) = theta_prop; 
                                                    else  % take the old tree
                                                       ALL_THETA_MCMC(i) = theta_lana;   
                                                    end
                                                    
                                             else  %%% FOR the Tree
                                                 
                                                 k_ratio_nibi = ALL_K_RATIO_MCMC(i);
                                                 theta_nibi = ALL_K_RATIO_MCMC(i);
                                                 freq_nibi = ALL_FREQ_MCMC(:,i);
                                                 
                                                   %%% get a new tree from the proposal distribution
                                                   Tree_proposed = fn_BKYF_proposal(Tree_lana,n);
                                                    
                                                   
                                                   %%% get the prior ratio 
                                                   %prior_ratio = fn_coalesce_prior_ratio(Tree_proposed,Tree_lana,n,theta_nibi,theta_nibi);
                                                 %%% get the likelihood ratio
                                                   [lik_num, log_10_lik_num] = fn_likelihood_general(subst_model,Tree_proposed,theta_nibi,Data,k_ratio_nibi,pu_py_ratio,freq_nibi);
                                                   [lik_den, log_10_lik_den] = fn_likelihood_general(subst_model,Tree_lana,theta_nibi,Data,k_ratio_nibi,pu_py_ratio,freq_nibi);
                                                     lik_ratio = 10^((log_10_lik_num - log_10_lik_den)*ep_t);   
                                                   
                                                     
                                                     
                                                   
                                                    A = lik_ratio; 

                                                    alpha = min(1,A);

                                                    u = unifrnd(0,1);

                                                    if u <= alpha  % consider the new tree
                                                       ALL_TREE_MCMC(:,:,i) = Tree_proposed; 
                                                    else  % take the old tree
                                                       ALL_TREE_MCMC(:,:,i) = Tree_lana;   
                                                    end
                                                   
                                                   
                                             end
                                 end  % end for b = 1:B
                             
                          end

                          
                  %%% moved particles
                  TREE_SAMPLES(:,:,m) = ALL_TREE_MCMC(n+1:end,:,end); 
                  THETA_SAMPLES(m) = ALL_THETA_MCMC(end);
                  FREQ_SAMPLES(:,m) = ALL_FREQ_MCMC(:,end);
              end 
      
      
      
      
      
   end
end


%%% FINAL WEIGHT
wt_final = wt(end,:);

%%% Estimates
    %%% 1.  THETA
      
            mean_theta = sum(wt_final.*THETA_SAMPLES);
            var_theta = sum(wt_final.*((THETA_SAMPLES - mean_theta).^2));
            
    
    
     
     
     
     
     
     
     
     
        
        
        
        
        
    case 4    %%% HKY85  
        
     
                   
%%%%%%%%%%%% CONSTRUCT THE FIRST ARRANGEMENT OF THE TREE %%%%%%%%%%%%

       
               
 

        
        %%%%%%%%%%%% MODEL SPECIFIC VARIABLES %%%%%%%%%%%%



%%% Number of parameters to sample
           
       B = 4; 
           %%%% 1. tree
           %%%% 2. Theta
           %%%% 3. k_ratio
           %%%% 4. nucleotide freq 
           
           
%%% STORAGE
            TREE_SAMPLES = zeros(n-1,6,N);
            THETA_SAMPLES = zeros(1,N);
            K_RATIO_SAMPLES = zeros(1,N);
            FREQ_SAMPLES = zeros(4,N);




%%% PRIORS DISTRIBUTION PARAMETERS
    %%% 1... THETA
        theta_max = 1;
        theta_min = 0.002;
        delta_theta = 0.5;  % less than or equal to theta_max 
        

    %%% 2... k_ratio
    
        k_ratio_max = 5;
        k_ratio_min = 0;
        delta_k_ratio = 3;  % less than or equal to theta_max 
        
        
            
    %%% 3... purine_pyrimidine ratio
        
        pu_py_ratio = 1;
      
        
    %%% 4... frequency get the empirical freq of nucleotides to fix the parameter
         
             n_A = length(find(Data_Matrix_Int == 1));
             n_C = length(find(Data_Matrix_Int == 2));
             n_G = length(find(Data_Matrix_Int == 3));
             n_T = length(find(Data_Matrix_Int == 4));
             
             
             f_A = n_A/(n_A + n_C + n_G + n_T);
             f_C = n_C/(n_A + n_C + n_G + n_T);
             f_G = n_G/(n_A + n_C + n_G + n_T);
             f_T = n_T/(n_A + n_C + n_G + n_T);
             
          emp_freq = [f_A f_C f_G f_T];
          w_peak1 = 10;
          

          
      %distances = seqpdist(Data_leaves,'Method','Jukes-Cantor','Alpha','DNA');
        [mat,distances] = fn_seqdistmatrix(Data_leaves,subst_model,emp_freq);
         UPGMA_arrangee =  fn_UPGMA_tree_Structure(distances,no_leaves);



              
    
for t = 1:2%T %T
   t 
    
    
    
    if t == 1
                  
            for d = 1:N
                 
                   %%% 1..... sample k_ratio 
                       k_ratio_sample = unifrnd(k_ratio_min,k_ratio_max);
                       K_RATIO_SAMPLES(d) = k_ratio_sample;
                       
                   %%% 2..... sample theta 
                       freq_sample = drchrnd(w_peak1*emp_freq,1);
                       FREQ_SAMPLES(:,d) = freq_sample;
          
                   
                   %%% 3..... sample theta
                       theta_sample = unifrnd(theta_min,theta_max);
                       THETA_SAMPLES(d) = theta_sample;
                   %%% 4..... sample tree and change its arm to unit of number of substitutions
                       %TreeU = fn_simu_tree(n,theta_sample);
                       TreeU = fn_simu_tree_UPGMA(n,UPGMA_arrangee);
                       TREE_SAMPLES(:,:,d) = TreeU(n+1:end,:);
                  
                   %%% ..... sample substitution model parameters
            end

   
    
    
    else 
        
            ep_change = Epsil(t) - Epsil(t-1);
            ep_t = Epsil(t);

            %%% Computation of the weights for each sample

                wt_lana = wt(t-1,:);
                wt_unnorm = zeros(1,N);

                for a = 1:N
                    
                    %%% 1. Get the theta
                        theta_sam = THETA_SAMPLES(a);
                    
                    %%% 2. Get the subst model parameters
                        k_part = K_RATIO_SAMPLES(a);
                    
                    %%% 3. freq
                        freq_sam = FREQ_SAMPLES(:,a);
                    
                    %%% 4. Get the tree
                    tree_i = [TREE_HEAD; TREE_SAMPLES(:,:,a)];
                    
                    

                    %%% 5. Calculate the Likelihood
                    [lik,log_10_lik] = fn_likelihood_general(subst_model,tree_i,theta_sam,Data,k_part,pu_py_ratio,freq_sam);


                    wt_lana_i = wt_lana(a);
                    
                    wt_unnorm(a) = (wt_lana_i)*(10)^(log_10_lik*ep_change);
                    %wt_unnorm(a) = wt_lana_i*(lik)^(ep_change);
                end

                
              %%% Normalization of weights 
                wt_loni = wt_unnorm/sum(wt_unnorm);
                wt(t,:) = wt_loni;

              %%% Resampling step

                  neff = 1/sum(wt_loni.^2);

                  if neff < N/10

                      rr = randsample([1:N],N,true,wt_loni);
                      
                     %%% 1 Resample the trees
                           TREE_SAMPLES = TREE_SAMPLES(:,:,rr);
                      
                     %%% 2 Resample theta
                      
                      wt(t,:) = 1/N;
                  end


                  
                  
        %%% MOVE THE SAMPLES

              for m = 1:N
                  
                  %%% Get the jth parameter
                  Tree_init = [TREE_HEAD; TREE_SAMPLES(:,:,m)];
                  theta_j = THETA_SAMPLES(m);
                  k_value = K_RATIO_SAMPLES(m);
                  freq_val = FREQ_SAMPLES(:,m);
                  
                      %%% MCMC step
                        
                          ALL_TREE_MCMC = zeros(2*n-1,6,N_mcmc + 1);
                             ALL_TREE_MCMC(:,:,1) = Tree_init;
                             
                          ALL_THETA_MCMC = zeros(1,N_mcmc + 1);
                             ALL_THETA_MCMC(1) = theta_j;
                            
                          ALL_K_RATIO_MCMC = zeros(1,N_mcmc + 1);
                             ALL_K_RATIO_MCMC(1) = k_value;
                             
                          ALL_FREQ_MCMC = zeros(4,N_mcmc + 1);  
                             ALL_FREQ_MCMC(:,1) = freq_val; 
                          
                          for i = 2:N_mcmc + 1
                               
                              Tree_lana = ALL_TREE_MCMC(:,:,i-1);
                              theta_lana = ALL_THETA_MCMC(i-1);
                              k_ratio_lana = ALL_K_RATIO_MCMC(i-1);
                              freq_lana = ALL_FREQ_MCMC(:,i-1);
                                 for b = 1:B
                                             if b == 1 %%% FOR k_ratio
                                                 
                                                 
                                                   %%% proposal for k_ratio
                                                    k_ratio_prop = fn_theta_proposal(k_ratio_lana,k_ratio_max,theta_min,delta_k_ratio);
                                                 
                                                   
                                                    %%% get the likelihood ratio
                                                   [lik_num, log_10_lik_num] = fn_likelihood_general(subst_model,Tree_lana,theta_lana,Data,k_ratio_prop,pu_py_ratio,freq_lana);
                                                   [lik_den, log_10_lik_den] = fn_likelihood_general(subst_model,Tree_lana,theta_lana,Data,k_ratio_lana,pu_py_ratio,freq_lana);
                                                      
                                                   %%% get the acceptance ratio

                                                    A = 10^((log_10_lik_num - log_10_lik_den)*ep_t);

                                                    alpha = min(1,A);

                                                    u = unifrnd(0,1);

                                                    if u <= alpha  % consider the new tree
                                                       ALL_K_RATIO_MCMC(i) = k_ratio_prop; 
                                                    else  % take the old tree
                                                       ALL_K_RATIO_MCMC(i) = k_ratio_lana;   
                                                    end
                                                 
                                               
                                             elseif b == 2  %%% FOR FREQUENCY
                                                   
                                                     %%% proposal for k_ratio
                                                    freq_prop = drchrnd(w_peak1*freq_lana',1);
                                                 
                                                   
                                                    k_ratio_nibi = ALL_K_RATIO_MCMC(i); 
                                                    
                                                    %%% prior ratio
                                                    prior_num = dirpdftunji(freq_prop,w_peak1*emp_freq);
                                                    prior_den = dirpdftunji(freq_lana',w_peak1*emp_freq);
                                                      prior_ratio = prior_num/prior_den;
                                                     
                                                    
                                                    %%% get the likelihood ratio
                                                   [lik_num, log_10_lik_num] = fn_likelihood_general(subst_model,Tree_lana,theta_lana,Data,k_ratio_nibi,pu_py_ratio,freq_prop);
                                                   [lik_den, log_10_lik_den] = fn_likelihood_general(subst_model,Tree_lana,theta_lana,Data,k_ratio_nibi,pu_py_ratio,freq_lana);
                                                       
                                                   likelihood_ratio = 10^((log_10_lik_num - log_10_lik_den)*ep_t);
                                                     
                                                     %%% proposal ratio
                                                      prop_num = dirpdftunji(freq_prop,w_peak1*freq_lana');
                                                      prop_den = dirpdftunji(freq_lana',w_peak1*freq_prop);
                                                         prop_ratio = prop_num/prop_den; 
                                                      
                                                  %%% get the acceptance ratio

                                                    A = likelihood_ratio*prior_ratio*prop_ratio;

                                                    alpha = min(1,A);

                                                    u = unifrnd(0,1);

                                                    if u <= alpha  % consider the new tree
                                                       ALL_FREQ_MCMC(:,i) = freq_prop'; 
                                                    else  % take the old tree
                                                       ALL_FREQ_MCMC(:,i) = freq_lana;   
                                                    end
                                                 
                                        
                                                 
                                                 
                                             elseif b == 3  %%% FOR THETHA
                                                 
                                                 k_ratio_nibi = ALL_K_RATIO_MCMC(i);
                                                 freq_nibi = ALL_FREQ_MCMC(:,i);
                                                   
                                                   %%% proposal for theta
                                                    theta_prop = fn_theta_proposal(theta_lana,theta_max,theta_min,delta_theta);
                                               
                                                  
                                                  %%% get the prior ratio 
                                                   prior_ratio = fn_coalesce_prior_ratio(Tree_lana,Tree_lana,n,theta_prop,theta_lana);
                                                 %%% get the likelihood ratio
                                                   [lik_num, log_10_lik_num] = fn_likelihood_general(subst_model,Tree_lana,theta_prop,Data,k_ratio_nibi,pu_py_ratio,freq_nibi);
                                                  [lik_den, log_10_lik_den] = fn_likelihood_general(subst_model,Tree_lana,theta_lana,Data,k_ratio_nibi,pu_py_ratio,freq_nibi);
                                                     lik_ratio = 10^((log_10_lik_num - log_10_lik_den)*ep_t);   
                                                   
                                                     
                                                    
                                                   
                                                    A = lik_ratio*prior_ratio; 

                                                    alpha = min(1,A);

                                                    u = unifrnd(0,1);

                                                    if u <= alpha  % consider the new tree
                                                       ALL_THETA_MCMC(i) = theta_prop; 
                                                    else  % take the old tree
                                                       ALL_THETA_MCMC(i) = theta_lana;   
                                                    end
                                                    
                                             else  %%% FOR the Tree
                                                 
                                                 k_ratio_nibi = ALL_K_RATIO_MCMC(i);
                                                 theta_nibi = ALL_K_RATIO_MCMC(i);
                                                 freq_nibi = ALL_FREQ_MCMC(:,i);
                                                 
                                                   %%% get a new tree from the proposal distribution
                                                   Tree_proposed = fn_BKYF_proposal(Tree_lana,n);
                                                    
                                                   
                                                   %%% get the prior ratio 
                                                   %prior_ratio = fn_coalesce_prior_ratio(Tree_proposed,Tree_lana,n,theta_nibi,theta_nibi);
                                                 %%% get the likelihood ratio
                                                   [lik_num, log_10_lik_num] = fn_likelihood_general(subst_model,Tree_proposed,theta_nibi,Data,k_ratio_nibi,pu_py_ratio,freq_nibi);
                                                   [lik_den, log_10_lik_den] = fn_likelihood_general(subst_model,Tree_lana,theta_nibi,Data,k_ratio_nibi,pu_py_ratio,freq_nibi);
                                                     lik_ratio = 10^((log_10_lik_num - log_10_lik_den)*ep_t);   
                                                   
                                                     
                                                     
                                                   
                                                    A = lik_ratio; 

                                                    alpha = min(1,A);

                                                    u = unifrnd(0,1);

                                                    if u <= alpha  % consider the new tree
                                                       ALL_TREE_MCMC(:,:,i) = Tree_proposed; 
                                                    else  % take the old tree
                                                       ALL_TREE_MCMC(:,:,i) = Tree_lana;   
                                                    end
                                                   
                                                   
                                             end
                                 end  % end for b = 1:B
                             
                          end

                          
                  %%% moved particles
                  TREE_SAMPLES(:,:,m) = ALL_TREE_MCMC(n+1:end,:,end); 
                  THETA_SAMPLES(m) = ALL_THETA_MCMC(end);
              end 
      
      
      
      
      
   end
end


%%% FINAL WEIGHT
wt_final = wt(end,:);

%%% Estimates
    %%% 1.  THETA
      
            mean_theta = sum(wt_final.*THETA_SAMPLES)
            
            
    %%% 2.  K_ratio
         mean_k_ratio = sum(wt_final.*K_RATIO_SAMPLES)
         
     %%% 2.  Estimates relating to the TREE....e.g time to MRCA

        
        
    
     
     
     
     
     
     
        
        
end

        
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            


    