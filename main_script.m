
clear all; close all;  clc;
rng('default'); rng(3000); %300



%%% Load the input data
load('TDS_n_matrix')
load('TDS_NN_matrix')    
variant_count = TDS_n_matrix;
total_count = TDS_NN_matrix;
alpha = 2; rep = 20;

[matrix_Z_est, matrix_W_est,p_o_est] = fn_haplo_inf(TDS_n_matrix,TDS_NN_matrix,rep,alpha);

disp('This is the estimated matrix of haplotypes:')
disp(matrix_Z_est)
disp('This is the estimated matrix of proportions:')
disp(matrix_W_est)
disp('This is the estimated p_o:')
disp(p_o_est)