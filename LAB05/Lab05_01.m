    clc; clear; close all;
    p = 0.5;                                                
    % theoretical mean and variance: E[X]=p, Var[X]=p*(1-p)
    E_theory   = p;
    Var_theory = p*(1-p);
    fprintf('Theoretical: E[X]=%.4f, Var[X]=%.4f\n', E_theory, Var_theory);
    
    % simulate for Nsim = 10, 100, 1000
    Nsim_vals = [10, 100, 1000];
    
    for n = 1:length(Nsim_vals)
        Nsim = Nsim_vals(n);
        X    = double(rand(1,Nsim) < p);   
        fprintf('Nsim=%4d: sim E[X]=%.4f, sim Var[X]=%.4f\n', Nsim, mean(X), var(X));
    end