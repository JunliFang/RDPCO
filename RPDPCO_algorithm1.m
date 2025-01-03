close all;
clear all;
clc;
n = 5; % Source dimension
m = 2; % Compressed dimension or I/Q dimension

D0 = 10;
P0 = 10;
epsilon0 = 0.3;
p_h0=0.5;
p_h1=0.5;
iter_lambda=50;

count=0;
sigma_v = rand(m,1);

feasible_num =0;
feasible_max =1;
%the maximum of feasible points 
tic

while(1)
    c = 2*randn(n,1);
    Dk0 = randn(n,m);
    Ek0 = randn(m,n);
    V = [c/norm(c), randn(n, n-1)];
    Q = gramschmidt(V); 
    Sigma_hat = Q*diag([2; rand(m-1,1)+ones(m-1,1); zeros(n-m, 1)])*Q';
    sigma0 = rand(m,1);
    
    [Dk2,Ek2] = designED(c,m,n,Sigma_hat,sigma0,Dk0,Ek0);
    
    sigma_v0 = rand(m,1);
    [feasible_DE, sigma_v1, rate_0, distortion_0, perception_0, classification_0] = ...
        inner_point_RDPCO_m1(Ek2,Dk2,c,m,n,p_h0,p_h1,D0,P0,epsilon0,sigma_v0,1);
    
    if feasible_DE == 1

        feasible_num = feasible_num +1
    
        Dk = Dk2;
        Ek = Ek2;
        sigma_v = rand(m,1);
        
        [feasible, sigma_v1, rate_0, distortion_0, perception_0, classification_0] = ...
                inner_point_RDPCO_m1(Ek,Dk,c,m,n,p_h0,p_h1,D0,P0,epsilon0,sigma_v,iter_lambda);
    
    end

    if feasible_num == feasible_max
         break;
    end
end       

