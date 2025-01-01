close all;
clear all;
clc;
n = 5; % Source dimension
m = 2; % Compressed dimension or I/Q dimension

D1 = 4.5;
P0 = 4.1;
epsilon0 = 0.2;
p_h0=0.5;
p_h1=0.5;
count_max =50;% control the optimazation iteration
iter_lambda=50;
i_max = 20;%length of D0


count=0;
sigma_v = rand(m,1);

feasible_num =0;
feasible_max =100;
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
        inner_point_RDPCO_m1(Ek2,Dk2,c,m,n,p_h0,p_h1,D1,P0,epsilon0,sigma_v0,1);

    if feasible_DE == 1

        feasible_num = feasible_num +1
    
        for i =1:i_max

            D0 = D1 + 0.5*(i-1);
            Dk = Dk2;
            Ek = Ek2;
            sigma_v = rand(m,1);

            [feasible, sigma_v1, rate_0, distortion_0, perception_0, classification_0] = ...
                    inner_point_RDPCO_m1(Ek,Dk,c,m,n,p_h0,p_h1,D0,P0,epsilon0,sigma_v,iter_lambda);
    
            A(1,i,feasible_num)=rate_0;
            A(2,i,feasible_num)=distortion_0;
            A(3,i,feasible_num)=perception_0;
            A(4,i,feasible_num)=classification_0;
                        

        end

    end 

    if feasible_num == feasible_max
         break;
    end
end

toc

SaveFileName=['m_' num2str(m) 'Data_20231213_Vary_D_' 'D0_' num2str(D1) '_P0_' num2str(P0) ...
    '_epsilon0_' num2str(epsilon0)  '.mat']

save(SaveFileName)


% rate = rate_sum./num_valid;
% distortion = distortion_sum./num_valid;
% perception = perception_sum./num_valid;
% classification = classification_sum./num_valid;
