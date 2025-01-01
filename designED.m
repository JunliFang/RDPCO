function [Dk1,Ek1]=designED(c,m,n,Sigma_hat,sigma_v,Dk,Ek)

    
    S = eye(n)-Sigma_hat;

    Dk1 = Dk;
    Ek1 = Ek;
    alpha = 0.0001;
    lambda_d1 = 0.5;
    lambda_p1 = 0.5;
    iter=20000; 
    epsilon = 1e-5;
    
    for i=1:iter
    
        A = (S+S')*Dk*diag(sigma_v) + 2*Dk*diag(sigma_v)*Dk'*Dk*diag(sigma_v);
        F = 2*Dk*diag(sigma_v)*Dk'*Dk*diag(sigma_v);
        BD = (Dk*Ek*c-c)*c'*Ek';
        BE = Dk'*(Dk*Ek*c-c)*c';
    
        Dk1 = Dk - alpha*(lambda_d1*A+lambda_p1*F+lambda_d1*BD);
        Ek1 = Ek - alpha*(lambda_p1*BE)-Ek/n*ones(n,1)*ones(1,n);
        f1(i) = lambda_d1*0.5*norm(S+Dk1*diag(sigma_v)*Dk1',"fro")^2+lambda_p1*0.5*norm(c-Dk1*Ek1*c)^2+lambda_p1*0.5*norm(Dk1*diag(sigma_v)*Dk1',"fro")^2;
    
        if sqrt(norm(A+F+BD,"fro")^2+norm(BE,"fro")^2)<epsilon
            break;
        else
            Dk = Dk1;
            Ek = Ek1;
        end
        if i>10000
            alpha=0.1*alpha;
        end
    end
% figure();
% iteration=1:iter;
% semilogy(iteration,f1,'-');