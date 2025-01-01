function [feasibleFin,sigma_v,rateFin,distortionFin,perceptionFin,classificationFin]=inner_point_RDPCO_m1(Ek,Dk,c,m,n,p_h0,p_h1,D0,P0,epsilon0,sigma_v,iter_lambda)


alpha_0 = 0.00005;
cc = 1;
iter=20000;
% iter=20;
epsilon = 1e-5;
t=0.01;
u=2;
condition =3;

TernimateCondition=0.05;

feasibleFin=0;
rateFin=0;
distortionFin=0;
perceptionFin=0;
classificationFin=0;
weigh_D = 1/log(max(D0));
weigh_P = 1/log(max(P0));
% weigh_C = 1;
weigh_C = -1/log(sqrt(p_h0*p_h1));
gamma_k = Dk*Ek*Ek'*Dk'-eye(n);

for l=1:iter_lambda
    alpha = alpha_0;
    %     if l>=5
    %         alpha = 0.5*alpha_0;
    %     end
    %     if l>=10
    %         alpha = 0.1*alpha_0;
    %     end
    %     if l>=15
    %         alpha = 0.05*alpha_0;
    %     end
    if l>=5
        alpha = 0.5*alpha_0;
        if l>=10
            alpha = 0.1*alpha_0;
        end
    end

    for h=1:iter

        distortion_d = -diag(Dk'*Dk);
        perception_d = -(diag(Dk'*Dk)-diag(Dk'*generalize_inv(sigma_v,Dk,Ek)*Dk));
        classification_d = diag(Dk'*generalize_inv(sigma_v,Dk,Ek)*Dk*Ek*c*c'*Ek'*Dk'*generalize_inv(sigma_v,Dk,Ek)*Dk);
        rate_d = -1./(sigma_v.*(sigma_v+1));

        distortion = trace(Dk*Ek*Ek'*Dk'+Dk*diag(sigma_v)*Dk')+0.5*c'*Ek'*Dk'*Dk*Ek*c-trace(Ek*(2*eye(n)+c*c')*Dk)+n+0.5*norm(c)^2;
        H = Dk*Ek*Ek'*Dk'+Dk*diag(sigma_v)*Dk';
        [Q,L] = eig(H);
        H = Q*sqrt(L)*Q';
        perception = p_h1*norm(Dk*Ek*c-c)^2+norm(H-eye(n),"fro")^2;
        classification = c'*Ek'*Dk'*generalize_inv(sigma_v,Dk,Ek)*Dk*Ek*c+8*log(epsilon0/(p_h1*p_h0)^(1/2));
        rate = sum(log(1+1./sigma_v));

        in_p_r =  rate_d *t;
        in_p_d = -weigh_D*distortion_d/(D0-distortion);
        in_p_p = -weigh_P*perception_d/(P0-perception);
        in_p_c = -weigh_C*classification_d/classification;

        sigma_v1 = sigma_v - alpha*(in_p_r+in_p_d+in_p_p+in_p_c);
        % sigma_v1


        f2(cc) = t*rate -weigh_D*log(-(distortion-D0))-weigh_P*log(-(perception-P0))-weigh_C*log(-(classification-epsilon0));
        classification_err = (p_h0*p_h1)^(1/2)*exp(-1/8*c'*Ek'*Dk'*generalize_inv(sigma_v,Dk,Ek)*Dk*Ek*c);

        if sqrt(norm(in_p_r+in_p_d+in_p_p+in_p_c,"fro")) <= epsilon
            break;
        else
            sigma_v = max(sigma_v1,0.0001);
        end

        cc = cc+1;

        if h==15000
            alpha=0.1*alpha;
        end

    end
    %     iteration=1:cc-1;
    %     semilogy(iteration,f2,'-');
    % ;
%     l
%     C1=(distortion-D0)
%     C2=(perception-P0)
%     C3=(classification_err-epsilon0)
    if (distortion-D0)/D0<=TernimateCondition &&...
            (perception-P0)/P0<=TernimateCondition ...
            && (classification_err-epsilon0)/epsilon0<=TernimateCondition
        %可行的话
        %         Inform=['l:'  num2str(l)  ' Feasible Rate:'  num2str(rate) '  distortion:'  num2str(distortion)...
        %             '  perception: ' num2str(perception)  '  classification: ' num2str(classification_err)  'D0: '  num2str(D0) 'P0: '  num2str(P0) 'c0: '  num2str(epsilon0)]
        feasibleFin=1;
        rateFin=rate;
        distortionFin=distortion;
        perceptionFin=perception;
        classificationFin=classification_err;
    else
        %不可行的话
        %         Inform=['l: '  num2str(l)  ' Infeasible  Rate:'  num2str(rate) '  distortion:'  num2str(distortion)...
        %             '  perception:' num2str(perception)  '  classification:' num2str(classification_err)  'D0: '  num2str(D0) 'P0: '  num2str(P0) 'c0: '  num2str(epsilon0)]
        feasibleFin=0;
        break;
    end
    %     ;

    if condition/t<0.01
        feasibleFin=1;
        rateFin=rate;
        distortionFin=distortion;
        perceptionFin=perception;
        classificationFin=classification_err;
        break;
    else
        t = u*t;
    end

end


% if (distortion-D0)<=TernimateCondition && (perception-P0)<=TernimateCondition && (classification-epsilon0)<=TernimateCondition
%     feasible=1;
% else
%     feasible=0;
% end
