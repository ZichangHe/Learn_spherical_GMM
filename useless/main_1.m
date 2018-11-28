clear,clc
%% generate GMM variables
%mu = [1 2 3 4;-3 -1 2 2 ;2 1 3 3 ];
k=2; % num of component
d=3; % dimension of variables
mu = rand(k,d);
sigma = [1 1 1]; % shared diagonal covariance matrix
weight=[0.3 0.7];
gm = gmdistribution(mu,sigma,weight)
mu=mu'; %to fit [mu_1 mu_2 ... mu_n]
rng('default'); % For reproducibility
[Y, compIdx] = random(gm,1e+5);
% %% verfication
M2_t=cov(Y); % theory
M2_m=zeros(size(Y,2)); 
for i=1:size(Y,1)
    M2_m=M2_m+Y(i,:)'*Y(i,:);
end
M2_m=M2_m/size(Y,1);
M2_m=M2_m-mean(Y)'*mean(Y); % M2_t = M2_m
%% generate GMM moments: 2-order
M1=mean(Y);
sigma2 = min(eig(cov(Y)));
M2=M2_m-sigma2*eye(size(Y,2));
%M22=zeros(d);
% for n=1:k
%     M22=M22+weight(n)*mu(:,n)*mu(:,n)'; %% M22 £¡= M2
% end
[Uw, Lw, Vw] = svd(M2); 
W = Uw(:,1:k)* sqrt(pinv(Lw(1:k,1:k)));  %% check W'*M2*W = I
%% generate GMM moments: 3-order (whiten Y)
Y=Y*W; %% Whitening!
% M3_a=zeros(size(Y,2),size(Y,2),size(Y,2));
% for i=1:size(Y,1)
%     a=tensor(Y(i,:)'*Y(i,:));
%     b=tensor(Y(i,:));
%     temp=double(ttt(a,b));
%     temp=squeeze(temp); 
%     M3_a=M3_a+temp;
% end
% M3_a=M3_a/size(Y,1);
% M3_b=zeros(size(Y,2),size(Y,2),size(Y,2));
% for i=1:size(Y,2)
%     e=zeros(size(Y,2),1);
%     e(i)=1;
%     %
%     a=tensor(mean(Y)'*e');
%     b=tensor(e);
%     temp1=double(ttt(a,b));
%     temp1=squeeze(temp1);
%     %
%     a=tensor(e*mean(Y));
%     b=tensor(e);
%     temp2=double(ttt(a,b));
%     temp2=squeeze(temp2);
%     %
%     a=tensor(e*e');
%     b=tensor(mean(Y)); 
%     temp3=double(ttt(a,b));
%     temp3=squeeze(temp3); 
%     %
%     M3_b=M3_b+temp1+temp2+temp3;
% end
% M3_Y=M3_a-sigma2*M3_b; 
M3_aa=zeros(size(Y,2),size(Y,2),size(Y,2));
for t = 1:size(Y,1)
    for i = 1:size(Y,2)
        for j = 1:size(Y,2)
            for l = 1:size(Y,2)
                M3_aa(i,j,l)=M3_aa(i,j,l)+Y(t,i)*Y(t,j)*Y(t,l);
            end
        end
    end
end
M3_aa=M3_aa/size(Y,1); % M3_a = M3_aa
M3_bb=zeros(size(Y,2),size(Y,2),size(Y,2));
for t=1:size(Y,2)
    e=zeros(size(Y,2),1);
    e(t)=1;
    for i=1:size(Y,2)
        for j=1:1:size(Y,2)
            for l=1:1:size(Y,2)
                M3_bb(i,j,l)=M3_bb(i,j,l)+M1(i)*e(j)*e(l)+e(i)*M1(j)*e(l)+e(i)*e(j)*M1(l); % M3_b = M3_bb
            end
        end
    end
end
M3_Y=M3_aa-sigma2*M3_bb;
%% generate GMM moments: 3-order (whiten mu)
% [muu weightt]=eig(M2)
% mu_whiten = W'*muu;
% M3_mu=zeros(k,k,k);
% for n=1:k
%     for i=1:k
%         for j=1:k
%             for l=1:k
%                 M3_mu(i,j,l)=M3_mu(i,j,l)+weightt(n)*mu_whiten(i,n)*mu_whiten(j,n)*mu_whiten(l,n);
%             end
%         end
%     end
% end
%% Power Iteration for M3
%linewidth = 2;
% TOL = 1e-6; Maxiter = 10000;
% V_est = zeros(k,k); Lambda = zeros(k,1);
% for i = 1:k
%     %v_old = W'*mu(:,i) + 0.3 * rand(size(W'*mu(:,i)));%rand(k,1);
%      v_old = 10*rand(k,1);
%      v_old = v_old./norm(v_old);
% %     v1est_handle = vectarrow(origin, v_old, linewidth, '-', 'r'); 
% %             str = sprintf('Iter %d',0);
% %             text(v_old(1)-0.1,v_old(2)+0.01*0,str);
% %             hold on;
%     for iter = 1 : Maxiter    
%         %v_new = (Y'* ((Y* v_old).* (Y* v_old)))/size(Y,1); %%!!
%         v_new = double(ttv(tensor(M3), {v_old v_old},[1 2]));
% %         if i >1 
% %             for j = 1: i-1
% %                 v_new = v_new - (V_est(:,j) * (v_old'*V_est(:,j))^2)* Lambda(j);
% %             end          
% %         end
%         lambda = norm(v_new);
%         v_new = v_new./lambda;
%         % plot v_new at each iteration
% %         if i == 1
% %             v1est_handle = vectarrow(origin, v_new, linewidth, '-', 'r'); 
% %             str = sprintf('Iter %d',iter);
% %             text(v_new(1)-0.1,v_new(2)+0.02*iter,str);
% %             hold on;
% %         end
%         if norm(v_old - v_new) < TOL
%             fprintf('First %d eigenvector updates converged at iteration %d. \n', i, iter);
%             V_est(:,i) = v_new;
%             Lambda(i,1) = lambda;
%             break;
%         end
%         v_old = v_new;
%     end
% end 
% %% unwhitening/ recovering
% weight_est=(1./Lambda).^2; weight_est=weight_est/sum(weight_est);
% mu_est = pinv(W') * V_est * diag(Lambda);
% mu_est = rdivide(mu_est,repmat(sum(mu_est,1),size(mu,1),1));
%% Poly form of Moment
f_c = [];    
mpol x 2;
options.add_orth=1;
m = 3; % third dimensional tensor
p = 2; % Z-eigenvalue
f_A = M3_Y(1,1,1)*x(1)^3 + (M3_Y(1,2,1)+M3_Y(1,1,2)+M3_Y(2,1,1))*x(1)^2*x(2) + (M3_Y(1,2,2)+M3_Y(2,1,2)+M3_Y(2,2,1))*x(1)*x(2)^2 + M3_Y(2,2,2)*x(2)^3;
[lmd, eigvec, info] = AReigSTensors(f_A, f_c, mpol(x), m, p, options);

[weight_est mu_est]=para_recover(lmd,VecM,W)
%weight_est=(1./lmd).^2; %weight_est=weight_est/sum(weight_est);
%mu_est = pinv(W') * eigvec * lmd;