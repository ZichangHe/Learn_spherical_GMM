clear,clc
%% generate GMM variables
%mu = [1 2 3 4;-3 -1 2 2 ;2 1 3 3 ];
k=2; % num of component
d=3; % dimension of variables
dist = 40; % Distance of Gaussian!!
mu = -dist+(dist+dist)*rand(k,d);
sigma = diag([2 2 2]); % or [2 2 2]: shared diagonal covariance matrix
weight=[0.3 0.7];
gm = gmdistribution(mu,sigma,weight)
mu=mu'; %to fit [mu_1 mu_2 ... mu_n]
rng('default'); % For reproducibility
[Y, compIdx] = random(gm,1e+5);
% Y=[];
% for i=1:k %%%% Generate Y seperately
%     mean = mu(:,i)';
%     covariance = sigma;
%     n=1e+5;
%     mvn = mvnrnd(mean, covariance, n);
%     Y = [Y ; mvn];
% end
%% generate GMM moments: 1,2-order and Whitening matrix
M1=sum(Y)/size(Y,1);
M2_raw=zeros(size(Y,2)); 
for i=1:size(Y,1)
    M2_raw=M2_raw+Y(i,:)'*Y(i,:);
end
M2_raw=M2_raw/size(Y,1);
M2_cov=M2_raw-(sum(Y)/size(Y,1))'*(sum(Y)/size(Y,1)); % M2_raw-mean(Y)'*mean(Y)
%sigma2 = min(eig(cov(Y)));
sigma2_est = min(eig(M2_cov));
M2=M2_raw-sigma2_est*eye(size(Y,2));
[Uw, Lw, Vw] = svd(M2); 
W = Uw(:,1:k) * sqrt(pinv(Lw(1:k,1:k)));  %% check W'*M2*W = I
%% generate GMM moments: 3-order (whiten Y)
Y_whit=Y*W; %% Whitening!
M3_aa=zeros(size(Y_whit,2),size(Y_whit,2),size(Y_whit,2));
for t = 1:size(Y_whit,1)
    for i = 1:size(Y_whit,2)
        for j = 1:size(Y_whit,2)
            for l = 1:size(Y_whit,2)
                M3_aa(i,j,l)=M3_aa(i,j,l)+Y_whit(t,i)*Y_whit(t,j)*Y_whit(t,l);
            end
        end
    end
end
M3_aa=M3_aa/size(Y_whit,1); 
M3_bb=zeros(size(Y_whit,2),size(Y_whit,2),size(Y_whit,2));
M1_whit=M1*W;
for t=1:size(mu,1)
    e=zeros(1,size(mu,1));
    e(t)=1;
    e_whit=e*W;
    for i=1:size(Y_whit,2)
        for j=1:1:size(Y_whit,2)
            for l=1:1:size(Y_whit,2)
                M3_bb(i,j,l)=M3_bb(i,j,l)+M1_whit(i)*e_whit(j)*e_whit(l) + e_whit(i)*M1_whit(j)*e_whit(l) + e_whit(i)*e_whit(j)*M1_whit(l); 
            end
        end
    end
end
M3=M3_aa-sigma2_est*M3_bb;
%% Power Iteration (with/ without) constructing Tensor
TOL = 1e-8;                 % Convergence delta
maxiter = 10000;              % Maximum number of power step iterations
V_est = zeros(k,k);         % Estimated eigenvectors for tensor V
lambda = zeros(k,1);        % Estimated eigenvalues for tensor V
for i = 1:k
    % Generate initial random vector
    v_old = rand(k,1);
    v_old = v_old./norm(v_old);
    for iter = 1 : maxiter
        % Perform multilinear transformation
        %v_new = (Y_whit'* ((Y_whit* v_old) .* (Y_whit* v_old)))/size(Y,1);
        %v_new = v_new - sigma2_est * (W' * M1' * dot((W*v_old),(W*v_old)));
        %v_new = v_new - sigma2_est * (2 * W' * W * v_old * ((W'*M1')' * (v_old)));
        v_new = double(ttv(tensor(M3), {v_old v_old},[1 2]));
        % Defaltion
        if i > 1 
            for j = 1:i-1
                v_new = v_new - (V_est(:,j) * (v_old'*V_est(:,j))^2)* lambda(j); 
            end          
        end
        % Compute new eigenvalue and eigenvector
        l = norm(v_new);
        v_new = v_new./norm(v_new);
        % Check for convergence
        if norm(v_old - v_new) < TOL
            fprintf('Eigenvector %d converged at iteration %d. \n', i, iter);
            V_est(:,i) = v_new;
            lambda(i,1) = l;
            break;
        end
        v_old = v_new;
    end
end
%% unwhitening/ recovering
weight_est=(1./lambda).^2; %weight_est=weight_est/sum(weight_est);
mu_est = pinv(W') * V_est * diag(lambda);
%% Poly form of Moment
f_c = [];    
mpol x 2;
options.add_orth=1;
m = 3; % third dimensional tensor
p = 2; % Z-eigenvalue
f_A = M3(1,1,1)*x(1)^3 + (M3(1,2,1)+M3(1,1,2)+M3(2,1,1))*x(1)^2*x(2) + (M3(1,2,2)+M3(2,1,2)+M3(2,2,1))*x(1)*x(2)^2 + M3(2,2,2)*x(2)^3;
[lmd, eigvec, info] = AReigSTensors(f_A, f_c, mpol(x), m, p, options);
weight_est2=(1./lmd).^2; %weight_est=weight_est/sum(weight_est);
mu_est2 = pinv(W') * eigvec' * diag(lmd); %estimation results are completely the same!