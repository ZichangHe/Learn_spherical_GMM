%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Power Method for Tensor Decomposition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clear; clc;
%%% (1) Exchangeable Hidden Variable Model
%%%                  h
%%%                / | \ 
%%%              A/ A| A\
%%%              /   |   \
%%%            x1   x2    x3 
%%%  
%%% (1a) 
%%% Emission Matrix A \in R^{d \times k}: 
%%% E[x1|h] = A*h, E[x2|h] = A*h, E[x3|h] = A*h.
%%%
%%% h is a vecotr which represents the hidden state
%%% Observations x1,x2, and x3 \in R^{d}
%%% d is the number of possible observation
%%% "one hot encoding":
%%% x \in {0,1}^d, just one entry which corresponds to the observation
%%% is equal to 1 and the rest are zero
%%%
%%% Hidden Variable h \in 2 dimensional basis ={[0;1],[1;0]}. k = 2

d = 3;k = 2;
n = 100000;%% the number of trials
rand('seed',1);
h = zeros(n,2);
%%% Matrices of samples, each row corresponds to one sample of x1, x2, x3
x1_samples= zeros(n,d); x2_samples= zeros(n,d); x3_samples= zeros(n,d); 
A_true = rand(d, k); 
A_true = rdivide(A_true,repmat(sum(A_true,1),d,1));

%%%
%%% (1b) Generate Synthetic Samples
%%%
for t = 1 : n
    % generate h for this sample
    this_h_category = (rand()>0.5) + 1;
    h (t,this_h_category) = 1;
    % generate x1 for this sample | h
    this_transition_cum = cumsum(A_true(:,this_h_category));
    this_x_category = find(this_transition_cum > rand(),1);
    x1_samples(t,this_x_category) = 1;
    % generate x2 for this sample | h
    this_x_category = find(this_transition_cum > rand(),1);
    x2_samples(t,this_x_category) = 1;
    % generate x3 for this sample | h
    this_x_category = find(this_transition_cum > rand(),1);
    x3_samples(t,this_x_category) = 1;
end

%%%                                       n
%%% (2) Moments: Second Order Moment M2 = \sum x1(t,:) \otimes x2(t,:) 
%%%                                      t=1
%%%
%%% (3) Whitening matrix W such that W' M2 W = I;
%%%  Anandkumar et al. (2014) page 21
%%% 
fprintf('The second order moment: ');
M2 = x1_samples'*x2_samples/n
%%% one way of computing Whiteinig Matrix: W = U D^(-0.5)
%%% U \in R^{d \times k} is the matrix of orthonormal eigenvectors of M2
%%% D \in R^{k \times k} is the diagonal matrix of positive eigenvalues of M2
[Uw, Lw, Vw] = svd(M2); 
fprintf('The second order moment singular values:'); Lw
W = Uw(:,1:k)* sqrt(pinv(Lw(1:k,1:k)));
[Uw_true, Lw_true, Vw_true] = svd(A_true*A_true');
W_true = Uw_true(:,1:k)* sqrt(pinv(Lw_true(1:k,1:k)));
V_true = W_true' * A_true;
linewidth = 3;
figure()
subplot(1,2,1)
origin = zeros(d,1);
for ind_k = 1: k
    htrue_handle = vectarrow(origin, A_true(:,ind_k), linewidth, '-', 'g'); 
    hold on;
end
title('Original Bases a1, a2')
subplot(1,2,2)
origin = zeros(k,1);
for ind_k = 1: k
    Vtrue_handle = vectarrow(origin', V_true(:,ind_k)', linewidth, '-', 'r'); 
    hold on;
end
title('Whitened Bases v1, v2')
set(gcf, 'color','w');
%%%
%%%                                       n
%%% Third Order Moment of raw data M3 = \sum x1(t,:) \otimes x2(t,:) \otimes x3(t,:) 
%%%                                      t=1
%%% (3a) Whiten the data
%%%

y1 = x1_samples * W; %%% W^T*x1
y2 = x2_samples * W; %%% W^T*x2
y3 = x3_samples * W; %%% W^T*x3

%%%                                       n
%%% (4) Moments: Third Order Moment T = \sum y1(t,:) \otimes y2(t,:) \otimes y3(t,:)
%%%                                      t=1
%%% Which is equal to T = M3(W,W,W)
T = zeros(k,k,k);
for t = 1 : n
    for l = 1 : k
        for j = 1 : k
            for i = 1 : k
                T(i,j,l) = T(i,j,l) + y1(t,i)*y2(t,j)*y3(t,l);
            end
        end
    end
end
fprintf('The third order moment: '); T = T / n %% Empirical Average 
%%%                             k
%%% (5) Power Iteration:  T = \sum \lambda_i V(:,i) \otimes V(:,i) \otimes V(:,i)
%%%                            i=1

figure()
origin = zeros(k,1);
v1true_handle = vectarrow(origin, V_true(:,1), linewidth, '-', 'g'); 
hold on;

linewidth = 2;
TOL = 1e-6; Maxiter = 10;
V_est = zeros(k,k); Lambda = zeros(k,1);
for i = 1:k
    v_old = V_true(:,1)+ 0.3* rand(size(V_true(:,1)))%rand(k,1);
    v_old = v_old./norm(v_old);
    v1est_handle = vectarrow(origin, v_old, linewidth, '-', 'r'); 
            str = sprintf('Iter %d',0);
            text(v_old(1)-0.1,v_old(2)+0.01*0,str);
            hold on;
    for iter = 1 : Maxiter    
        v_new = (y1'* ((y2* v_old) .* (y3* v_old)))/n;
        if i > 1 
            for j = 1: i-1
                v_new = v_new - (V_est(:,j) * (v_old'*V_est(:,j))^2) * Lambda(j);
            end          
        end
        lambda = norm(v_new);
        v_new = v_new./norm(v_new);
        % plot v_new at each iteration
        if i == 1
            v1est_handle = vectarrow(origin, v_new, linewidth, '-', 'r'); 
            str = sprintf('Iter %d',iter);
            text(v_new(1)-0.1,v_new(2)+0.02*iter,str);
            hold on;
        end
        if norm(v_old - v_new) < TOL
            fprintf('First %d eigenvector updates converged at iteration %d. \n', i, iter);
            V_est(:,i) = v_new;
            Lambda(i,1) = lambda;
            break;
        end
        v_old = v_new;
    end
end
%%%
%%% (6) Unwhitening and Compare A_true and V_est
%%%
%%%



display('True parameter A_true vs Estimated parameter A_est: ');
A_est = pinv(W') * V_est * diag(Lambda);
A_est = rdivide(A_est,repmat(sum(A_est,1),d,1))
A_true
figure()
origin = zeros(d,1);
for ind_k = 1: k
    htrue_handle = vectarrow(origin, A_true(:,ind_k), linewidth, '-', 'g'); 
    hold on;
end

for ind_k = 1: k
    hest_handle = vectarrow(origin, A_est(:,ind_k), linewidth, ':', 'r'); 
    hold on;
end
legend([htrue_handle,hest_handle], 'True A','Estimated A', 'Location','North')
legend('boxoff')
set(gcf, 'color','w');
