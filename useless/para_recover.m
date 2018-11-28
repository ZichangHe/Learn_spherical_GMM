function [weight mu] = para_recover(lmd,VecM,W) 
%VecM=[vec1|vec2|...]

weight =(1./lmd).^2; % weight=weight/sum(weight);
mu = pinv(W') * VecM * diag(lmd);
mu = rdivide(mu,repmat(sum(mu,1),size(VecM,1),1));