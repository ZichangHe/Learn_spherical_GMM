clear,clc

%% generate GMM variables
mu = [1 2 3;-3 -5 2];
sigma = [1 1 1]; % shared diagonal covariance matrix
gm = gmdistribution(mu,sigma)
rng('default'); % For reproducibility
[Y,compIdx] = random(gm,1000);
% %% verfication
% M2_t=cov(Y); % theory
% M2_m=zeros(size(Y,2)); 
% for i=1:size(Y,1)
%     M2_m=M2_m+Y(i,:)'*Y(i,:);
% end
% M2_m=M2_m/size(Y,1);
% M2_m=M2_m-mean(Y)'*mean(Y); % test
%% generate GMM moments
M2=cov(Y)-eye(size(Y,2));
M3_a=zeros(size(Y,2),size(Y,2),size(Y,2));
for i=1:size(Y,1)
    a=tensor(Y(i,:)'*Y(i,:));
    b=tensor(Y(i,:));
    temp=double(ttt(a,b));
    temp=squeeze(temp); 
    M3_a=M3_a+temp;
end
M3_a=M3_a/size(Y,1);
% M3_aa=zeros(size(Y,2),size(Y,2),size(Y,2));
% for t = 1:size(Y,1)
%     for i = 1:size(Y,2)
%         for j = 1:size(Y,2)
%             for k = 1:size(Y,2)
%                 M3_aa(i,j,k)=M3_aa(i,j,k)+Y(t,i)*Y(t,j)*Y(t,k);
%             end
%         end
%     end
% end
% M3_aa=M3_aa/size(Y,1); % M3_a = M3_aa
M3_b=zeros(size(Y,2),size(Y,2),size(Y,2));
for i=1:size(Y,2)
    e=zeros(size(Y,2),1);
    e(i)=1;
    %
    a=tensor(mean(Y)'*e');
    b=tensor(e);
    temp1=double(ttt(a,b));
    temp1=squeeze(temp1);
    %
    a=tensor(e*mean(Y));
    b=tensor(e);
    temp2=double(ttt(a,b));
    temp2=squeeze(temp2);
    %
    a=tensor(e*e');
    b=tensor(mean(Y)); 
    temp3=double(ttt(a,b));
    temp3=squeeze(temp3); 
    %
    M3_b=M3_b+temp1+temp2+temp3;
end
M3=M3_a-M3_b;
%% Poly form of Moment
sym x1 x2 x3;
%f_A=M3_poly(Y,x);
%% Moment decomposition
  