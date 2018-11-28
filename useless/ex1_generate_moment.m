%% ex1
%%
%%
    clc;
    clear all;
    
    filename = strcat('result\','Result_ex1', '.txt');
    myfid = fopen(filename,'w+');
    
%%
%% initialization
%%  
    begintime = tic;
    m = 3; % third dimensional tensor
    n = 3; 
    p = 2; 
    if exist ('mpol') ~= 2
        error('GloptiPoly 3 is not properly installed')
    end
    mpol x 3;
    options = [];
    
%% generate the polynomail
  % generate GMM variables
    mu = [1 2 3 ;-3 -5 2];
    sigma = [1 1 1]; % shared diagonal covariance matrix
    gm = gmdistribution(mu,sigma);
    rng('default'); % For reproducibility
    [Y,compIdx] = random(gm,1000);
    f_A = M3_poly(Y,x);
    %%
%% solving process
%%
    %f_A = x(1)^4 + 2*x(2)^4 + 3*x(3)^4;
    f_c = [];
    profile on
    [lmd, eigvec,  info] = AReigSTensors(f_A, f_c, mpol(x), m, p, options);
    
%%
%% output information to my file
%%       
    printinform(myfid);
    
    timespan = toc(begintime);
    printeigpair(myfid,timespan,lmd,eigvec);

%%
%% post-process
%%   
    fclose(myfid);
    
    mset clear;
    clear all; 

 