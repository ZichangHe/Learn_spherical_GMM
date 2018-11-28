%% ex1
%%
%%
    clc;
    clear all;
    
    filename = strcat('result\','Result__ex1', '.txt');
    myfid = fopen(filename,'w+');
 
%%
%% initialization
%%  
    begintime = tic;
    m = 4;
    n = 3;
    p = 2;
    if exist ('mpol') ~= 2
        error('GloptiPoly 3 is not properly installed')
    end
    mpol x 3;
    options = [];
%%
%% solving process
%%
    f_A = x(1)^4 + 2*x(2)^4 + 3*x(3)^4;
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
 