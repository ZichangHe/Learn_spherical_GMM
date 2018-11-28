
    function [lmd, eigvec,  info] = AReigSTensors(f_A, f_c, x, m, p, options)
    
 %filename: AReigSTensors.m
%           To   compute all the real eigenpairs of symmetric tensor by
%                       max  f_A      s.t.   f_c
%Input:
%     f_A     --   objective function, homogenous polynomial of tensor A
%     f_c      --  initial constraints array;  if none, input []
%     x         --   mpol var x with degree n of polynomial f_A or tensor A
%     m        --   order of  tensor A
%     p         --   type of tensor desired (p = 2 for Z-eig and p = 4 for H-eig)
%     options
%              options.tol:          --  tolrence for two distinct eigenvalues, default = 1e-4
%              options.dt0:         --  initial gap of eigenvalue for lmd(k) and lmd(k+1), default=0.05
%              options.N0:          --  initial order of relaxation,  default = ceil((m+p-2)/2)
%              options.Nmax:     --  maximum of N, default = 11-n
%              options.printinf:   --  print detailed information or not. 1 yes 0 no, default = 0
%              options.maxeignum:   -- maximum number of eigenvalues desired (with mutltiplicity)
%
%Output:
%      lmd      --  eigenvalues
%      eigvec  --  eigenvector array
%      info      --  a struct array containing basic information such as
%            info.success      -- 1 if success; 0 if terminate by Nmax; -1 if  infeasible
%            info.time           -- record computational time
%            info.Nrecord     -- record the order N used to comput lmd(i)
%
% Useful plugins including
%         glopotipoly 3
%         sedumi_1_3   
%
% ref:      C. Cui, Y. Dai and J. Nie: All real eigenvalues of symmetric tensors, 
%            SIAM J. Matrix  Anal. Appl. 35 (2014) 1582¨C1601
% @authors: Chunfeng Cui, Yu-Hong Dai, and Jiawang Nie

%%
%% initialize output information
%%

    algstart = tic;
    success = 1;
    lmd = []; 
    eigvec = []; 
    Nrecord = [];
    
%%
%% default information 
%%

    if nargin < 5
        error('Error: you should input at least  5 parameters \n');
    end
    n = length(x);
    if nargin < 6
        options = [];
    end
    if  isfield(options, 'tol')
        tol = options.tol;
    else
        tol = 1e-4;
    end
    if  isfield(options, 'dt0')
        dt0 = options.dt0;
    else
        dt0 = 0.05;
    end
    if  isfield(options, 'N0')
        N0 = options.N0;
    else
        N0 = ceil((m+p-2)/2);
    end
    if  isfield(options, 'Nmax')
        Nmax = options.Nmax;
    else
        Nmax = 10; %Relaxation order
    end
    
    if  isfield(options, 'add_orth')
        add_orth = options.add_orth;
        eps_orth = 1e-6; %threshold for orthogonality check
    else
        add_orth = 0;
    end
    
    if  isfield(options, 'printinf')
        printinf = options.printinf;
    else
        printinf = 1; 
    end
    
%%
%% define the constraints
%%
    minors = []; 
    sumidx = [];
    for i = 1 : n
        for j = i+1 : n
            minors = [minors; diff(f_A,x(i))*x(j)^(p-1)-diff(f_A,x(j))*x(i)^(p-1)];
            sumidx = [sumidx; i+j];
        end
    end

    for t = 1 : 2*n-3
        rowI =  sumidx == t+2;
        f_c = [f_c, sum( minors(rowI)) == 0];
    end
    
    f_c = [f_c, sum(x.^p)-1==0];
    if (p == 4) || (p == 2 && mod(m,2) == 0)
        f_c = [f_c, sum(x) >= 0];
    elseif p==2 && mod(m,2)~=0
        f_c = [f_c, f_A >=0];
    end
    
%%
%% firstly,  solve a hierarchy of SDP and get the largest eigvalue lmd1
%%
    
    N = N0; 
    localstart = tic;
    
    % transform the polynomial problem to SDP form by GloptiPoly3    
    if isfield(options, 'infty') && options.infty
        options.maxeignum = 1;
        s = rng;
        rng(s);
        f_AInfty = randn(n,1)'*x;
        P1 = msdp(max(f_AInfty), f_c, N);
    else
        P1 = msdp(max(f_A), f_c, N);
    end
    
    % solve by GloptiPoly3 - sedumi
    [sta, lmdp] = msol(P1);

    if printinf
        fprintf('\n to compute lmd 1:\n');
        fprintf('\n N = %d \n sta = %d \n', N, sta); 
        fprintf(' lmd = %f \n time = %f \r\n', lmdp, toc(localstart));
    end

    if sta == -1;  
        success = -1; 
        info.success = success;
        info.time = toc(algstart);
        info.Nrecord = Nrecord;
        return;
    end

    while sta ~=1 && N<Nmax
        N = N+1;        % increase N
        localstart = tic;
        if isfield(options, 'infty') && options.infty
            s = rng;
            rng(s);
            f_AInfty = randn(n,1)'*x;
            P1 = msdp(max(f_AInfty), f_c, N);
        else
            P1 = msdp(max(f_A), f_c, N);
        end
        [sta, lmdp] = msol(P1);
        if printinf
            fprintf('\n to compute lmd 1:\n');
            fprintf('\n N = %d \n sta = %d \n', N, sta); 
            fprintf(' lmd = %f \n time = %f\n\n', lmdp, toc(localstart));
        end
    end
    
    if sta~=1   
        if sta ==0
            success = 0;
            Nrecord = [Nrecord N];
            lmd = [lmd lmdp];
            info.success = success;
            info.time = toc(algstart);
            info.Nrecord = Nrecord;
            return;       
        else
            success = -1;
            info.success = success;
            info.time = toc(algstart);
            info.Nrecord = Nrecord;
            return;    
        end
    else
        %  update the output infomation
        Nrecord = [Nrecord N];
        xk = double(x);
        tk = size(xk, 3); 
        for k = 1 : tk
            eigvec = [eigvec; xk(:,:,k)'];
            lmd = [lmd lmdp];
            if printinf
                fprintf('\n output eigpair for lmd 1:\r\n');
                fprintf(' eigenvalue = \n%10.4f \n',lmdp);
                fprintf(' eigenvector =\n');
                for ki = 1:n
                     fprintf('%10.4f ', xk(ki,1,k));
                end
                fprintf('\n\n');
            end
        end
    end
    
%%
%% two steps to compute lmdk: 
%% 1, compute dt; 2, solve a hierarchy of SDP to get lmdk
%%

    while success 
        eignum = numel(unique(lmd)); % eignum = length(lmd);
        if isfield(options, 'maxeignum') && eignum >= options.maxeignum
            break;
        end
        lmdold = lmd(end);    
        decresedt = 1;  
        N = N0; 
        dt =dt0; 
        
%% 
%% first, decrease dt until lmdk = chklmdk, where chklmdk is from a hierachy of SDP 
%%
        
        while decresedt  && dt>=tol   % check dt
            localstart = tic;
            Pchk = msdp(min(f_A), [f_c  f_A >= lmdold - dt], N);
            [sta, chklmdk] = msol(Pchk);
            if printinf
                fprintf('\n to compute chklmd %d:\n', eignum);
                fprintf('\n N = %d \n sta = %d \n', N, sta); 
                fprintf(' lmd = %10.4f \n time = %10.4f \r\n', chklmdk, toc(localstart));
            end    

            if sta == -1
                error('Error: unable to solve SDP\n');
            end
            while sta ~= 1 && N<Nmax % increase N          
                N = N+1;
                localstart = tic;
                Pchk = msdp(min(f_A), [f_c  f_A>= lmdold - dt], N);
                [sta, chklmdk] = msol(Pchk);
                if printinf
                    fprintf('%%*********************************************\r\r\n\n');
                    fprintf('\n to compute chklmd %d:\n', eignum);
                    fprintf('\n N = %d \n sta = %d \n', N, sta); 
                    fprintf(' lmd = %10.4f \n time = %10.4f \r\n', chklmdk, toc(localstart));
                end 
            end

            gap = lmdold - chklmdk;

            if gap <= -10*tol % return back
                error(' reterive back some eigenvalues. \r\n');
            elseif gap <= tol      
                decresedt = 0;
            else             
                dt = dt/5;        % dt = min(dt/5, abs(lmdold - lmd_chk)+1e-4 );
                if printinf 
                    fprintf(' dt = %10.4f\n N = %d \n\n', dt,N);
                end
            end
        end %while decrease dt
%%
%% second, solve a hierachy of SDP to obtain lmdk. 
%%
        N = N0; 
        localstart = tic; 
        if add_orth 
            f_c_orth = [f_c, (lmdold - dt)-f_A>=0];
            for i_orth = 1:size(eigvec,1)
                f_c_orth = [f_c_orth, -eps_orth<=eigvec(i_orth,:)*x, eps_orth>=eigvec(i_orth,:)*x];
            end
            Pk = msdp(max(f_A), f_c_orth, N);
        else
            Pk = msdp(max(f_A), [f_c, (lmdold - dt)-f_A>=0], N);
        end
        [sta, lmdp] = msol(Pk);
        if printinf
            fprintf('\n to compute lmd %d:\r\n', eignum+1);
            fprintf('\n N = %d \n sta = %d \n', N, sta); 
            fprintf(' lmd = %10.4f \n time = %10.4f \r\n', lmdp, toc(localstart));
        end 
        if sta == -1;  
            success = -1; 
            info.success = success;
            info.time = toc(algstart);
            info.Nrecord = Nrecord;
            return;
        end

        while sta ~= 1 && N<Nmax
            N = N+1; 
            localstart = tic; 
            if add_orth 
                f_c_orth = [f_c,(lmdold- dt)-f_A>=0];
                for i_orth = 1:size(eigvec,1)
                    f_c_orth = [f_c_orth, eps_orth>=eigvec(i_orth,:)*x,-eps_orth<=eigvec(i_orth,:)*x];
                end
                P = msdp(max(f_A), f_c_orth, N);
            else
                P = msdp(max(f_A), [f_c,(lmdold- dt)-f_A>=0], N);
            end
            [sta, lmdp] = msol(P);
            if printinf
                fprintf('\n to compute lmd %d :\n', eignum+1);
                fprintf('\n N = %d \n sta = %d \n', N, sta); 
                fprintf(' lmd = %10.4f \n time = %10.4f \n\n', lmdp, toc(localstart));
            end 
        end
        %  update the output infomation
        if sta ~=1
            if sta == 0
                success = 0;
                Nrecord = [Nrecord N];
                lmd = [lmd lmdp];
                info.success = success;
                info.time = toc(algstart);
                info.Nrecord = Nrecord;
                return;   
            else
                success = -1;
                info.success = success;
                info.time = toc(algstart);
                info.Nrecord = Nrecord;
                return; 
            end
        else            
            xk = double(x);
            tk = size(xk, 3); 
            for k = 1 : tk
                eigvec = [eigvec; xk(:,:,k)'];
                lmd = [lmd lmdp];
                if printinf
                    fprintf('\n output  lmd %d :\n', eignum+1);
                    fprintf(' eigenvalue = \n %10.4f \n',lmdp);
                    fprintf(' eigenvector =\n');
                    for ki = 1:n
                         fprintf('%10.4f ', xk(ki,1,k));
                    end
                    fprintf('\n\n');
                end
            end  
        end
    end
    info.success = success;
    info.time = toc(algstart);
    info.Nrecord = Nrecord;
    end
 %%**********************************************************************