function f_A = symtensor2poly(A,x,m,n)
%% file name: symtensor2poly.m
%       Given a symmetric tensor A, compute its corresponding polynpmial f_A
%
% Input:  A	 --   m-degree n-dimensional symmetric tensor
%         x  --   n-dimensioal mpol vector
% output: f_A  -- m-degree homogenous polynomial coressponding to A
%
%%
   
    f_A = 0;
    alpha = enumhomopoly(n,m);
    for kalpha = 1:size(alpha,1)
        s = alpha(kalpha,:);
        As = symtensor2vec(s,m,n);
        pA = 1;
        for ipA=1:n
            pA = pA*x(ipA)^s(ipA);
        end
        f_A = f_A + coeff(s,m,n)*A(As)*pA;
    end
end

%%
%%
function alpha = enumhomopoly(dim,degree)
%% enmurate all  (k = nchoosek(dim+degree, dim)) homogenous polynomial terms
% input:     dim   
%            degree
% output:    alpha  --  2 dimensional array in  k x n.
%            for example m=3, n =2
%                        alpha = [3 0; 2 1; 1 2; 0 3]
%%

    alpha = [];
    s = [];
    while 1
        while length(s)~= (dim-1)
            s=[s 0];
        end
        s = [s degree-sum(s)];
        alpha = [alpha; s];  
        
        smin = degree-sum(s(1:end-1));
        while (s(end) == smin)||(length(s)== dim)
            s(end)=[];
            if isempty(s)
                return;
            end
            smin = degree-sum(s(1:end-1));
        end
        s(end)=s(end) + 1;
    end
    
end

%%
%%
function As = symtensor2vec(s,m,n)
%% 
%input:  s  --  n-dimensional vector, with s(i) denote the degree of x(i)
%        m  -- degree of homogenous polynomial, sum of s equals to m
%        n  -- dimension of s
%output: As -- a positive integer denotes the label for transform A to vector
%        for example   m=3,n=2
%         s = [3 0]; A(1 1 1) As = 1 + 0 + 0 + 0=1; 
%         s = [2 1]; A(1 1 2) As = 1 + 0 + 0 + 4 = 5;
%         s = [1 2]; A(1 2 2) As = 1 + 0 + 2 + 4 = 7;
%         s = [0 3]; A(2 2 2) As = 1 + 1 + 2 + 4 = 8;
%        its r.p.t. vector number is 1 + (3-1) = 3
%%

    if (length(s)~= n) || (sum(s)~=m)
        error('Error in compute\n');
    end
    As = 1;
    kAs = 0;
    for iAs = n:-1:1
        for Asj = 1:s(iAs)
            As = As + (iAs-1)*n^kAs;
            kAs = kAs + 1;
        end
    end
end

%%
%%
function coeff = coeff(s,m,n)
%%
% input: s -- n dimension vector, with s(i) denote the degree of x(i)
%        m -- degree of homogenous polynomial, sum of s equals to m
%        n -- dimension of s
%output: coeff -- the coeff in (e'x)^m
% for example m=3,n=2, (e'x)^m = x1^3 + 3x1^2x2 + 3x1x2^2 + x2^3
% s = [3 0]  coeff = 1; s = [2 1] coeff = 3
%%

    if (length(s)~= n) || (sum(s)~=m)
        error('Error in compute\n');
    end
    coeff = 1;
    localm = m;
    for ci = 1:n-1
        si = s(ci);
        if si ~= 0
            coeff = coeff*nchoosek(localm, si);
            localm = localm - si;
        end
    end
end
%%         
                 
                


            
