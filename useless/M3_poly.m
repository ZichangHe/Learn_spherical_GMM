function f_A = M3_poly(Y,x)
        k=size(Y,2);
        M2=cov(Y)-eye(size(Y,2));
        [Uw, Lw, Vw] = svd(M2); 
        W = Uw(:,1:k)* sqrt(pinv(Lw(1:k,1:k)));
        f_A = 0;
        Y_whiten=Y*W;
        %Y_whiten=Y;
        sigma2 = data2variance(Y_whiten);
        for i_smp = 1:size(Y_whiten,1)
            y = Y_whiten(i_smp,:);
            f_A = f_A + (y*x)^3-3*sigma2*(y*x)*sum(x.^2);
        end
%       f_A = sum((Y*x).^3)-3*sigma2*sum(Y*x)*sum(x.^2);
        f_A = f_A/size(Y,1); 

        function sigma2 = data2variance(Y)
%            M1 = mean(Y,1);
%            M2 = (Y'*Y)/size(Y,1) - M1'*M1;
            M2 = cov(Y);
            sigma2 = min(eig(M2));
        end

        
    end