function printeigpair(myfid,timespan,lmd,eigvec)
% print cpu time, all eigvalues, all eigenvectors
    %% print to myfid file
    % print time 
    fprintf(myfid, '\n\ncpu time = %fs\r\n\n\n', timespan);  
    
    % print eigenvalues
    fprintf(myfid, 'eigvalues = \r\n');
    for i = 1:length(lmd)
        fprintf(myfid, '%10.4f',lmd(i));
    end
    fprintf(myfid,'\r\n\n');
    
    % print eigenvectors
    fprintf(myfid, '\neigvectors = \r\n');
    eigth = 1;
    for i = 1:size(eigvec,1)
        if i>1 && lmd(i)<lmd(i-1)
            eigth = eigth+1;
        end
        fprintf(myfid,' %d %10.4f : ',eigth, lmd(i));
        for j = 1:size(eigvec,2)
            fprintf(myfid, '%10.4f', eigvec(i,j));
        end
        fprintf(myfid, '\r\n\n');
    end
    fprintf(myfid,'\r\n\n');  
    
  %% print to command window  
        % print time 
    fprintf('\ncpu time = %fs\r\n', timespan);  
    
    % print eigenvalues
    fprintf('eigvalues = \r\n');
    for i = 1:length(lmd)
        fprintf('%10.4f',lmd(i));
    end
    fprintf('\r\n');
    
    % print eigenvectors
    fprintf('eigvectors = \r\n');
    eigth = 1;
    for i = 1:size(eigvec,1)
        if i>1 && lmd(i)<lmd(i-1)
            eigth = eigth+1;
        end
        fprintf(' %d %10.4f : ',eigth, lmd(i));
        for j = 1:size(eigvec,2)
            fprintf('%10.4f', eigvec(i,j));
        end
        fprintf('\n');
    end
    fprintf('\r\n\n');  
end