%-----------------------------------------------------------------------------------
%
% BLAS-3 version of GE with pivoting. Based on matrix-matrix multiply and
% solving triangular system with multiple right hand sides. No error checking
% is performed. Block size is nu, system size is n x n. A is overwritten.
% 
% This calls the function LUmatvec, which does LU factorization on the block column.
% 
% This function only does the factorization, not the actual solves for the
% solution vector.
%
%-----------------
% Randall Bramley
% Department of Computer Science
% Indiana University, Bloomington
%-----------------
% Started: Thu Nov 13 09:28:53 EST 1997
% Modified: Wed 21 Oct 2009, 10:58 AM
% Modified: Mon 25 Jun 2012, 02:20 PM formatting
% Last Modified: Mon 13 Oct 2014, 09:58 AM
%-----------------------------------------------------------------------------------

function [A, piv, errflag] = LUmatmat(A, nu);

    [m,n] = size(A);
    piv = 1:n;       % must return *something* for piv if early termination occurs
    piv = piv(:);    % make sure piv is a column vector

    nu = round(nu);  % Guard against dumb things
    nu = max(nu,1);  % Guard against profoundly dumb things
    errflag = 0;

    if (m ~= n)
        disp('This block LU is only for square matrices, but matrix A is')
        disp(['of size ', num2str(m), ' by ',  num2str(n)])
        errflag = -1;
        return;
    end 

    for k=1:nu:n
        cols = k:min(k+nu-1,n); % indices for current block column
        nub = min(nu,n-k+1); % number of columns in current block column
        % apply previous transformations to current block column
        A(k:n,cols) = A(k:n,cols) - A(k:n,1:k-1)*A(1:k-1,cols);
        % LU factorization on current block column
        [A(k:n,cols),piv(cols),errflag]  = LUmatvec(A(k:n,cols));

        if (errflag ~= 0)
            disp('Evil things happened on a LUmatvec call');
            errflag = errflag + k-1;
            disp(['errflag is nonzero and equal to value ', num2str(errflag)]) ;
            return
        end

        piv(cols) = piv(cols) + k - 1; % Adjust pivots to full matrix indexing

        other_cols = [1:k-1  k+nub:n]; % Apply pivots to rest of matrix
        for j = cols
          A([piv(j) j], other_cols) = A([j piv(j)], other_cols);
        end 

        if (k+nub <= n) % Update trailing part of matrix

            tail = k+nub:n; % Indices for trailing part of matrix

            % Update top block of nub rows
            A(cols,tail) = A(cols,tail) - A(cols,1:k-1)*A(1:k-1,tail);

            % Solve triangular system for top block of num rows
            A(cols,tail) = (eye(nub)+tril(A(cols,cols),-1)) \ A(cols,tail);

        end  % if (k+nub <= n)
    end  % for k=1:nu:n

return
