%===================================================================================
% Matrix-vector version of Gaussian elimination. This one does not use any
% triangular solves. Matrix A is overwritten by the combined L/U data. No pivoting
% is performed.
% 
% Parameter SMALL is the cut-off for safe inversion of a scalar. errflag is 0 if
% all's OK, otherwise has the row number for which the pivot element was too small.
% 
% This version is for m x n (possibly rectangular) matrices, with no pivoting.
%===================================================================================

function [A, pivots, errflag] = LUmatvec(A);

    [m, n] = size(A);
    lda = m;
    p = min([m,n]);
    SMALL = 1.0e3*eps;
    errflag = 0;
    pivots = 1:p;
    
    for k = 1:n-1   % proceed column by column
        A(k:m, k)   = A(k:m, k) - A(k:m, 1:k-1)*A(1:k-1, k);  % mat-vec multiply

        % Find the pivot row, bail out if pivot entry too small
        [s, I] = max(abs(A(k:m, k)));
        pivots(k) = I(1)+k-1;
        if (abs(A(k, k)) < SMALL)
            errflag = k;
            return;
        end
        % Swap pivot row with row k
        if (pivots(k) ~= k)
            A([ k pivots(k)], :) = A([pivots(k) k], :);
        end
        A(k, k+1:n) = A(k, k+1:n) - A(k, 1:k-1)*A(1:k-1, k+1:n);  % mat-vec multiply
        A(k+1:m, k) = A(k+1:m, k)/A(k, k); % scale the subdiagonal for next column

    end 
    % Next is only done if m > n
    A(n:m, n) = A(n:m, n) - A(n:m, 1:n-1)*A(1:n-1, n);
    % Next is only done if n+1 > m
    A(n+1:m, n) = A(n+1:m, n)/A(n, n);

return

