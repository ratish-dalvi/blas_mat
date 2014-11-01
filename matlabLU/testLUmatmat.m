%------------------------------------------------------------------------------------
% Test driver script for LUmatmat. Vary n and nu to check different configurations.
% The Boolean diagdom = true will create a system that does not require pivoting,
% but the array piv needs to be set to 1:n if that is done.
%------------------------------------------------------------------------------------

n        =   8;     % matrix order
nu       =   3;     % blocksize
diagdom  = false;   % make system diagonally dominant

A   = randn(n, n);
if diagdom
    for k = 1:n
        A(k,k) = 1000;
    end
end 

Aorig = A;
t0 = clock;
    [A, piv, errflag] = LUmatmat(A, nu);
t1 = clock;
timing = etime(t1, t0);
if (errflag ~= 0)
    disp(sprintf('Bad news from LUmatmat; errflag = %d', errflag))
    disp('no point in continuing; try again with a different system')
    return
end
disp(sprintf('Time required: %g seconds', timing));

% Can compute how many flops are taken by looking at the algorithm, 
% but it is known to be around (2/3)n^3, so just use that.
mflops = (0.666666666666666666666666667e-6)*(n^3);
mfloprate = mflops/timing;
disp(sprintf('Computational rate: %g Mflop/second', mfloprate));

% Correctness check of the factorization
L = eye(n) + tril(A, -1);
U = triu(A);

for k = 1:length(piv)
    Aorig([piv(k) k], :) = Aorig([k piv(k)], :);
end  

% Use the 1, infinity, or Frobenius norm for matrices. Do not use the 2-norm for
% large matrices. To reflect the extra numerical loss of precision caused by larger
% values of n. Because the one-norm will grow O(n^2), scale by that.
err = norm(Aorig - L*U, 1)/n^2;
disp(sprintf('Scaled 1-norm of PLU - A: %e', err));

