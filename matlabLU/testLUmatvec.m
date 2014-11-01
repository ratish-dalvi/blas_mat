%------------------------------------------------------------------------------------
% Test driver script for LUmatvec.
%------------------------------------------------------------------------------------

m           =   8;     % order of matrix 
n           =   5;     % blocksize
diagdom     = false;   % make system diagonally dominant

A   = randn(m, n);
if diagdom
    for k = 1:n
        A(k,k) = 1000;
    end
end 

Aorig = A;
t0 = clock;
    [A, piv, errflag] = LUmatvec(A);
t1 = clock;
timing = etime(t1, t0);
if (errflag ~= 0)
    disp(sprintf('Bad news from LUmatvec; errflag = %d', errflag))
    disp('no point in continuing; try again with a different system')
    return
end
disp(sprintf('Time required: %g seconds', timing));

% Can compute how many flops are taken by looking at the algorithm, 
% but it is known to be around (2/3)n^3, so just use that.
mflops = (0.666666666666666666666666667e-6)*(n^3);
mfloprate = mflops/timing;
disp(sprintf('Computational rate: %g Mflop/second', mfloprate));

% Correctness checking of the factorization
L = tril(A, -1);
for i = 1:min([m,n])
    L(i,i) = 1.0;
end
U = triu(A);
U = U(1:n,1:n);   % Force U to be square

for k = 1:length(piv)
    Aorig([piv(k) k], :) = Aorig([k piv(k)], :);
end  

% Use the 1, infinity, or Frobenius norm for matrices. Do not use the 2-norm for
% large matrices. To reflect the extra numerical loss of precision caused by larger
% values of n. Because the one-norm will grow O(n^2), scale by that.
err = norm(Aorig - L*U, 1)/n^2;
disp(sprintf('Scaled 1-norm of PLU - A: %e', err));

