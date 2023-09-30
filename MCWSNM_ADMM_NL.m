function  [Z, sigma] =  MCWSNM_ADMM_NL( Y, SigmaRow,SigmaCol, Par )
% This routine solves the following weighted nuclear norm optimization problem with column weights,
%
% min_{X, Z} ||W(Y-X)||_F^2 + ||Z||_w,*  s.t.  X = Z
%
% Inputs:
%        Y        -- 3p^2 x M dimensional noisy matrix, D is the data dimension, and N is the number of image patches.
%        NSig   -- 3p^2 x 1 dimensional vector of weights
%        Par     -- structure of parameters
% Output:
%        Z        -- 3p^2 x M dimensional denoised matrix
%        sigma -- the noise standard deviation

% tol = 1e-8;
if ~isfield(Par, 'maxIter')
    Par.maxIter = 10;
end
if ~isfield(Par, 'rho')
    Par.rho = 1;
end
if ~isfield(Par, 'mu')
    Par.mu = 1;
end
if ~isfield(Par, 'display')
    Par.display = true;
end
%% Initializing optimization variables
% intialize

% Initializing optimization variables
% Intialize the weight matrix W
if Par.lambda2 == 0
mNSig = min(SigmaRow);
W = (mNSig+eps) ./ (SigmaRow+eps);
else
% W = (mNSig+eps) ./ (SigmaRow+eps);
mNSig = sqrt(min(SigmaRow));
W1 = mNSig ./ (sqrt(SigmaRow)+eps);
% W1 = exp( - Par.lambda*mean(SigmaRow, 2));
W2 = 1 ./ (sqrt(SigmaCol) + eps); 
end

Z = zeros(size(Y));
A = zeros(size(Y));
%% Start main loop
iter = 0;
PatNum       = size(Y,2);
C = Par.Constant*sqrt(PatNum)*(mNSig^2);
% C = Par.Constant*sqrt(PatNum);
p = Par.p;
while iter < Par.maxIter
    iter = iter + 1;
    
    % update X, fix Z and A
    % min_{X} ||W * Y - W * X||_F^2 + 0.5 * rho * ||X - Z + 1/rho * A||_F^2
    if Par.lambda2 == 0 
        X = diag(1 ./ (W.^2 + 0.5 * Par.rho)) * (diag(W.^2) * Y + 0.5 * Par.rho * Z - 0.5 * A);
        else
        AA = diag(W1.^2);
        W2inv = diag(1./(W2.^2));
        BB = 0.5 * Par.rho * W2inv;
        CC = diag(W1.^2) * Y + 0.5 * ( Par.rho * Z - A) * W2inv;
        X = sylvester(AA, BB, CC);
    end
    
    % update Z, fix X and A
    % min_{Z} ||Z||_*,w + 0.5 * rho * ||Z - (X + 1/rho * A)||_F^2
    X_A = X + A/Par.rho;
	[U, SigmaY, V] = svd(full(X_A), 'econ');
	s = diag( SigmaY ); % diaglog matrix -> vector
	s1= zeros( size(s) );
	SIGMA = diag( SigmaY ); 
%         SIGMA = sqrt(s.^2+Par.lambda*PatNum);
	for i=1:4
		W_vec = C./(SIGMA.^(1/p)+eps);
		s1 = solve_Lp_w(s, W_vec, p, Par.J, Par.lambda );
		SIGMA = s1;
	end
	SigmaX = diag(s1);
	Z = U*SigmaX*V';

    % update the multiplier A, fix Z and X
    A = A + Par.rho * (X - Z);
    Par.rho = min(1e4, Par.mu * Par.rho);
end
sigma = sqrt( mean( reshape(mean((Y - Z).^2, 2)', [Par.ps2 Par.ch])) )';
return;
