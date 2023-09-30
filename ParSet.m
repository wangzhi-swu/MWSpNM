function Par = ParSet(imOri,i)
[Par.m,Par.n,Par.p] = size(imOri);
Par.image = i;
Par.win       =   40;   % Non-local patch searching window
Par.delta     =   0;  % Parameter between each iter
Par.Constant  =   2 * sqrt(2);   % Constant num for the weight vector
Par.Innerloop =   2;   % InnerLoop Num of between re-blockmatching
Par.ps        =   7;   % Patch size
Par.step      =   4;
Par.nlsp  =  70;      
Par.Iter      =   2; % total iter numbers
Par.p = 0.9;
Par.J = 4; 
Par.display = true;
% Par.method = 'WNNM_ADMM';
Par.method  = 'MCWSNM_ADMM';
Par.model   = '2';
Par.maxIter = 10;
Par.model   = 1;
Par.I = imOri;

% this parameter is not finally determined yet
Par.mu = 1+eps;
Par.rho = 0.5;
Par.lambda = 23.5;
Par.lambda2 = 11;
Par.p2 = 0.9;
