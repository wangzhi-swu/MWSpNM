function [ Y_hat, W_hat, Sigma] = MCWSNM_ADMM_NL_Estimation(  NL_mat, SigmaRow,SigmaCol, CurPat, Par )

Y_hat = zeros(size(CurPat));
W_hat = zeros(size(CurPat));
Sigma = zeros(Par.ch, length(Par.SelfIndex));
for i    =   1 : length(Par.SelfIndex) % For each keypatch group
    Y    =   CurPat(:, NL_mat(1:Par.nlsp,i)); % Non-local similar patches to the keypatch
    mY   =   repmat(mean( Y, 2 ),1,Par.nlsp);
    Y    =   Y-mY;
    [X, Sigma(:, i)] 	=   MCWSNM_ADMM_NL(  Y, SigmaRow(:, i),SigmaCol(:, NL_mat(1:Par.nlsp,i)), Par); % WSNM Estimation
    Y_hat(:,NL_mat(1:Par.nlsp,i))  =  Y_hat(:,NL_mat(1:Par.nlsp,i))+X+mY;
    W_hat(:,NL_mat(1:Par.nlsp,i))  =  W_hat(:,NL_mat(1:Par.nlsp,i))+ones(Par.ps2ch, Par.nlsp);
end
end

