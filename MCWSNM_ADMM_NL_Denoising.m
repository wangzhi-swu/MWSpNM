function [E_Img, Par]   =  MCWSNM_ADMM_NL_Denoising( N_Img, O_Img, Par )
E_Img           = N_Img;   % Estimated Image
[h, w, ch]  = size(E_Img);
Par.h = h;
Par.w = w;
Par.ch = ch;
Par = SearchNeighborIndex( Par );
% noisy image to patch
NoiPat =	Image2Patch( N_Img, Par );
Par.TolN = size(NoiPat, 2);
Par.Sigma = sqrt(mean(Par.nSig0.^2));
Sigma = ones(Par.ch, length(Par.SelfIndex));
for iter = 1 : Par.Iter
    Par.iter = iter;
    % iterative regularization
    E_Img =	E_Img + Par.delta * (N_Img - E_Img);
    % image to patch
    CurPat =	Image2Patch( E_Img, Par );
    % estimate local noise variance
    for c = 1:Par.ch
        if(iter == 1)
            TempSigmaRow = Par.lambda * Par.nSig0(c) * Sigma(c, :);
        else
            TempSigmaRow = Par.lambda * Sigma(c, :);
        end
        SigmaRow((c-1)*Par.ps2+1:c*Par.ps2, :) = repmat( TempSigmaRow , [Par.ps2, 1]);
    end
    SigmaCol = Par.lambda2 * sqrt(abs(repmat(Par.Sigma^2, 1, size(CurPat,2)) - mean((NoiPat - CurPat).^2))); % Estimated Local Noise Level
    if (mod(iter-1, Par.Innerloop) == 0)
        Par.nlsp = Par.nlsp - 10;  % Lower Noise level, less NL patches
        NL_mat  =  Block_Matching_Real(CurPat, Par);% Caculate Non-local similar patches for each
%         if iter == 1 && Par.lambda ~= 0
%             SigmaRow = ones(Par.ps2ch, length(Par.SelfIndex));
%         end
%     if iter == 1 && Par.lambda2~=0
%         SigmaCol = Par.Sigma * ones(size(SigmaCol));
%     end
    end
    [Y_hat, W_hat, Sigma]  =  MCWSNM_ADMM_NL_Estimation( NL_mat, SigmaRow,SigmaCol, NoiPat, Par );  % Estimate all the patches
    E_Img = PGs2Image(Y_hat, W_hat, Par);
    E_Img(E_Img>255)  =  255;
    E_Img(E_Img<0.0)  =  0.0;
    %%
    PSNR  =  csnr( O_Img, E_Img, 0, 0 );
    SSIM  =  cal_ssim( O_Img, E_Img, 0, 0 );
    fprintf( 'Iter = %2.0f, PSNR = %2f, SSIM = %2f \n', iter, PSNR, SSIM );
    Par.PSNR(iter, Par.image)   =   PSNR;
    Par.SSIM(iter, Par.image)   =   SSIM;
end
return;





