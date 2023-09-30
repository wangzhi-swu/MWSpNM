clear;

GT_Original_image_dir = 'Real_ccnoise_denoised_part/';
GT_fpath   = fullfile(GT_Original_image_dir, '*mean.png');
TT_Original_image_dir = 'Real_ccnoise_denoised_part/';
TT_fpath   = fullfile(TT_Original_image_dir, '*real.png');
GT_im_dir  = dir(GT_fpath);
TT_im_dir  = dir(TT_fpath);
im_num     = length(TT_im_dir);
% record all the results in each iteration
Par.PSNR = zeros(2, im_num, 'single');
Par.SSIM = zeros(2, im_num, 'single');
for i = [3]
    Par.I = double( imread(fullfile(GT_Original_image_dir, GT_im_dir(i).name)) );
    Par = ParSet(Par.I,i);
    S = regexp(GT_im_dir(i).name, '\.', 'split');
    [h, w, ch] = size(Par.I);
    Par.nim = double( imread(fullfile(TT_Original_image_dir, TT_im_dir(i).name)) );
    for c = 1:ch
        Par.nSig0(c,1) = NoiseEstimation(Par.nim(:, :, c), Par.ps);
    end
    fprintf('The noise levels are %2f, %2f, %2f. \n', Par.nSig0(1), Par.nSig0(2), Par.nSig0(3) );
    PSNR = csnr( Par.nim, Par.I, 0, 0 );
    SSIM = cal_ssim( Par.nim, Par.I, 0, 0 );
    fprintf('The initial value of PSNR = %2f, SSIM = %2f \n', PSNR,SSIM);
    [im_out, Par] = MCWSNM_ADMM_NL_Denoising( Par.nim, Par.I, Par ); 
    im_out(im_out>255) = 255;
    im_out(im_out<0.0) = 0.0;
    % calculate the PSNR
    Par.PSNR(Par.Iter, Par.image)  =  csnr( im_out, Par.I, 0, 0 );
    Par.SSIM(Par.Iter, Par.image)  =  cal_ssim( im_out, Par.I, 0, 0 );
    % imname = sprintf(['../../cc_Results/' Par.method '_NL_CC15_' Par.model '_Oite' num2str(Par.Iter) '_Iite' num2str(Par.maxIter) '_rho' num2str(rho) '_mu' num2str(Par.mu) '_lambda' num2str(lambda) '_' TT_im_dir(i).name]);
    % imwrite(im_out/255, imname);
    fprintf('%s : PSNR = %2f, SSIM = %2f \n',TT_im_dir(i).name, max(Par.PSNR(:, Par.image)),max(Par.SSIM(:, Par.image)) );
end
fprintf("============================================\n")
Par.PSNR(Par.Iter+1,:) = max(Par.PSNR);
Par.SSIM(Par.Iter+1,:) = max(Par.SSIM);
Par.PSNR(Par.Iter+1,im_num+1) = mean(Par.PSNR(Par.Iter+1,:));
Par.SSIM(Par.Iter+1,im_num+1) = mean(Par.SSIM(Par.Iter+1,:));
fprintf(' PSNR = %2.4f, SSIM = %2.4f. \n', Par.PSNR(Par.Iter+1,im_num+1),Par.SSIM(Par.Iter+1,im_num+1));
% results = [Par.PSNR(3,:)',Par.SSIM(3,:)'];
fprintf("============================================\n")
%             mPSNR    = mean(Par.PSNR,2);
%             [~, idx] = max(mPSNR);
%             PSNR = Par.PSNR(idx,:);
%             SSIM = Par.SSIM(idx,:);
%             mSSIM=mean(SSIM,2);
%             fprintf('The best PSNR result is at %d iteration. \n',idx);
%             fprintf('The average PSNR = %2.4f, SSIM = %2.4f. \n', mPSNR(idx),mSSIM);
%             name = sprintf(['../../' Par.method '_NL_CC' num2str(im_num) '_' Par.model '_Oite' num2str(Par.Iter) '_Iite' num2str(Par.maxIter) '_rho' num2str(rho) '_mu' num2str(Par.mu) '_lambda' num2str(lambda) '.mat']);
%             save(name,'PSNR','SSIM','mPSNR','mSSIM');
%             imwrite(im_out/255, imname);
%%
Par.PSNR(Par.Iter+1,:) = max(Par.PSNR);
Par.SSIM(Par.Iter+1,:) = max(Par.SSIM);
Par.PSNR(Par.Iter+1,im_num+1) = mean(Par.PSNR(Par.Iter+1,:));
Par.SSIM(Par.Iter+1,im_num+1) = mean(Par.SSIM(Par.Iter+1,:));
fprintf(' PSNR = %2.4f, SSIM = %2.4f. \n', Par.PSNR(Par.Iter+1,im_num+1),Par.SSIM(Par.Iter+1,im_num+1));
results = [Par.PSNR(3,:)',Par.SSIM(3,:)'];

