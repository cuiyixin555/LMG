function [k, lambda_lmg, lambda_grad, S] = blind_deconv_main(blur_B, k, ...
                                    lambda_lmg, lambda_grad, threshold, opts)
% derivative filters
dx = [-1 1; 0 0];
dy = [-1 0; 1 0];
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2013-08-11
H = size(blur_B,1);    W = size(blur_B,2);
blur_B_w = wrap_boundary_liu(blur_B, opt_fft_size([H W]+size(k)-1));
blur_B_tmp = blur_B_w(1:H,1:W,:);
Bx = conv2(blur_B_tmp, dx, 'valid');
By = conv2(blur_B_tmp, dy, 'valid');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for iter = 1:opts.xk_iter
    %%
   if lambda_lmg==0
       %% L0 deblurring
       S = L0Restoration(blur_B, k, lambda_grad, 2.0);
   else
       S = L0_LMG_deblur(blur_B_w, k, lambda_lmg, lambda_grad, 2.0);
       S = S(1:H,1:W,:);
   end
   %% Necessary for refining gradient ???
  [latent_x, latent_y, threshold]= threshold_pxpy_v1(S,max(size(k)),threshold); 
  %
  k = estimate_psf(Bx, By, latent_x, latent_y, 2, size(k));
  %%
  fprintf('pruning isolated noise in kernel...\n');
  CC = bwconncomp(k,8);
  for ii=1:CC.NumObjects
      currsum=sum(k(CC.PixelIdxList{ii}));
      if currsum<.1 
          k(CC.PixelIdxList{ii}) = 0;
      end
  end
  k(k<0) = 0;
  k=k/sum(k(:));
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%
  %% Parameter updating
  if lambda_lmg~=0;
      lambda_lmg = max(lambda_lmg/1.1, 1e-4);
  else
      lambda_lmg = 0;
  end
%   lambda_lmg = lambda_lmg/1.1;  %% for natural images
  if lambda_grad~=0;
      lambda_grad = max(lambda_grad/1.1, 1e-4);
  else
      lambda_grad = 0;
  end
  %
  figure(1); 
  S(S<0) = 0;
  S(S>1) = 1;
  subplot(1,3,1); imshow(blur_B,[]); title('Blurred image');
  subplot(1,3,2); imshow(S,[]);title('Interim latent image');
  subplot(1,3,3); imshow(k,[]);title('Estimated kernel');
  drawnow;
end;
k(k<0) = 0;  
k = k ./ sum(k(:));  
end
