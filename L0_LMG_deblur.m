function S = L0_LMG_deblur(Im, kernel, lambda, wei_grad, kappa)
if ~exist('kappa','var')
    kappa = 2.0;
end
S = Im;
betamax = 1e5;
fx = [1, -1];
fy = [1; -1];
[N,M,D] = size(Im);
sizeI2D = [N,M];
otfFx = psf2otf(fx,sizeI2D);
otfFy = psf2otf(fy,sizeI2D);
%%
KER = psf2otf(kernel,sizeI2D);
Den_KER = abs(KER).^2;
%%
Denormin2 = abs(otfFx).^2 + abs(otfFy ).^2;
if D>1
    Denormin2 = repmat(Denormin2,[1,1,D]);
    KER = repmat(KER,[1,1,D]);
    Den_KER = repmat(Den_KER,[1,1,D]);
end
Normin1 = conj(KER).*fft2(S);
%% pixel sub-problem
%%
patch_size = 35; %% Fixed size!
mybeta_pixel = lambda/(graythresh((S).^2));
for iter = 1:4
    [J,A] = LMG(S, patch_size);
    t = 2 - J; t2 = lambda/(2*mybeta_pixel);
    t3 = abs(t) - t2; t3(t3<0) = 0;
    u = sign(t).*t3;
    clear t;
    %% Gradient sub-problem
    alpha3 =mybeta_pixel*2; 
    for i = 1:4
        [M,N] = size(A);
        sparse_eye = speye(M,N);
        subsitute_I = (mybeta_pixel*(A'*A) + alpha3*sparse_eye) \ (mybeta_pixel*A'*(2-u(:)) + alpha3*(S(:)));
        subsitute_I = reshape(subsitute_I,size(u));
        beta = 2*wei_grad;
        while beta < betamax
            h = [diff(S,1,2), S(:,1,:) - S(:,end,:)];
            v = [diff(S,1,1); S(1,:,:) - S(end,:,:)];
            th = h.^2<wei_grad/beta;
            tv = v.^2<wei_grad/beta;
            h(th)=0; v(tv)=0;
            clear t;
            Normin2 = [h(:,end,:) - h(:, 1,:), -diff(h,1,2)];
            Normin2 = Normin2 + [v(end,:,:) - v(1, :,:); -diff(v,1,1)];
            FS = (Normin1 + beta*fft2(Normin2) + alpha3*fft2(subsitute_I))./(Den_KER + beta*Denormin2 + alpha3);
            S = real(ifft2(FS));
            beta = beta*kappa;
        end
        alpha3 = alpha3*4;
    end
     mybeta_pixel = mybeta_pixel*4;
end
end

