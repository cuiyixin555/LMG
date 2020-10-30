function latent = deconv_RL_sat(im,psf,reg_strength)

    blurred = im;
%     blurred = imread(im_path);
%     blurred = im2double(blurred);
%     psf = imread(psf_path);
    psf = double(psf);
    if size(psf,3) == 3
        psf = mean(psf,3);
    end
    psf = psf / sum(psf(:));
    
    if ~exist('reg_strength','var')
        reg_strength = 0.004;
    end
    
    blurred_padded = wrap_boundary_liu(blurred, [size(blurred,1) size(blurred,2)] + size(psf) - 1);
    
    % add outlier
%     blurred_padded(100:103, 100:103,:) = 1;
    
    % mask
    mask = zeros(size(blurred_padded));
    mask(1:size(blurred,1), 1:size(blurred,2), :) = 1;

    % set options for deconvolution
    opt.mask = mask;
    opt.outlier_handling = 'cho';  % 'none', 'cho', 'whyte'
    opt.reg_strength = reg_strength;
    opt.sps_prior_exponent = 0.5; % 1: total variation,  less than 1: sparse prior used by Levin et al.
    opt.prevent_ringing = true; % Whyte et al.'s trick to prevent ringing artifacts from saturated pixels
    latent = deconv_lucy_outlier(blurred_padded, psf, opt);
    
    latent = latent(1:size(blurred,1), 1:size(blurred,2), :);

%     figure, imshow(latent);
end


function latent = deconv_lucy_outlier(blurred, psf, options)

    if ~exist('options', 'var')
        options = struct([]);
    end
    
    if ~isfield(options, 'mask')
        options.mask = 1;
    end
    
    if ~isfield(options, 'n_iters')
        options.n_iters = 50;
    end
    
    if ~isfield(options, 'reg_strength')
        options.reg_strength = 0.001;
    end
    
    if ~isfield(options, 'sps_prior_exponent')
        options.sps_prior_exponent = 1;
    end
    
    if ~isfield(options, 'outlier_handling')
        options.outlier_handling = 'none';
    end
    
    if ~isfield(options, 'prevent_ringing')
        options.prevent_ringing = false;
    end
    
    psf_f = psf2otf(psf, size(blurred));

    latent = blurred;

    if options.prevent_ringing
        dilate_radius = 3;
        dilate_filter = bsxfun(@plus,(-dilate_radius:dilate_radius).^2,(-dilate_radius:dilate_radius)'.^2) <= eps+dilate_radius.^2;
        smooth_filter = fspecial('gaussian', [21 21], 3);
        psf_dilate_f = psf2otf(double(psf~=0), size(blurred));
    end
    
    if strcmp(options.outlier_handling, 'cho')
        P_in = 0.9;
    end
    
    if options.reg_strength ~= 0
        dxf = [0 -1 1];
        dyf = dxf';
    end

    for iter=1:options.n_iters
        b = real(ifft2(fft2(latent) .* psf_f));

        % outlier weight
        if strcmp(options.outlier_handling, 'cho') && iter > 1
            err = b - blurred;
            sigma2 = max(latent, 1/255);
            err_dist = 1./sqrt(2*pi*sigma2).*exp(-err.^2./(2*sigma2));
            weight = err_dist.*P_in ./ (err_dist.*P_in + (1-P_in));
            weight(b>1) = 0;

            weight = weight .* options.mask;
        elseif strcmp(options.outlier_handling, 'whyte')
            a = 50;
            b_sat = b - 1/a * log(1 + exp(a*(b - 1)));
            b_sat_p = 1 - 1/a * (a*exp(a*(b-1)))./(1+exp(a*(b-1)));
            b = b_sat;
            weight = options.mask .* b_sat_p;
        else
            weight = options.mask;
        end
        
        if options.prevent_ringing
            S_mask = double(imdilate(latent >= 0.9, dilate_filter));
            V_mask = 1 - min(real(ifft2(fft2(S_mask).*psf_dilate_f)), 1);

            uu = blurred ./ max(b, eps) .* weight - weight;
            ratio_u = uu .* V_mask + 1;
            ratio_s = uu + 1;

            ratio_u = real(ifft2(fft2(ratio_u) .* conj(psf_f)));
            ratio_s = real(ifft2(fft2(ratio_s) .* conj(psf_f)));

            blend_weights = imfilter(S_mask, smooth_filter);
            update_ratio = ratio_u + (ratio_s - ratio_u) .* blend_weights;
        else
            update_ratio = blurred ./ max(b, eps) .* weight + 1 - weight;
            update_ratio = real(ifft2(fft2(update_ratio) .* conj(psf_f)));
        end
            
        if options.reg_strength ~= 0
            alpha = options.sps_prior_exponent;
            dx  = imfilter(latent,dxf,'same','circular');
            dy  = imfilter(latent,dyf,'same','circular');
            threshold = 0.01;

            sps_x = zeros(size(latent));
            m = dx >= threshold;
            sps_x(m) = alpha*dx(m).^(alpha-1);
            m = dx <= -threshold;
            sps_x(m) = -alpha*(-dx(m)).^(alpha-1);
            m = dx < threshold & dx > -threshold;
            sps_x(m) = alpha * threshold^(alpha-2) * dx(m);
            
            sps_y = zeros(size(latent));
            m = dy >=  threshold;
            sps_y(m) = alpha*dy(m).^(alpha-1);
            m = dy <= -threshold;
            sps_y(m) = -alpha*(-dy(m)).^(alpha-1);
            m = dy < threshold & dy > -threshold;
            sps_y(m) = alpha * threshold^(alpha-2) * dy(m);
            
            sps = imfilter(sps_x,dxf,'same','circular','conv') + imfilter(sps_y,dyf,'same','circular','conv');
            update_ratio = update_ratio ./ (1 + options.reg_strength * sps);
        end
        
        update_ratio(update_ratio < 0) = 0;
        
        latent = latent .* update_ratio;
        latent(latent<0) = 0;

%         figure(5), imshow(latent), title(sprintf('%d', iter));
        drawnow;
    end

%     figure(2), imshow(blurred);
end


function latent = deconv_gaussian(blurred, psf, alpha)
    lap = [0 -1 0;-1 4 -1;0 -1 0];
    R = psf2otf(lap, size(blurred));
    B = fft2(blurred);
    K = psf2otf(psf, size(blurred));
    L = conj(K).*B ./ (conj(K).*K + alpha.*R);
    latent = real(ifft2(L));
end
