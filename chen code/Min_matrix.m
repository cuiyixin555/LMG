function [ min_mat ] = Min_matrix( I,patch_size )
%% Codes provided by Chen, Liang. Please cite the following paper
%% if it is used
% @inproceedings{chen2019blind,
%   title={Blind Image Deblurring With Local Maximum Gradient Prior},
%   author={Chen, Liang and Fang, Faming and Wang, Tingting and Zhang, Guixu},
%   booktitle={Proceedings of the IEEE Conference on Computer Vision and Pattern Recognition},
%   pages={1742--1750},
%   year={2019}
% }

%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
[M, N, C] = size(I);
J_index = zeros(M, N); % Create empty index matrix

% Test if patch size has odd number
if ~mod(numel(patch_size),2) % if even number
    error('Invalid Patch Size: Only odd number sized patch supported.');
end
% I = padarray(I, [floor(patch_size./2) floor(patch_size./2)], 'replicate');
padsize = floor(patch_size/2);

h = ceil(patch_size/2);
for m = 1:M
        for n = 1:N
            patch = I(max(1,m-padsize):min(M,m+padsize),max(1,n-padsize):min(N,n+padsize));
            [h1,h2] = size(patch);
            tmp = max(patch, [], 3);
            [~, tmp_idx] = max(tmp(:));
            ori_i = h - (patch_size-h1);
            ori_j = h - (patch_size-h2);
            if ori_i ~= h && m > h
                ori_i = h1+1-ori_i;
            end
            if ori_j ~= h && n > h
                ori_j = h2+1-ori_j;
            end
            J_need = ceil(tmp_idx/h1);
            I_need = tmp_idx - (J_need-1)*h1;
            
            i_quote = m + I_need - ori_i;
            j_quote = n + J_need - ori_j;
            
            J_index(m,n) = (j_quote-1)*M + i_quote;
        end
end
sparse_index_row = (1:M*N)';
ss = ones(M,N);
ss = ss(:);
sparse_index_col = J_index(:);
min_mat = sparse(sparse_index_row, sparse_index_col,ss, M*N, M*N, M*N);
end


