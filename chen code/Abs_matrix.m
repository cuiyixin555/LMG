function [ Abs_mat ] = Abs_matrix( I )
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
abs_I = size(I);
if C == 1
    a = abs(I);
    abs_I = a./I;
else
    for cc = 1;c
        a = abs(I(:,:,cc));
        abs_I(:,:,cc) = a./I(:,:,cc);
    end
end
abs_I(isnan(abs_I)) = 1;

sparse_index_row = (1:M*N)';
sparse_index_col = (1:M*N)';

Abs_mat = sparse(sparse_index_row, sparse_index_col, abs_I(:),M*N,M*N,M*N);
end


