function [ output_img,A ] = LMG( img,patch_size )
%% Codes provided by Chen, Liang. Please cite the following paper
%% if it is used
% @inproceedings{chen2019blind,
%   title={Blind Image Deblurring With Local Maximum Gradient Prior},
%   author={Chen, Liang and Fang, Faming and Wang, Tingting and Zhang, Guixu},
%   booktitle={Proceedings of the IEEE Conference on Computer Vision and Pattern Recognition},
%   pages={1742--1750},
%   year={2019}
% }
%UNTITLED8 此处显示有关此函数的摘要
%   此处显示详细说明
[M,N] = size(img);
[px_mat,py_mat] = gen_partialmat(M,N); 
px = reshape(px_mat*img(:),[M,N]);
py = reshape(py_mat*img(:),[M,N]);
abs_x_mat = Abs_matrix(px);
abs_y_mat = Abs_matrix(py);

max_tv_mat = Max_matrix(abs(px) + abs(py),patch_size);
output_img =  max_tv_mat*(abs_x_mat*px_mat + abs_y_mat*py_mat)*img(:);
output_img = reshape(output_img,[M,N]);
A = max_tv_mat*(abs_x_mat*px_mat + abs_y_mat*py_mat);
end

