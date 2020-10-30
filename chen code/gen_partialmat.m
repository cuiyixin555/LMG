function [px_mat,py_mat] = gen_partialmat(im_row,im_col)

rowpx_seq = [];
colpx_seq = [];
valuepx_seq = [];
rowpy_seq = [];
colpy_seq = [];
valuepy_seq = [];

% rowpxx_seq = [];
% colpxx_seq = [];
% valuepxx_seq = [];
% rowpxy_seq = [];
% colpxy_seq = [];
% valuepxy_seq = [];
% rowpyy_seq = [];
% colpyy_seq = [];
% valuepyy_seq = [];
for j=1:im_col
    for i=1:im_row
        ind = (i+(j-1)*im_row);

        if(i==1)
            rowpy_seq = [rowpy_seq,ind,ind];
            colpy_seq = [colpy_seq,ind,ind+1];
            valuepy_seq = [valuepy_seq,-1,1];        
        else
            rowpy_seq = [rowpy_seq,ind,ind];
            colpy_seq = [colpy_seq,ind-1,ind];
            valuepy_seq = [valuepy_seq,-1,1];
        end
        
        if(j==1)
            rowpx_seq = [rowpx_seq,ind,ind];
            colpx_seq = [colpx_seq,ind,ind+im_row];
            valuepx_seq = [valuepx_seq,-1,1];       
        else
            rowpx_seq = [rowpx_seq,ind,ind];
            colpx_seq = [colpx_seq,ind,ind-im_row];
            valuepx_seq = [valuepx_seq,1,-1];
        end
        
%         if(i==1)
%             rowpyy_seq = [rowpyy_seq,ind,ind,ind];
%             colpyy_seq = [colpyy_seq,ind,ind+1,ind+2];
%             valuepyy_seq = [valuepyy_seq,1,-2,1];        
%         elseif(i==2)
%             rowpyy_seq = [rowpyy_seq,ind,ind,ind];
%             colpyy_seq = [colpyy_seq,ind-1,ind,ind+1];
%             valuepyy_seq = [valuepyy_seq,1,-2,1];             
%         else
%             rowpyy_seq = [rowpyy_seq,ind,ind,ind];
%             colpyy_seq = [colpyy_seq,ind-2,ind-1,ind];
%             valuepyy_seq = [valuepyy_seq,1,-2,1];
%         end
% 
%         if(j==1)
%             rowpxx_seq = [rowpxx_seq,ind,ind,ind];
%             colpxx_seq = [colpxx_seq,ind,ind+im_row,ind+2*im_row];
%             valuepxx_seq = [valuepxx_seq,1,-2,1];      
%         elseif(j==2)
%             rowpxx_seq = [rowpxx_seq,ind,ind,ind];
%             colpxx_seq = [colpxx_seq,ind-im_row,ind,ind+im_row];
%             valuepxx_seq = [valuepxx_seq,1,-2,1];             
%         else
%             rowpxx_seq = [rowpxx_seq,ind,ind,ind];
%             colpxx_seq = [colpxx_seq,ind-2*im_row,ind-im_row,ind];
%             valuepxx_seq = [valuepxx_seq,1,-2,1];
%         end
        
%         if(i==1 && j==1)
%             rowpxy_seq = [rowpxy_seq,ind,ind,ind,ind];
%             colpxy_seq = [colpxy_seq,ind,ind+1,ind+im_row,ind+im_row+1];
%             valuepxy_seq = [valuepxy_seq,1,-1,-1,1];
%         elseif(i==1)
%             rowpxy_seq = [rowpxy_seq,ind,ind,ind,ind];
%             colpxy_seq = [colpxy_seq,ind-im_row,ind-im_row+1,ind,ind+1];
%             valuepxy_seq = [valuepxy_seq,1,-1,-1,1];            
%         elseif(j==1)
%             rowpxy_seq = [rowpxy_seq,ind,ind,ind,ind];
%             colpxy_seq = [colpxy_seq,ind-1,ind,ind-1+im_row,ind+im_row];
%             valuepxy_seq = [valuepxy_seq,1,-1,-1,1];     
%         else
%             rowpxy_seq = [rowpxy_seq,ind,ind,ind,ind];
%             colpxy_seq = [colpxy_seq,ind-1-im_row,ind-im_row,ind-1,ind];
%             valuepxy_seq = [valuepxy_seq,1,-1,-1,1];               
%         end
        
    end
end

px_mat = sparse(rowpx_seq,colpx_seq,valuepx_seq);
py_mat = sparse(rowpy_seq,colpy_seq,valuepy_seq);
% pxx_mat = sparse(rowpxx_seq,colpxx_seq,valuepxx_seq);
% pxy_mat = sparse(rowpxy_seq,colpxy_seq,valuepxy_seq);
% pyy_mat = sparse(rowpyy_seq,colpyy_seq,valuepyy_seq);
end