function b=idst(a,n)
% IDST Inverse discrete sine transform (Used in Poisson reconstruction)
%
%    X = IDST(Y) inverts the DST transform, returning the
%    original vector if Y was obtained using Y = DST(X).
%    X = IDST(Y,N) pads or truncates the vector Y to length N
%    before transforming.
%    If Y is a matrix, the IDST operation is applied to each column.

if nargin == 1
    if min(size(a)) == 1
        n = length(a);
    else
        n = size(a,1);
    end
end

n2 = n + 1;            
b = 2/n2*dst(a,n);
end
