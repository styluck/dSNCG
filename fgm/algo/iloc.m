function [dim1, dim2] =iloc(i, ro)
dim1 = mod(i, ro);
if dim1 == 0
    dim1 = ro;
    dim2 = floor(i/ro);
else
    dim2 = floor(i/ro)+1;
end
