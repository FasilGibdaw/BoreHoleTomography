function L2=diff_ver(row,col)
% function to compute 1st order difference along the
%vertical, Note that: pixels are assigned column wise
c=zeros(1,row*col);
c(1)=1;c(col+1)=-1;
L2=tril(toeplitz(c));