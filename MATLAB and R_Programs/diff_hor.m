function L1=diff_hor(row,col)
%function to compute 1st order difference along the horizontal
L=diag(ones(col,1)) - diag(ones(col-1,1),-1);
% L=zeros(col,col);
% L(1,1)=1;
% for i=2:col
%     L(i,i-1:i)=[-1,1];
% end
nn=tril(repmat(L,row));
 L1=triu(nn,-1);