# function to calculate the first order difference matrix horizontally
# row is the number of the pixel along the vertical direction
# col is the number of pixel along the horizontal durection
# the total pixel will be row*col
diff_hor=function(row,col){
  c=rep(0,row*col);c[1]=1;c[2]=-1
L=toeplitz(c);L[upper.tri(L)]=0;b=1;cc=0
  for (i in 2:row){L[col+b,col+cc]=0;b=col+b;cc=col+cc}
  return(list(L))}