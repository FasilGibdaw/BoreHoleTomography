# function to calculate the first order difference matrix vertically
# row is the number of the pixel along the vertical direction
# col is the number of pixel along the horizontal durection
# the total pixel will be row*col
diff_ver=function(row,col){
  c=rep(0,row*col);c[1]=1;c[col+1]=-1
  L=toeplitz(c);L[upper.tri(L)]=0
  return(list(L))}