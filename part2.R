
#test1: Basic (trivial, unoprtimized) matrix multiply code 
#basic_dgemm.c
dyn.load("basic_dgemm.so")
set.seed(123)
A = replicate(10, rnorm(10)) 
set.seed(234)
B = replicate(10, rnorm(10))
C = matrix(as.double(0), nrow = nrow(A), ncol = ncol(A))

matmatmult1 <- function(m,a, b, c)
{
  .C("square_dgemm", m=m, a = as.double(a), b = as.double(b),
     c = as.double(c))$c
}

start.time <- Sys.time()
result1 <- matmatmult1(100,A,B,C)
end.time <- Sys.time()
time_ls[i] = start.time-end.time
#test2: A very simple blocked implementation
#blocked_dgemm.c  with default BLOCK_SIZE ((unsigned) 16)
dyn.load("blocked_dgemm.so")


#test3: A wrapper for calling BLAS library
#wrap_dgemm.c
dyn.load("wrap_dgemm.so")
matmatmult1 <- function(m,a, b, c)
{
  .C("square_dgemm", m=m, a = as.double(a), b = as.double(b),
     c = as.double(c))$c
}

#plot
library(ggplot2)
library(reshape2)
#1
flops = read.csv('basic_block_compare.csv',header=F)/10^6
sizes = seq(100,1500,100)
perf.dt = data.frame(matrix_size=sizes,basic_mul=flops[,1],block_mul=flops[,2])
nperf.dt = melt(perf.dt,id='matrix_size')
ggplot(data=nperf.dt,aes(x=matrix_size,y=value,colour=variable))+
  geom_line()+scale_y_continuous('Performance (Mflops/s)')+
  ggtitle("Basic V.S Block Matrix-Matrix Multiplication")

#2
blocks = read.csv('block_size_compare.csv',header=F)
blocks[,3] = blocks[,3]/10^6
names(blocks) = c('BlockSize','MatrixSize','Performance')
ggplot(blocks, aes(x = MatrixSize, y = Performance)) +
  geom_line(aes(color = factor(BlockSize)))+scale_y_continuous('Performance (Mflops/s)')+
  ggtitle("Different Block Matrix-Matrix Multiplications")
