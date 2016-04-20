%start script
load('../share/local_uniform_2d.mat')
pkg load image %statistics
tic
gwr(Data)
toc
