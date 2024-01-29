
library(huge)
L<-huge.generator(n=20, d=11, graph="random", prob =0.2)

X<-L$data
Rho<-cor(X) # correlation matrix of the data
plot(L) # True graph

L10<-huge.generator(n=10, d=11, graph="random", prob =0.2)
L100<-huge.generator(n=100, d=11, graph="random", prob =0.2)
L1000<-huge.generator(n=1000, d=11, graph="random", prob =0.2)

X10<-L10$data
Rho<-cor(X10) # correlation matrix of the data
plot(L10) # True graph

X100<-L100$data
Rho<-cor(X100) # correlation matrix of the data
plot(L100) # True graph

X1000<-L1000$data
Rho<-cor(X1000) # correlation matrix of the data
plot(L1000) # True graph