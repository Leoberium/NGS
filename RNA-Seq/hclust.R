e = matrix(rnorm(2000,0:10/5),ncol=10,byrow = T)
dim(e)
e[1:10,]

c = cor(t(e),m='sp')
dim(c)
h = hclust(as.dist(1-c))
plot(h)
cl = cutree(h,k = 4)
cl
table(cl)

z = t(scale(t(e)))
plot(apply(z[cl==1,],2,mean))
par(mfrow=c(2,2))
for(i in 1:4)
boxplot(z[cl==i,])
