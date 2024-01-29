library(huge) # Load the package huge

L = huge.generator(n=200,d=200,graph="hub") # Generate data with hub structures

Sigma.hat <- cor(X<-L$data)
image(Sigma.hat) # Each block on the plot corresponds to a hub

X = L$data; X.pow = X^3/sqrt(15) # Power Transformation
X.npn = huge.npn(X.pow) # Non paranormal - we transform back to normal

out.mb = huge(X.pow,nlambda=30) # Estimate the solution path
out.npn = huge(X.npn,nlambda=30)
huge.roc(out.mb$path,L$theta) # Plot the ROC curve
huge.roc(out.npn$path,L$theta)
mb.stars = huge.select(out.mb,criterion="stars", stars.thresh=0.05) # Select the graph using StARS
npn.stars = huge.select(out.npn,criterion="stars",stars.thresh=0.05)
mb.ric = huge.select(out.mb) # Select the graph using RIC
npn.ric = huge.select(out.npn)

plot(mb.stars)
plot(out.npn)

#========================= Let's change number of replicates to 50===============

L = huge.generator(n=50,d=200,graph="hub") # Generate data with hub structures

Sigma.hat <- cor(X<-L$data)
image(Sigma.hat) # Each block on the plot corresponds to a hub

X = L$data; X.pow = X^3/sqrt(15) # Power Transformation
X.npn = huge.npn(X.pow) # Non paranormal - we transform back to normal

out.mb = huge(X.pow,nlambda=30) # Estimate the solution path
out.npn = huge(X.npn,nlambda=30)
huge.roc(out.mb$path,L$theta) # Plot the ROC curve
huge.roc(out.npn$path,L$theta)
mb.stars = huge.select(out.mb,criterion="stars", stars.thresh=0.05) # Select the graph using StARS
npn.stars = huge.select(out.npn,criterion="stars",stars.thresh=0.05)
mb.ric = huge.select(out.mb) # Select the graph using RIC
npn.ric = huge.select(out.npn)