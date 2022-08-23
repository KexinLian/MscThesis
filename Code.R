# LW
library(ggplot2)
library(grid)
library(aplore3)
library(clusterPower)
expit = function(x){
  return(1/(1+1/exp(x)))
}
library(R2jags)

# the combination function in the ggplot2
vp.layout <- function(x, y) viewport(layout.pos.row=x, layout.pos.col=y)
arrange <- function(..., nrow=NULL, ncol=NULL, as.table=FALSE) {
  dots <- list(...)
  n <- length(dots)
  if(is.null(nrow) & is.null(ncol)) { nrow = floor(n/2) ; ncol = ceiling(n/nrow)}
  if(is.null(nrow)) { nrow = ceiling(n/ncol)}
  if(is.null(ncol)) { ncol = ceiling(n/nrow)}
  ## NOTE see n2mfrow in grDevices for possible alternative
  grid.newpage()
  pushViewport(viewport(layout=grid.layout(nrow,ncol) ) )
  ii.p <- 1
  for(ii.row in seq(1, nrow)){
    ii.table.row <- ii.row
    if(as.table) {ii.table.row <- nrow - ii.table.row + 1}
    for(ii.col in seq(1, ncol)){
      ii.table <- ii.p
      if(ii.p > n) break
      print(dots[[ii.table]], vp=vp.layout(ii.table.row, ii.col))
      ii.p <- ii.p + 1
    }
  }
}

# uniform priors on logit scale
random_uniform = runif(5000,0,1)
data_random_uniform = data.frame(random_uniform)
p = ggplot(data = data_random_uniform,aes(x = random_uniform))
p1 = p+geom_density()+ggtitle("Uniform priors")
random_uniform_logit = log(random_uniform/(1-random_uniform))
data_random_uniform_logit = data.frame(random_uniform_logit)
p = ggplot(data = data_random_uniform_logit,aes(x = random_uniform_exp))
p2 = p+geom_density()+ggtitle("Uniform priors on logit scale")
arrange(p1,p2)

#random_normal = rnorm(5000,0,10^6)
#data_random_normal = data.frame(random_normal)
#p = ggplot(data = data_random_noraml,aes(x = random_normal))
#p3 = p+geom_density()+ggtitle("Normal priors with large variance 10^6")
#random_normal_logit = log(random_normal/(1-random_normal))
#data_random_normal_logit = data.frame(random_normal_logit[-is.na(random_normal_logit)])
#p = ggplot(data = data_random_normal_logit,aes(x = random_normal_logit[-is.na(random_normal_logit)]))
#p4 = p+geom_density()+ggtitle("Normal priors with large variance  on logit scale")


# 3.1 Working with chdage dataset: use glm, get results 
chdage = aplore3::chdage
chd = numeric(100)
chd[aplore3::chdage$chd=='Yes']=1

alpha = rnorm(10000,0,10)
beta = rnorm(10000,0,10)
pi60.glm = exp(alpha+beta*20)/(1+exp(alpha+beta*20))
age.scale = (aplore3::chdage$age-mean(aplore3::chdage$age))/sd(aplore3::chdage$age)
input.chdage = data.frame(x = age.scale,y = chd)

# glm method
chd.glm = glm(y~x,family=binomial,data = input.chdage)
summary(chd.glm)
#             Estimate Std. Error z value Pr(>|z|)    
#(Intercept)  -0.3868     0.2397  -1.613    0.107    
#x             1.3001     0.2820   4.610 4.02e-06 ***
c_alpha1 = c(-0.3868-1.96*0.2397/sqrt(100),-0.3868+1.96*0.2397/sqrt(100))
c_beta1 = c(1.3001-1.96*0.2820/sqrt(100), 1.3001+1.96*0.2820/sqrt(100))

# 3.2 Induced priors
# 3.2.3 N(0,S^2/b)

# The induced priors for christensen et al.

library(ggplot2)
library(grid)
S_2 = var(age.scale)

# The function of generating induced priors for different b values
# draw the normal distribution of the priors
induced_priors_chris = function(b){
  induced.prior.prob.chdage = matrix(data = NA,nrow = 100, ncol = 1000)
  induced.prior.x      = matrix(data = NA,nrow = 100, ncol = 1000)
  set.seed(1)

  alpha = rnorm(1000,0,sqrt(S_2/b))
  beta = rnorm(1000,0,sqrt(S_2/b))
  for(i in 1:100){
    induced.prior.prob.chdage[i,] = expit(alpha+beta*age.scale[i])
    induced.prior.x[i,] = expit(alpha+beta*x[i])
  }
  hist(induced.prior.prob.chdage,xlab = expression(pi(x)),main = paste("Normal(0,",expression(s^2/b),") priors","b =", b))
  hist(induced.prior.x,xlab = expression(pi(x)), ylab="frequency",main = paste("random normal variables","Normal(0,",expression(s^2/b),") priors","b =", b))
}

par(mfrow = c(5,2))
induced_priors_chris(2.133)
# try different value of b
induced_priors_chris(1.25)
induced_priors_chris(1)
induced_priors_chris(0.75)
induced_priors_chris(0.5)
induced_priors_chris(0.25)


# b = 0.7
b = 0.75
# Compare the induced priors with uniform priors
par(mfrow = c(1,2))
induced.prior.prob.chdage = matrix(data = NA,nrow = 100, ncol = 1000)
induced.prior.x      = matrix(data = NA,nrow = 100, ncol = 1000)
x = rnorm(100,0,1)
set.seed(1)
alpha = rnorm(1000,0,sqrt(S_2/b))
beta = rnorm(1000,0,sqrt(S_2/b))
for(i in 1:100){
  induced.prior.prob.chdage[i,] = expit(alpha+beta*age.scale[i])
  induced.prior.x[i,] = expit(alpha+beta*x[i])
}
hist(induced.prior.prob.chdage,xlab = expression(pi(x)),main = "read dataset")
hist(induced.prior.x,xlab = expression(pi(x)), ylab="frequency",main = "random uniform priors")


# 3.2.1 N(0,25^2)
logit.normal.sigma = function(sigma){
set.seed(1)
induced.prior.normal = matrix(data = NA,nrow = 100, ncol = 1000)
induced.prior.random = matrix(data = NA,nrow = 100, ncol = 1000)
set.seed(1)
alpha = rnorm(1000,0,sigma)
beta = rnorm(1000,0,sigma)
x = rnorm(1000,0,1)
for(i in 1:100){
  induced.prior.normal[i,] = expit(alpha+beta*age.scale[i])
  induced.prior.random[i,] = expit(alpha+beta*x[i])
}
# maybe it is not a good choice
hist(induced.prior.normal,xlab = expression(pi(x)), ylab=("frequency"),main = paste("sigma=",sigma,"from original dataset"))
hist(induced.prior.random,xlab = expression(pi(x)), ylab=("frequency"),main = paste("sigma=",sigma,"from random variable N(0,1)"))
}

par(mfrow = c(1,2))
logit.normal.sigma(25)
logit.normal.sigma(5)
logit.normal.sigma(1)
logit.normal.sigma(0.5)
logit.normal.sigma(0.25)

logit.normal.sigma(2.133)
scale = c(2.5,1,0.5,0.25,0.05,0.01)

# 0.74 0.3
# peak

par(mfrow=c(2,9))

for (i in 1:6){
induced.prior.cauchy = matrix(data = NA,nrow = 100, ncol = 1000)
induced.prior.random = matrix(data = NA,nrow = 100, ncol = 1000)
set.seed(999)
alpha = rcauchy(1000,0,scale[i])
beta = rcauchy(1000,0,scale[i])
x = rnorm(100,0,0.5)
for(i in 1:100){
  induced.prior.cauchy[i,] = expit(alpha+beta*age.scale[i])
  induced.prior.random[i,] = expit(alpha+beta*x[i])
}
hist(induced.prior.cauchy,xlab = expression(pi(x)), ylab=("frequency"),main = ("original data"))
hist(induced.prior.random,xlab = expression(pi(x)), ylab=("frequency"),main = ("random normal variables"))
}

# 3.2.2 Cauchy 
# par(mfrow = c(2,3))
cauchy.induced = function(scale){
set.seed(1)
induced.prior.cauchy = matrix(data = NA,nrow = 100, ncol = 1000)
set.seed(1)
sd = sd(age.scale)
alpha = rcauchy(1000,0,scale)
beta = rcauchy(1000,0,scale)
x = rnorm(100,0,0.5)
for(i in 1:100){
  induced.prior.cauchy[i,] = expit(alpha+beta*age.scale[i])
  induced.prior.random[i,] = expit(alpha+beta*x[i])
}
hist(induced.prior.cauchy,xlab = expression(pi(x)), ylab=("frequency"),main = paste('Cauchy( 0,',scale,') original dataset'))
hist(induced.prior.random,xlab = expression(pi(x)), ylab=("frequency"),main = paste('Cauchy( 0,',scale,') random variable'))
}

set.seed(1)
induced.prior.cauchy = matrix(data = NA,nrow = 100, ncol = 1000)
set.seed(1)
sd = sd(age.scale)
alpha = rcauchy(1000,0,10)
beta = rcauchy(1000,0,2.5/sd)
x = rnorm(100,0,0.5)
for(i in 1:100){
  induced.prior.cauchy[i,] = expit(alpha+beta*age.scale[i])
  induced.prior.random[i,] = expit(alpha+beta*x[i])
}
hist(induced.prior.cauchy,xlab = expression(pi(x)), ylab=("frequency"),main = paste('alpha ~ Cauchy(0,10), beta~ Cauchy(0,2.5/sd) original dataset'))
hist(induced.prior.random,xlab = expression(pi(x)), ylab=("frequency"),main = paste('alpha ~ Cauchy(0,10), beta~ Cauchy(0,2.5/sd) random variable'))

par(mfrow = c(1,2))
sd = sd(age.scale)

cauchy.induced(2.5/sd)
cauchy.induced(2.5)
cauchy.induced(1)
cauchy.induced(0.5)
cauchy.induced(0.25)

# 3.2.3 N(0,S^2/b)
# induced priors
# b = 100

induced_priors_chris(0.7)

# 3.2.4 N(0,pi^2/(3*(p+1)))
par(mfrow = c(1,2))
p = 1
set.seed(1)
induced.prior.pi = matrix(data = NA,nrow = 100, ncol = 1000)
alpha = rnorm(1000,0,sqrt(pi^2/(3*(p+1))))
beta = rnorm(1000,0,sqrt(pi^2/(3*(p+1))))
x = rnorm(100,0,0.5)
for(i in 1:100){
  induced.prior.pi[i,] = expit(alpha+beta*age.scale[i])
  induced.prior.random[i,] = expit(alpha+beta*x[i])
}
hist(induced.prior.pi,xlab = expression(pi(x)), ylab=("frequency"),main = ("N(0,pi^2/3) priors original dataset"))
hist(induced.prior.random,xlab = expression(pi(x)), ylab=("frequency"),main = ("N(0,pi^2/3) priors random variable"))


#3.2.5 CMPs
cmp.induced = function(alpha1,beta1){
pi40 = rbeta(1000,alpha1,beta1)
pi60 = rbeta(1000,alpha1,beta1)
x40 = (40-mean(aplore3::burn1000$age))/sd(aplore3::burn1000$age)
x60 = (60-mean(aplore3::burn1000$age))/sd(aplore3::burn1000$age)
m1 = pi40/(1-pi40)
m2 = pi60/(1-pi60)
beta = (log(m2)-log(m1))/(x60-x40)
alpha = m2 - x60*beta
x = rnorm(100,0,1)
induced.prior.cmp = matrix(data = NA,nrow = 100, ncol = 1000)
induced.prior.random = matrix(data = NA,nrow = 100, ncol = 1000)


for(i in 1:100){
  induced.prior.cmp[i,] = expit(alpha+beta*age.scale[i])
  induced.prior.random[i,] = expit(alpha+beta*x[i])
}

hist(induced.prior.cmp,xlab = expression(pi(x)), ylab=("frequency"),main = paste('beta(',alpha1,',',beta1,')',"CMPs method original dataset"))
hist(induced.prior.random,xlab = expression(pi(x)), ylab=("frequency"),main = paste('beta(',alpha1,',',beta1,')',"CMPs method random variable"))
}
par(mfrow = c(1,2))
cmp.induced(10,10)
cmp.induced(5,5)
cmp.induced(2.5,2.5)
cmp.induced(1,1)
cmp.induced(0.5,0.5)


par(mfrow = c(5,2))
logit.normal.sigma(25)
cauchy.induced(2.5)
induced_priors_chris(0.5)

p = 1
set.seed(1)
induced.prior.pi = matrix(data = NA,nrow = 100, ncol = 1000)
alpha = rnorm(1000,0,sqrt(pi^2/(3*(p+1))))
beta = rnorm(1000,0,sqrt(pi^2/(3*(p+1))))
x = rnorm(100,0,0.5)
for(i in 1:100){
  induced.prior.pi[i,] = expit(alpha+beta*age.scale[i])
  induced.prior.random[i,] = expit(alpha+beta*x[i])
}
hist(induced.prior.pi,xlab = expression(pi(x)), ylab=("frequency"),main = ("N(0,pi^2/3) priors original dataset"))
hist(induced.prior.random,xlab = expression(pi(x)), ylab=("frequency"),main = ("N(0,pi^2/3) priors random variable"))

cmp.induced(1,1)



# 3.3 Bayesian Analysis
# 3.3.1 Using jags and dnorm() priors on alpha and beta and should get values close to the glm
# ***************** Bayesian Inference *******************
num.chains   <- 3
adapt.length <- 1000
burnin       <- 1000
chain.length <- 5000

#num.iteration = 1000+5000

#-- Wide Normal Prior Case ---------------------
beta.sd <- 25
beta.prec <- 1/beta.sd^2

#--- Data file ---
wide.data <- list(n=100,chd=chd,age.scale=age.scale,beta.prec=beta.prec)

#--- Initial Values for Parameters
set.seed(1)
wide.inits <- list()
for(i in 1:num.chains) {
  # get "incompatable" initial values...?
  # b0 <- rnorm(n=1,mean=0,sd=beta.sd)
  # b1 <- rnorm(n=1,mean=0,sd=beta.sd)
  b0 <- rnorm(n=1,mean=0,sd=10)
  b1 <- rnorm(n=1,mean=0,sd=10)
  wide.inits[[i]] <- list(b0=b0,b1=b1)
}


wide.model <- "model {
 # priors
 b0 ~ dnorm(0,beta.prec)
 b1 ~ dnorm(0,beta.prec)

 # obs'n model
 for(i in 1:n) {
  logit(p[i]) <- b0+b1*age.scale[i]
  chd[i] ~ dbern(p[i])
 }
}
"
cat(wide.model,file="D:/Study/Postgraduate/paper/wide_model.txt")

#-- run jags
# Initialize the model
wide.init.JAGS <- jags.model( file = "D:/Study/Postgraduate/paper/wide_model.txt",
                              data=wide.data, inits=wide.inits, n.chains=num.chains,
                              n.adapt=adapt.length)

# specify parameters to monitor
wide.params <- c("b0","b1",'p')

# Run JAGS, using coda.samples a wrapper for jags.samples
wide.samples <- coda.samples(model=wide.init.JAGS,
                             variable.names=wide.params,
                             n.iter=burnin+chain.length)
logit.normal.matrix = rbind(wide.samples[[1]],wide.samples[[2]],wide.samples[[3]])

# plot(wide.samples)
wide.results <- summary(window(wide.samples,start=(burnin+1)))
print(wide.results$statistics)


# 2
#          Mean         SD     Naive SE Time-series SE
# b0     -0.39638007 0.23967968 0.0019569764   0.0027077664
# b1      1.34492515 0.28589835 0.0023343502   0.0038322077

# 2.5
#         Mean         SD     Naive SE Time-series SE
# b0     -0.39741821 0.24442517 0.0019957232   0.0028026210
# b1      1.35304930 0.28536679 0.0023300101   0.0037032616


# 25^2
#           Mean         SD     Naive SE Time-series SE
# b0     -0.39675681 0.24366670 0.0018161844   0.0025323714
# b1      1.34214038 0.29026355 0.0021634967   0.0034952621

c_alpha2 = c(wide.results$statistics[1,1]-1.96*wide.results$statistics[1,2]/sqrt(100), wide.results$statistics[1,1]+1.96*wide.results$statistics[1,2]/sqrt(100))
c_beta2 = c(wide.results$statistics[2,1]-1.96*wide.results$statistics[2,2]/sqrt(100), wide.results$statistics[2,1]+1.96*wide.results$statistics[2,2]/sqrt(100))
# > c_alpha2
# [1] -0.4445155 -0.3489981
# > c_beta2
# [1] 1.285249 1.399032


# 3.3.2 use jags and dcauchy(); use parameters for the Cauchy described in Seaman et al.

#------------------- Cauchy Prior Case ---------------------
beta.sd <- 2.5

#--- Data file ---
# Use age.std instead of age.std.5 for comparability
Cauchy.data <- list(n=100,
                    chd=chd,
                    #age.std.5=age.std.5,
                    age.scale=age.scale,
                    beta.sd=beta.sd)

#--- Initial Values for Parameters
set.seed(1)
Cauchy.inits <- list()
for(i in 1:num.chains){
  b0 <- rnorm(n=1,mean=0,sd=1)
  b1 <- rnorm(n=1,mean=0,sd=1)
  #--- occasionally getting initial values too "extreme"??
  #b0 <- rcauchy(n=1,location=0,scale=beta.scale)
  #b1 <- rcauchy(n=1,location=0,scale=beta.scale)
  Cauchy.inits[[i]] <- list(b0=b0,b1=b1)
}

Cauchy.model <- "model {
 # priors
 b0 ~ dt(0, pow(beta.sd,-2), 1)
 b1 ~ dt(0, pow(beta.sd,-2), 1)

 # obs'n model
 for(i in 1:n) {
  logit(p[i]) <- b0+b1*age.scale[i]
  chd[i] ~ dbern(p[i])
 }
}
"
cat(Cauchy.model,file="D:/Study/Postgraduate/paper/Cauchy_model.txt")

#-- run jags
# Initialize the model
Cauchy.init.JAGS <-  jags.model( file = "D:/Study/Postgraduate/paper/Cauchy_model.txt",
                                 data=Cauchy.data, inits=Cauchy.inits, n.chains=num.chains,
                                 n.adapt=adapt.length)

# specify parameters to monitor
Cauchy.params <- c("b0","b1",'p')

# Run JAGS, using coda.samples a wrapper for jags.samples
Cauchy.samples <- coda.samples(model=Cauchy.init.JAGS,
                               variable.names=Cauchy.params,
                               n.iter=burnin+chain.length)
logit.cauchy.matrix = rbind(Cauchy.samples[[1]],Cauchy.samples[[2]],Cauchy.samples[[3]])

# plot(Cauchy.samples)
Cauchy.results <- summary(window(Cauchy.samples,start=(burnin+1)))
print(Cauchy.results$statistics)
# Sample size per chain = 6000 
#          Mean         SD     Naive SE Time-series SE
# b0     -0.38620428 0.23880099 0.0017799175   0.0023270513
# b1      1.31665219 0.28539928 0.0021272406   0.0028968774
c_alpha3 = c(Cauchy.results$statistics[1,1]-1.96*Cauchy.results$statistics[1,2]/sqrt(100), Cauchy.results$statistics[1,1]+1.96*Cauchy.results$statistics[1,2]/sqrt(100))
c_beta3 = c(Cauchy.results$statistics[2,1]-1.96*Cauchy.results$statistics[2,2]/sqrt(100), Cauchy.results$statistics[2,1]+1.96*Cauchy.results$statistics[2,2]/sqrt(100))
c_alpha3
c_beta3
# > c_alpha3
# [1] -0.4330093 -0.3393993
# > c_beta3
# [1] 1.260714 1.372590

# 3.3.3 Use Christensen et al method on the CHD data.
# decide the value of b
#------------------- Christensen et al Prior Case ---------------------
b       <- 0.5   # much like pi prior
beta.sd <- 1/b # S^2 = 1
beta.prec <- 1/beta.sd

#--- Data file ---
Chris.data <- list(n=100,chd=chd,age.scale=age.scale,
                   beta.prec=beta.prec)

#--- Initial Values for Parameters
set.seed(1)
Chris.inits <- list()
for(i in 1:num.chains) {
  b0 <- rnorm(n=1,mean=0,sd=beta.sd)
  b1 <- rnorm(n=1,mean=0,sd=beta.sd)
  Chris.inits[[i]] <- list(b0=b0,b1=b1)
}


Chris.model <- "model {
 # priors
 b0 ~ dnorm(0,beta.prec)
 b1 ~ dnorm(0,beta.prec)

 # obs'n model
 for(i in 1:n) {
  logit(p[i]) <- b0+b1*age.scale[i]
  chd[i] ~ dbern(p[i])
 }
}
"
cat(Chris.model,file="D:/Study/Postgraduate/paper/Chris_model.txt")

#-- run jags
# Initialize the model
Chris.init.JAGS <- jags.model(file = "D:/Study/Postgraduate/paper/Chris_model.txt",
                              data=Chris.data, inits=Chris.inits, n.chains=num.chains,
                              n.adapt=adapt.length)

# specify parameters to monitor
Chris.params <- c("b0","b1",'p')

# Run JAGS, using coda.samples a wrapper for jags.samples
Chris.samples <- coda.samples(model=Chris.init.JAGS,
                              variable.names=Chris.params,
                              n.iter=burnin+chain.length)
#plot(Chris.samples)
Chris.results <- summary(window(Chris.samples,start=(burnin+1)))
print(Chris.results$statistics)
#            Mean         SD     Naive SE Time-series SE
# b0     -0.38419882 0.23943989 0.0019550185   0.0026583853
# b1      1.29612109 0.27568579 0.0022509651   0.0035095754c_alpha4 = c(Chris.results$statistics[1,1]-1.96*Chris.results$statistics[1,2]/sqrt(100), Chris.results$statistics[1,1]+1.96*Chris.results$statistics[1,2]/sqrt(100))
c_alpha4 = c(Chris.results$statistics[1,1]-1.96*Chris.results$statistics[1,2]/sqrt(100), Chris.results$statistics[1,1]+1.96*Chris.results$statistics[1,2]/sqrt(100))
c_beta4 = c(Chris.results$statistics[2,1]-1.96*Chris.results$statistics[2,2]/sqrt(100), Chris.results$statistics[2,1]+1.96*Chris.results$statistics[2,2]/sqrt(100))
c_alpha4
c_beta4
# > c_alpha4
# [1] -0.4311290 -0.3372686
# > c_beta4
# [1] 1.242087 1.350156

logit.chris.matrix = rbind(Chris.samples[[1]],Chris.samples[[2]],Chris.samples[[3]])

# 3.3.4 Use the normal(0,pi^2/(3*(p+1))) prior
#------------------- pi Prior Case ---------------------
p       <- 1
beta.sd <- pi/sqrt(3*(p+1))
beta.prec <- 1/beta.sd^2

#--- Data file ---
pi.data <- list(n=100,chd=chd,age.scale=age.scale,beta.prec=beta.prec)

#--- Initial Values for Parameters
set.seed(1)
pi.inits <- list()
for(i in 1:num.chains) {
  b0 <- rnorm(n=1,mean=0,sd=beta.sd)
  b1 <- rnorm(n=1,mean=0,sd=beta.sd)
  pi.inits[[i]] <- list(b0=b0,b1=b1)
}


pi.model <- "model {
 # priors
 b0 ~ dnorm(0,beta.prec)
 b1 ~ dnorm(0,beta.prec)

 # obs'n model
 for(i in 1:n) {
  logit(p[i]) <- b0+b1*age.scale[i]
  chd[i] ~ dbern(p[i])
 }
}
"
cat(pi.model,file="D:/Study/Postgraduate/paper/pi_model.txt")

#-- run jags
# Initialize the model
pi.init.JAGS <- jags.model( file = "D:/Study/Postgraduate/paper/pi_model.txt",
                            data=pi.data, inits=pi.inits, n.chains=num.chains,
                            n.adapt=adapt.length)

# specify parameters to monitor
pi.params <- c("b0","b1",'p')

# Run JAGS, using coda.samples a wrapper for jags.samples
pi.samples <- coda.samples(model=pi.init.JAGS,
                           variable.names=pi.params,
                           n.iter=burnin+chain.length)
#create the matrix
logit.pi.matrix = rbind(pi.samples[[1]],pi.samples[[2]],pi.samples[[3]])

#plot(pi.samples)
pi.results <- summary(window(pi.samples,start=(burnin+1)))
print(pi.results$statistics)
#             Mean         SD     Naive SE Time-series SE
# b0     -0.37260935 0.23307097 0.0019030165   0.0025556876
# b1      1.27957423 0.27146313 0.0022164872   0.0033466349
c_alpha5 = c(pi.results$statistics[1,1]-1.96*pi.results$statistics[1,2]/sqrt(100), pi.results$statistics[1,1]+1.96*pi.results$statistics[1,2]/sqrt(100))
c_beta5 = c(pi.results$statistics[2,1]-1.96*pi.results$statistics[2,2]/sqrt(100), pi.results$statistics[2,1]+1.96*pi.results$statistics[2,2]/sqrt(100))
c_alpha5
c_beta5
# > c_alpha5
# [1] -0.4182913 -0.3269274
# > c_beta5
# [1] 1.226367 1.332781

# 3.3.5 CMPs
#CMPs
x40 = (40-mean(aplore3::chdage$age))/sd(aplore3::chdage$age)
x60 = (60-mean(aplore3::chdage$age))/sd(aplore3::chdage$age)
modelString_cmp="
model{
  pi40 ~ dbeta(1,1)
  pi60 ~ dbeta(1,1)
  m1 = log(pi40/(1-pi40))
  m2 = log(pi60/(1-pi60))
  beta = (m2-m1)/(x60-x40)
  alpha = m1-beta*x40
  for (i in 1:100){
      y[i] ~ dbern(p[i])
      logit(p[i]) = alpha+beta*x[i]
  }
}"

writeLines(modelString_cmp, con='D:/Study/Postgraduate/paper/modelString_cmp.bug')
# create the input dataset
input.chdage.cmp = list('x' = age.scale,'y' = chd,'x40' = x40,'x60' = x60) 

chd.jags.cmp <- jags(data=input.chdage.cmp,model.file='D:/Study/Postgraduate/paper/modelString_cmp.bug',
                            param=c('alpha','beta','p'),
                            n.chains=3, n.iter=19000,n.burnin = 1000, n.thin=3)
chd.cmp.output = chd.jags.cmp$BUGSoutput$sims.matrix
chd.jags.cmp$BUGSoutput$summary
#               mean         sd         2.5%          25%          50%          75%        97.5%
# alpha     -0.39605972 0.23799752  -0.87264508  -0.55384034
# beta       1.27577720 0.27575341   0.76668393   1.08394484
# deviance 109.31022143 1.96202959 107.40460443 107.92941800chd.cmp.matrix = chd.cmp.output[,4:103]
c_alpha5 = c(-0.396-1.96*0.238/sqrt(100),-0.396+1.96*0.238/sqrt(100))
c_beta5 = c(1.276-1.96*0.276/sqrt(100), 1.276+1.96*0.276/sqrt(100))
# > c_alpha5
# [1] -0.442648 -0.349352
# > c_beta5
# [1] 1.221904 1.330096

# #------------------- CMP methods ---------------------
# 
# #--- Data file ---
# chdage.cmp.data = list('x' = age.scale,'y' = chd,'x40' = x40,'x60' = x60)
# 
# #--- Initial Values for Parameters
# set.seed(1)
# cmp.inits <- list()
# for(i in 1:num.chains) {
#   b1 = rnorm(n=1,1,1)
#   b0 = rnorm(n=1,1,1)
#   cmp.inits[[i]] <- list(b0=b0,b1=b1)
# }
# # every chain has an initial value
# 
# cmp.model <- "
# model{
#   for (i in 1:100){
#       y[i] ~ dbern(p[i])
#       logit(p[i]) = b0+b1*x[i]
#   }
#   pi40 ~ dbeta(1,1)
#   pi60 ~ dbeta(1,1)
#   m1 = log(pi40/(1-pi40))
#   m2 = log(pi60/(1-pi60))
#   b1 = (m2-m1)/(x60-x40)
#   b0 = (m1*x60-m2*x40)/(x60-x40)
# }"
# 
# cat(cmp.model,file="D:/Study/Postgraduate/paper/cmp_model.txt")
# 
# #-- run jags
# # Initialize the model
# cmp.init.JAGS <- jags.model(file = "D:/Study/Postgraduate/paper/cmp_model.txt",
#                             data=chdage.cmp.data, inits=cmp.inits, n.chains=num.chains,
#                             n.adapt=adapt.length)
# 
# # specify parameters to monitor
# cmp.params <- c("b0","b1","p")
# 
# # Run JAGS, using coda.samples a wrapper for jags.samples
# cmp.samples <- coda.samples(model=cmp.init.JAGS,
#                            variable.names=cmp.params,
#                            n.iter=burnin+chain.length)
# 
# plot(cmp.samples)
# cmp.results <- summary(window(cmp.samples,start=(burnin+1)))
# print(cmp.results$statistics)
# logit.cmps.matrix = rbind()
# 
# # output the results
# print(chd.jags.cmp)
# # n.sims = 3000 iterations saved
# #           mu.vect sd.vect    2.5%     25%     50%     75%   97.5%  Rhat n.eff
# # alpha     -0.394   0.234  -0.870  -0.546  -0.386  -0.236   0.048 1.002  2100
# # beta       1.276   0.277   0.761   1.080   1.263   1.451   1.848 1.001  3000
# 
# str(chd.jags.cmp$BUGSoutput$sims.matrix)
# 
# c_alpha5 = c(-0.394-1.96*0.234/sqrt(100),-0.394+1.96*0.234/sqrt(100))
# c_beta5 = c(1.276-1.96*0.277/sqrt(100), 1.276+1.96*0.277/sqrt(100))
# # > c_alpha5
# # [1] -0.439864 -0.348136
# # > c_beta5
# # [1] 1.221708 1.330292
# -0.222+1.323*x60

# > c_alpha5
# [1] -0.428768 -0.347232
# > c_beta5
# [1] 1.02694 1.11906

num.chains   <- 3
adapt.length <- 1000
burnin       <- 1000
chain.length <- 5000


# posterior distribution
# draw the prediction graph of the mle method
library(reshape2)

chd.model=glm(formula = chd~age.scale,family = binomial(link = 'logit'))
summary(chd.model)

fitted.prob.living = predict(chd.glm,type = "response")
plot(sort(age.scale),
     expit(-0.387+1.300*age.scale),type = 'l')
mtext('x60',side=1,line=x60)
age.order = (aplore3::chdage$age)


# N(0,sigma^2)
fitted.bayes1 = expit(-0.397+1.342*age.scale)
# Cauchy
fitted.bayes2 = expit(-0.386+1.317*age.scale)
# Chris
fitted.bayes3 = expit(-0.384+1.296*age.scale)
# N(0,pi^2/3)
fitted.bayes4 = expit(-0.373+1.280*age.scale)
# CMPs
fitted.bayes5 = expit(-0.396+1.276*age.scale)

# N(0,sigma^2)
lines(sort(age.scale),fitted.bayes1[order(age.scale)],lty = 2, col = 'red')
abline(h = expit(-0.397+1.342*x60),col = 'red')
# Cauchy
lines(sort(age.scale),fitted.bayes2[order(age.scale)],lty = 2, col = 'blue')
abline(h = expit(-0.386+1.317*x60),col = 'blue',lty = 2)
# Chris
lines(sort(age.scale),fitted.bayes3[order(age.scale)],lty = 2, col = 'orange')
abline(h = expit(-0.384+1.296*x60),col = 'orange',lty = 2)
# N(0,pi^2/3)
lines(sort(age.scale),fitted.bayes4[order(age.scale)],lty = 2, col = 'pink')
abline(h = expit(-0.373+1.280*x60),col = 'pink',lty = 2)
# CMPs
lines(sort(age.scale),fitted.bayes5[order(age.scale)],lty = 2, col = 'purple')
abline(h = expit(-0.396+1.276*x60),col = 'purple',lty = 2)
abline(v = x60,lty = 2)

axis(1,at=20:69,las = 3)
legend("topleft",legend = c("MLE","N(0,25^2)","Cauchy","N(0,S^2/b)(b=0.5)","N(0,pi^2/3)","CMPs"),
       col = c('black','red','blue','green','pink','purple'),lty = c(1,2,2,2,2,2),cex = 1)




par(mfrow=c(3,2))
plot(x = fitted.prob.living,y = fitted.bayes1,xlab = "CHD rate(MLE)",ylab = "CHD rate",pch = 20,main = expression(N(0,1.2^2)))
abline(a = 0 ,b = 1, col  = "blue")
plot(x = fitted.prob.living,y = fitted.bayes2,xlab = "CHD rate(MLE)",ylab = "CHD rate",pch = 20,main = "Cauchy(0,0.5)" )
abline(a = 0 ,b = 1, col  = "blue")
plot(x = fitted.prob.living,y = fitted.bayes3,xlab = "CHD rate(MLE)",ylab = "CHD rate",pch = 20,main = "N(0,S^2/b)" )
abline(a = 0 ,b = 1, col  = "blue")
plot(x = fitted.prob.living,y = fitted.bayes4,xlab = "CHD rate(MLE)",ylab = "CHD rate",pch = 20,main = "N(0,pi^2/3)" )
abline(a = 0 ,b = 1, col  = "blue")
plot(x = fitted.prob.living,y = fitted.bayes5,xlab = "CHD rate(MLE)",ylab = "CHD rate",pch = 20,main = "CMPs" )
abline(a = 0 ,b = 1, col  = "blue")

# show the fluctuation of b

data.mle.x = c()
data.normal.y = c()
data.cauchy.y = c()
data.chris.y = c()
data.pi.y = c()
data.cmps.y = c()

data.mle.x = rep(fitted.prob.living,18000)

posterior.survival.matrix.normal <- logit.normal.matrix[,3:102]
data.normal.y <- as.vector(t(posterior.survival.matrix.normal))

posterior.survival.matrix.cauchy <- logit.cauchy.matrix[,3:102]
data.cauchy.y <- as.vector(t(posterior.survival.matrix.cauchy))

posterior.survival.matrix.chris <- logit.chris.matrix[,3:102]
data.chris.y <- as.vector(t(posterior.survival.matrix.chris))

posterior.survival.matrix.pi <- logit.pi.matrix[,3:102]
data.pi.y <- as.vector(t(posterior.survival.matrix.pi))

data.cmps.y <- as.vector(t(chd.cmp.matrix))


# show the deviation of the Bayesian analysis

subset.values <- sort(sample(1:1800000,size=10000,replace=FALSE))
data.mle.x.subset = data.mle.x[subset.values]
data.normal.y.subset = data.normal.y[subset.values]
data.cauchy.y.subset = data.cauchy.y[subset.values]
data.chris.y.subset = data.chris.y[subset.values]
data.pi.y.subset = data.pi.y[subset.values]
data.cmps.y.subset = data.cmps.y[subset.values]

par(mfrow = c(2,3))
plot(data.mle.x.subset,data.normal.y.subset,pch = 20,xlab = "Estimate CHD rate (Frequencist)",ylab = "Estimate CHD rate (Bayesian)",main = "N(0,25^2)")
abline(a = 0 ,b = 1, col  = "blue")
plot(data.mle.x.subset,data.cauchy.y.subset,pch = 20,xlab = "Estimate CHD rate (Frequencist)",ylab = "Estimate CHD rate (Bayesian)",main = "Cauchy distribution")
abline(a = 0 ,b = 1, col  = "blue")
plot(data.mle.x.subset,data.chris.y.subset,pch = 20,xlab = "Estimate CHD rate (Frequencist)",ylab = "Estimate CHD rate (Bayesian)",main = "N(0,S^2/b),b = 0.5")
abline(a = 0 ,b = 1, col  = "blue")
plot(data.mle.x.subset,data.pi.y.subset,pch = 20,xlab = "Estimate CHD rate (Frequencist)",ylab = "Estimate CHD rate (Bayesian)",main = "N(0,pi^2/3)")
abline(a = 0 ,b = 1, col  = "blue")
plot(data.mle.x.subset,data.cmps.y.subset,pch = 20,xlab = "Estimate CHD rate (Frequencist)",ylab = "Estimate CHD rate (Bayesian)",main = "CMPs methods")
abline(a = 0 ,b = 1, col  = "blue")


# Chapter 4
# [Multivariate regression]
# dataset burn1000
library(aplore3)
str(aplore3::burn1000)

y = numeric(1000)
# Die or not
y[aplore3::burn1000$death=="Death"] = 1
# hopspital discharge status: 1 = Alive, 0 = Dead
# gender:1 = female, 0 = female
# race: 0 = white, 1 = non-white
# inh_inj: 1 = Yes, 0 = No
# flame: 1 = Yes, 0 = No
n      <- nrow(burn1000)
age    <- burn1000$age
facility <- burn1000$facility
tbsa   <- burn1000$tbsa
#gender <- burn1000$gender
#race   <- burn1000$race
#inh_inj <- burn1000$inh_inj
#flame <- burn1000$flame

# Make "death" binary 0,1
death <- numeric(n)
death[burn1000$death=="Dead"] <- 1

#--- standardize all  covariates
age.burn.scale = as.vector(scale(age))
tbsa.scale = as.vector(scale(tbsa))
facility.scale = as.vector(scale(facility))
gender =  numeric(1000)
gender[burn1000$gender=="Female"] = 1
# gender.scale = as.vector(scale(gender.scale))
race = numeric(1000)
race[burn1000$race=="Non-White"] = 1
# race.scale = as.vector(scale(race.scale))
inh_inj = numeric(1000)
inh_inj[burn1000$inh_inj=="Yes"] = 1
# inh_inj.scale = as.vector(scale(inh_inj.scale))
flame = numeric(1000)
flame[burn1000$flame=="Yes"] = 1
# flame.scale = as.vector(scale(flame.scale))

input.burn  <-  data.frame(death=death,age.scale=age.burn.scale,tbsa.scale = tbsa.scale,
                           facility.scale = facility.scale, gender = gender,
                           race = race, inh_inj = inh_inj,
                           flame = flame)

str(input.burn)
# 'data.frame':	1000 obs. of  8 variables:
# $ death         : num  0 0 0 0 0 0 0 0 0 0 ...
# $ age.scale     : num  -0.271 -1.27 -0.458 0.163 0.763 ...
# $ tbsa.scale    : num  0.616 -0.448 -0.605 -0.605 -0.396 ...
# $ facility.scale: num  -0.052 -0.9876 0.0415 -0.9876 -0.9876 ...
# $ gender        : num  0 1 1 0 0 0 1 1 0 0 ...
# $ race          : num  0 1 1 0 0 0 1 0 0 0 ...
# $ inh_inj       : num  0 0 0 0 0 0 0 0 0 0 ...
# $ flame         : num  1 0 0 0 1 0 0 1 0 1 ...
flame.glm = glm(death~facility.scale+age.scale+tbsa.scale+gender+race+inh_inj+flame,
                family = binomial,data = input.burn)
summary(flame.glm)
# Coefficients:
#                Estimate Std. Error z value Pr(>|z|)    
# (Intercept)     -4.6273     0.4457 -10.381  < 2e-16 ***
# facility.scale  -0.1762     0.1391  -1.266 0.205347    
# age.scale        2.0756     0.2174   9.546  < 2e-16 ***
# tbsa.scale       1.7411     0.1795   9.698  < 2e-16 ***
# gender           0.1531     0.3113   0.492 0.622991    
# race             0.7063     0.3109   2.272 0.023079 *  
# inh_inj          1.3409     0.3626   3.698 0.000217 ***
# flame            0.5829     0.3564   1.636 0.101927    
flame.glm = glm(death~input.burn$age.scale+input.burn$race+input.burn$tbsa.scale+input.burn$inh_inj,
                family = binomial,data = input.burn)
summary(flame.glm)
# Coefficients:
#                Estimate Std. Error z value Pr(>|z|)    
# (Intercept)    -4.1821     0.3464 -12.074  < 2e-16 ***
#   age.scale     2.0812     0.2091   9.954  < 2e-16 ***
#   race          0.6235     0.2989   2.086    0.037 *  
#   tbsa.scale    1.7247     0.1733   9.951  < 2e-16 ***
#   inh_inj       1.5231     0.3512   4.337 1.45e-05 ***
   
# 4.2 induced priors for christensen et al.
sd1 = 1
sd2 = 1
sd3 = 1
sd4 = 1

# Find the suitable b value in christensen et al. method
induced_priors_burn_chris = function(b){
  induced.prior.prob.burn = matrix(data = NA,nrow = 1000, ncol = 1000)
  
  set.seed(1)
  
  alpha = rnorm(1000,0,(sd1^2+sd^2+sd3^2+sd4^2)/b)
  beta1 = rnorm(1000,0,sqrt(sd1^2/b))
  beta2 = rnorm(1000,0,sqrt(sd2^2/b))
  beta3 = rnorm(1000,0,sqrt(sd3^2/b))
  beta4 = rnorm(1000,0,sqrt(sd4^2/b))
  for(i in 1:1000){
    induced.prior.prob.burn[i,] = expit(alpha+beta1*input.burn$age.scale[i]+beta2*input.burn$tbsa.scale[i]+beta3*input.burn$race.scale[i]+beta4*input.burn$inh_inj.scale[i])
  }
  hist(induced.prior.prob.burn,xlab = expression(pi(x)),main = b)
}
# b = 1000
par(mfrow = c(2,3))
induced_priors_burn_chris(4.5)
induced_priors_burn_chris(4)
induced_priors_burn_chris(3.5)
induced_priors_burn_chris(3)
induced_priors_burn_chris(2.5)
induced_priors_burn_chris(2)


# b = 3

# Other priors

# gai!
# 4.3 Induced priors
par(mfrow=c(2,2))

# 4.3.1 N(0,10^6)
par(mfrow = c(5,2))
set.seed(1)
logit.mul.normal = function(sigma){
induced.prior.normal = matrix(data = NA,nrow = 1000, ncol = 1000)
induced.prior.random = matrix(data = NA,nrow = 1000, ncol = 1000)
alpha = rnorm(1000,0,sigma)
beta1 = rnorm(1000,0,sigma)
beta2 = rnorm(1000,0,sigma)
beta3 = rnorm(1000,0,sigma)
beta4 = rnorm(1000,0,sigma)
x.age = rnorm(1000,0,1)
x.tbsa = rnorm(1000,0,1)
x.race = rnorm(1000,0,1)
x.inh_inj = rnorm(1000,0,1)
for(i in 1:1000){
  induced.prior.normal[i,] = expit(alpha+beta1*input.burn$age.scale[i]+beta2*input.burn$tbsa.scale[i]+
                                     beta3*input.burn$race.scale[i]+beta4*input.burn$inh_inj.scale[i]) 
  induced.prior.random[i,] = expit(alpha+beta1*x.age[i]+beta2*x.tbsa[i]+beta3*x.race[i]+beta4*x.inh_inj[i]) 
}
hist(induced.prior.normal,xlab = expression(pi(x)), ylab=("frequency"),main = paste("N(0,",sigma,"^2 )original dataset"))
hist(induced.prior.random,xlab = expression(pi(x)), ylab=("frequency"),main = paste("N(0,",sigma,"^2 )random variable"))
}


par(mfrow = c(4,2))
logit.mul.normal(10)

logit.mul.normal(2.5)
logit.mul.normal(1.75)
logit.mul.normal(1.5)
logit.mul.normal(1.25)
logit.mul.normal(1)
logit.mul.normal(0.75)
# b = 1


# hist(1/(1+1/exp(E)),xlab =expression(pi(60)),main = "histogram of induced prior (N(0,10^6))")

# 4.3.2 Cauchy priors
cauchy.mul.induced=function(scale){
sd1 = sd(input.burn$age.scale)
sd2 = sd(input.burn$tbsa.scale)
sd3 = sd(input.burn$race)
sd4 = sd(input.burn$inh_inj)
set.seed(1)
induced.prior.cauchy = matrix(data = NA,nrow = 1000, ncol = 1000)
induced.prior.random = matrix(data = NA,nrow = 1000, ncol = 1000)
alpha = rcauchy(1000,0,4*scale)
beta1 = rcauchy(1000,0,scale)
beta2 = rcauchy(1000,0,scale)
beta3 = rcauchy(1000,0,scale/(2*sd3))
beta4 = rcauchy(1000,0,scale/(2*sd4))
x.age = rnorm(1000,0,1)
x.tbsa = rnorm(1000,0,1)
x.race = rnorm(1000,0,1)
x.inh_inj = rnorm(1000,0,1)
for(i in 1:1000){
  induced.prior.cauchy[i,] = expit(alpha+beta1*input.burn$age.scale[i]+beta2*input.burn$tbsa.scale[i]+beta3*input.burn$race[i]+beta4*input.burn$inh_inj[i])
  induced.prior.random[i,] = expit(alpha+beta1*x.age[i]+beta2*x.tbsa[i]+beta3*x.race[i]+beta4*x.inh_inj[i]) 
  }
hist(induced.prior.cauchy,xlab = expression(pi(x)), ylab=("frequency"),main = paste('Cauchy distribution original dataset(scale =',scale,')'))
hist(induced.prior.random,xlab = expression(pi(x)), ylab=("frequency"),main = "Cauchy distribution random variable")
}

par(mfrow = c(5,2))
cauchy.mul.induced(0.1)
cauchy.mul.induced(0.25)
cauchy.mul.induced(0.5)
cauchy.mul.induced(1)
cauchy.mul.induced(2.5)




# 4.3.3 N(0,S^2/b)
# induced priors
# b = 0.1
# alpha = -5.505 , beta = 0.115
#1/(S_2/b)
#[1] 0.0007278577

# sd1 = sd2 = sd3 = sd4 = 1
b = 0.1
logit.chris.mul.induced = function(b){
induced.prior.chris = matrix(data = NA,nrow = 1000, ncol = 1000)
induced.prior.random = matrix(data = NA,nrow = 1000, ncol = 1000)
set.seed(1)
alpha = rnorm(1000,0,(1^2+1^2+1^2+1^2)/b)
beta1 = rnorm(1000,0,sqrt(1^2/b))
beta2 = rnorm(1000,0,sqrt(1^2/b))
beta3 = rnorm(1000,0,sqrt(1^2/b))
beta4 = rnorm(1000,0,sqrt(1^2/b))
x.age = rnorm(1000,0,1)
x.tbsa = rnorm(1000,0,1)
x.race = rnorm(1000,0,1)
x.inh_inj = rnorm(1000,0,1)

for(i in 1:1000){
  induced.prior.chris [i,] = expit(alpha+beta1*input.burn$age.scale[i]+beta2*input.burn$tbsa.scale[i]+
                                     beta3*input.burn$race[i]+beta4*input.burn$inh_inj[i])
  induced.prior.random [i,] = expit(alpha+beta1*x.age[i]+beta2*x.tbsa[i]+
                                     beta3*x.race[i]+beta4*x.inh_inj[i])
}
hist(induced.prior.chris,xlab = expression(pi(x)),main = paste("Christensen method Original dataset, b =",b))
hist(induced.prior.random,xlab = expression(pi(x)),main = paste("Christensen method Random variable, b =",b))
}

par(mfrow = c(5,2))
logit.chris.mul.induced(4)
logit.chris.mul.induced(3.5)
logit.chris.mul.induced(3)
logit.chris.mul.induced(2.5)
logit.chris.mul.induced(2)







# 4.3.4 N(0,pi^2/(3*(p+1)))
par(mfrow = c(1,2))
p = 4
set.seed(1)
alpha = rnorm(1000,0,pi^2/(3*(p+1)))
beta1 = rnorm(1000,0,pi^2/(3*(p+1)))
beta2 = rnorm(1000,0,pi^2/(3*(p+1)))
beta3 = rnorm(1000,0,pi^2/(3*(p+1)))
beta4 = rnorm(1000,0,pi^2/(3*(p+1)))
x.age = rnorm(1000,0,1)
x.tbsa = rnorm(1000,0,1)
x.race = rnorm(1000,0,1)
induced.prior.pi = matrix(data = NA,nrow = 1000, ncol = 1000)
induced.prior.random = matrix(data = NA,nrow = 1000, ncol = 1000)
induced.prior.pi
for(i in 1:1000){
  induced.prior.pi [i,] = expit(alpha+beta1*df.flame$age.scale[i]+beta2*df.flame$tbsa.scale[i]+beta3*df.flame$race[i]+beta4*df.flame$inh_inj[i])
  induced.prior.random [i,] = expit(alpha+beta1*x.age[i]+beta2*x.tbsa[i]+
                                      beta3*x.race[i]+beta4*x.inh_inj[i])
  }
hist(induced.prior.pi ,xlab = expression(pi(x)),main = "N(0,pi^2/(3*(p+1)) original dataset")
hist(induced.prior.random,xlab = expression(pi(x)),main = "N(0,pi^2/(3*(p+1)) random variable")


par(mfrow = c(2,2))


# N(0,sigma^2)
sigma = 10^3
induced.prior.normal = matrix(data = NA,nrow = 1000, ncol = 1000)
alpha = rnorm(1000,0,sigma)
beta1 = rnorm(1000,0,sigma)
beta2 = rnorm(1000,0,sigma)
beta3 = rnorm(1000,0,sigma)
beta4 = rnorm(1000,0,sigma)
x.age = rnorm(1000,0,1)
x.tbsa = rnorm(1000,0,1)
x.race = rnorm(1000,0,1)
x.inh_inj = rnorm(1000,0,1)
for(i in 1:1000){
  induced.prior.normal[i,] = expit(alpha+beta1*input.burn$age.scale[i]+beta2*input.burn$tbsa.scale[i]+
                                     beta3*input.burn$race[i]+beta4*input.burn$inh_inj[i]) 
}
hist(induced.prior.normal,xlab = expression(pi(x)), ylab=("frequency"),main = paste("N(0,",sigma,"^2 )original dataset"))

# Cauchy distribution
sd1 = sd(input.burn$age.scale)
sd2 = sd(input.burn$tbsa.scale)
sd3 = sd(input.burn$race)
sd4 = sd(input.burn$inh_inj)
set.seed(1)
scale = 0.25
induced.prior.cauchy = matrix(data = NA,nrow = 1000, ncol = 1000)
alpha = rcauchy(1000,0,4*scale)
beta1 = rcauchy(1000,0,scale)
beta2 = rcauchy(1000,0,scale)
beta3 = rcauchy(1000,0,scale/(2*sd3))
beta4 = rcauchy(1000,0,scale/(2*sd4))
x.age = rnorm(1000,0,1)
x.tbsa = rnorm(1000,0,1)
x.race = rnorm(1000,0,1)
x.inh_inj = rnorm(1000,0,1)
for(i in 1:1000){
  induced.prior.cauchy[i,] = expit(alpha+beta1*input.burn$age.scale[i]+beta2*input.burn$tbsa.scale[i]+beta3*input.burn$race[i]+beta4*input.burn$inh_inj[i])
}
hist(induced.prior.cauchy,xlab = expression(pi(x)), ylab=("frequency"),main = paste('Cauchy distribution original dataset(scale =',scale,')'))

#N(0,S^2/b)
b = 1
induced.prior.chris = matrix(data = NA,nrow = 1000, ncol = 1000)
set.seed(1)
alpha = rnorm(1000,0,(1^2+1^2+1^2+1^2)/b)
beta1 = rnorm(1000,0,sqrt(1^2/b))
beta2 = rnorm(1000,0,sqrt(1^2/b))
beta3 = rnorm(1000,0,sqrt(1^2/b))
beta4 = rnorm(1000,0,sqrt(1^2/b))
x.age = rnorm(1000,0,1)
x.tbsa = rnorm(1000,0,1)
x.race = rnorm(1000,0,1)
x.inh_inj = rnorm(1000,0,1)

for(i in 1:1000){
  induced.prior.chris [i,] = expit(alpha+beta1*input.burn$age.scale[i]+beta2*input.burn$tbsa.scale[i]+
                                     beta3*input.burn$race[i]+beta4*input.burn$inh_inj[i])
}
hist(induced.prior.chris,xlab = expression(pi(x)),main = paste("Christensen method Original dataset, b =",b))

# N(0,pi^2/3*(p+1))
p = 4
set.seed(1)
alpha = rnorm(1000,0,pi^2/(3*(p+1)))
beta1 = rnorm(1000,0,pi^2/(3*(p+1)))
beta2 = rnorm(1000,0,pi^2/(3*(p+1)))
beta3 = rnorm(1000,0,pi^2/(3*(p+1)))
beta4 = rnorm(1000,0,pi^2/(3*(p+1)))
x.age = rnorm(1000,0,1)
x.tbsa = rnorm(1000,0,1)
x.race = rnorm(1000,0,1)
induced.prior.pi = matrix(data = NA,nrow = 1000, ncol = 1000)
for(i in 1:1000){
  induced.prior.pi [i,] = expit(alpha+beta1*df.flame$age.scale[i]+beta2*df.flame$tbsa.scale[i]+beta3*df.flame$race[i]+beta4*df.flame$inh_inj[i])
}
hist(induced.prior.pi ,xlab = expression(pi(x)),main = "N(0,pi^2/(3*(p+1)) original dataset")


# 4.3 Regression coefficients
# 4.3.1 MLE
# The model is death ~ age+tbsa+race+inh_inj
burn.model1=glm(formula = death~facility.scale+age.scale+gender+race+tbsa.scale+inh_inj+flame,family = binomial(link = 'logit'),data = input.burn)
summary(burn.model1) 
# Coefficients:
#                 Estimate Std. Error z value Pr(>|z|)    
#   (Intercept)     -4.6273     0.4457 -10.381  < 2e-16 ***
#   facility.scale  -0.1762     0.1391  -1.266 0.205347    
#   age.scale        2.0756     0.2174   9.546  < 2e-16 ***
#   gender           0.1531     0.3113   0.492 0.622991    
#   race             0.7063     0.3109   2.272 0.023079 *  
#   tbsa.scale       1.7411     0.1795   9.698  < 2e-16 ***
#   inh_inj          1.3409     0.3626   3.698 0.000217 ***
#   flame            0.5829     0.3564   1.636 0.101927  


burn.model2=glm(formula = death~age.scale+tbsa.scale+race+inh_inj,family = binomial(link = 'logit'),data = input.burn)

summary(burn.model2)
# Coefficients:
#               Estimate Std. Error z value Pr(>|z|)    
# (Intercept)    -4.1821     0.3464 -12.074  < 2e-16 ***
#   age.scale     2.0812     0.2091   9.954  < 2e-16 ***
#   tbsa.scale    1.7247     0.1733   9.951  < 2e-16 ***
#   race          0.6235     0.2989   2.086    0.037 *  
#   inh_inj       1.5231     0.3512   4.337 1.45e-05 ***


# 4.3.2 N(0,10e-6)

# 10e-6 1

# build the model
# According to the David W.HOSMER Applied Logistic Regression, p116
# The model is: death ~ age

# ***************** Bayesian Inference *******************
num.chains   <- 3
adapt.length <- 1000
burnin       <- 1000
chain.length <- 5000

#num.iteration = 1000+5000

#-- Wide Normal Prior Case ---------------------
beta.sd <- 1000
beta.prec <- 1/beta.sd^2

#--- Data file ---
wide.mul.data <- list(n=1000,death=input.burn$death,age.scale=input.burn$age.scale,tbsa.scale=input.burn$tbsa.scale,
                  race= input.burn$race,inh_inj=input.burn$inh_inj,
                  beta.prec=beta.prec)

#--- Initial Values for Parameters
set.seed(1)
wide.mul.inits <- list()
for(i in 1:num.chains) {
  # get "incompatable" initial values...?
  # b0 <- rnorm(n=1,mean=0,sd=beta.sd)
  # b1 <- rnorm(n=1,mean=0,sd=beta.sd)
  b0       = rnorm(n=1,0,1)
  b.age    = rnorm(n=1,0,1)
  b.tbsa   = rnorm(n=1,0,1)
  b.race   = rnorm(n=1,0,1)
  b.inh_inj = rnorm(n=1,0,1)
  wide.mul.inits[[i]] <- list(b0=b0,b.age=b.age,b.tbsa=b.tbsa,
                          b.race= b.race,b.inh_inj=b.inh_inj)
}

wide.mul.model <- "model {
 # priors
  b0       ~ dnorm(0,beta.prec)
  b.age    ~ dnorm(0,beta.prec)
  b.tbsa   ~ dnorm(0,beta.prec)
  b.race   ~ dnorm(0,beta.prec)
  b.inh_inj ~ dnorm(0,beta.prec)
  
 # obs'n model
 for(i in 1:n) {
  logit(p[i]) <- b0 + b.age*age.scale[i] + b.tbsa*tbsa.scale[i] +
      b.race*race[i] + b.inh_inj*inh_inj[i]
  death[i] ~ dbern(p[i])
 }
}
"
cat(wide.mul.model,file="D:/Study/Postgraduate/paper/wide_mul_model.txt")

#-- run jags
# Initialize the model
wide.mul.init.JAGS <- jags.model( file = "D:/Study/Postgraduate/paper/wide_mul_model.txt",
                              data=wide.mul.data, inits=wide.mul.inits, n.chains=num.chains,
                              n.adapt=adapt.length)

# specify parameters to monitor
wide.mul.params <- c("b0","b.age",'b.tbsa',"b.race",'b.inh_inj','p')

# Run JAGS, using coda.samples a wrapper for jags.samples
wide.mul.samples <- coda.samples(model=wide.mul.init.JAGS,
                             variable.names=wide.mul.params,
                             n.iter=burnin+chain.length)
logit.mul.normal.matrix = rbind(wide.mul.samples[[1]],wide.mul.samples[[2]],wide.mul.samples[[3]])



#plot(wide.samples)
wide.mul.results <- summary(window(wide.mul.samples,start=(burnin+1)))
print(wide.mul.rsults$statistics)
#                Mean           SD     Naive SE Time-series SE
# b.age      2.1271131589 0.2139639860 1.747009e-03   1.050825e-02
# b.inh_inj  1.5409683527 0.3604160544 2.942785e-03   5.800296e-03
# b.race     0.6386155012 0.3006732073 2.454986e-03   6.356027e-03
# b.tbsa     1.7650514179 0.1747424659 1.426766e-03   5.988582e-03
# b0        -4.2650133957 0.3581931086 2.924634e-03   1.667765e-02

# 4.3.3 Cauchy priors

#--- Data file ---
cauchy.mul.data <- list(n=1000,death=input.burn$death,age.scale=input.burn$age.scale,tbsa.scale=input.burn$tbsa.scale,
                      race= input.burn$race,inh_inj=input.burn$inh_inj)

#--- Initial Values for Parameters
set.seed(1)
cauchy.mul.inits <- list()
for(i in 1:num.chains) {
  # get "incompatable" initial values...?
  # b0 <- rnorm(n=1,mean=0,sd=beta.sd)
  # b1 <- rnorm(n=1,mean=0,sd=beta.sd)
  b0       = rnorm(n=1,0,1)
  b.age    = rnorm(n=1,0,1)
  b.tbsa   = rnorm(n=1,0,1)
  b.race   = rnorm(n=1,0,1)
  b.inh_inj = rnorm(n=1,0,1)
  cauchy.mul.inits[[i]] <- list(b0=b0,b.age=b.age,b.tbsa=b.tbsa,
                              b.race= b.race,b.inh_inj=b.inh_inj)
}

cauchy.mul.model <- "model {
 # priors
  b0       ~ dt(0,0.1,1)
  b.age    ~ dt(0,0.25/(2*1),1)
  b.tbsa   ~ dt(0,0.25/(2*1),1)
  b.race   ~ dt(0,0.25,1)
  b.inh_inj ~ dt(0,0.25,1)
  
 # obs'n model
 for(i in 1:n) {
  logit(p[i]) <- b0 + b.age*age.scale[i] + b.tbsa*tbsa.scale[i] +
      b.race*race[i] + b.inh_inj*inh_inj[i]
  death[i] ~ dbern(p[i])
 }
}
"
cat(cauchy.mul.model,file="D:/Study/Postgraduate/paper/wide_mul_model.txt")

#-- run jags
# Initialize the model
cauchy.mul.init.JAGS <- jags.model(file = "D:/Study/Postgraduate/paper/wide_mul_model.txt",
                                  data=cauchy.mul.data, inits=cauchy.mul.inits, n.chains=num.chains,
                                  n.adapt=adapt.length)

# specify parameters to monitor
cauchy.mul.params <- c("b0","b.age",'b.tbsa',"b.race",'b.inh_inj','p')

# Run JAGS, using coda.samples a wrapper for jags.samples
cauchy.mul.samples <- coda.samples(model=cauchy.mul.init.JAGS,
                                 variable.names=cauchy.mul.params,
                                 n.iter=burnin+chain.length)
logit.mul.cauchy.matrix = rbind(cauchy.mul.samples[[1]],cauchy.mul.samples[[2]],cauchy.mul.samples[[3]])

#plot(cauchy.mul.samples)
cauchy.mul.results <- summary(window(cauchy.mul.samples,start=(burnin+1)))
print(cauchy.mul.results$statistics)
#                Mean           SD       Naive SE Time-series SE
# b.age      2.0794544789 0.2067852258 0.001541286072 0.005459798560
# b.inh_inj  1.4739107169 0.3487317349 0.002599292884 0.004772989345
# b.race     0.5907813715 0.2950602283 0.002199249093 0.004949659365
# b.tbsa     1.7459820755 0.1737440249 0.001295011501 0.003199165025
# b0        -4.1649477935 0.3399591567 0.002533905947 0.009962411631


# 4.3.4 N(0, S2/b) Christensen et al.
# calculate the regression result on different values of b
# "Naive" b value
b.sq <- 3
Chris.sd <- 1/sqrt(b.sq)
beta.prec <- 1/Chris.sd^2
# b.sq <- 3*(p+1)/pi^2  # note this yields sd = 0.811
# s2_1 = var(input.burn$age.scale)
# s2_2 = var(input.burn$tbsa.scale)
# s2_3 = 1# dichotomous corvariates
# s2_4 = 1# dichotomous corvariates
# s2_0 ?

#--- Data file ---
# what if change the beta.prec

chris.mul.data <- list("death" = input.burn$death,
                       "age.scale" = input.burn$age.scale,
                       "tbsa.scale" = input.burn$tbsa.scale,
                       "race" = input.burn$race,
                       "inh_inj" = input.burn$inh_inj,
                       "beta.prec" = beta.prec,
                       "n"=1000)

#--- Initial Values for Parameters
set.seed(1)
chris.mul.inits <- list()
for(i in 1:num.chains) {
  # get "incompatable" initial values...?
  # b0 <- rnorm(n=1,mean=0,sd=beta.sd)
  # b1 <- rnorm(n=1,mean=0,sd=beta.sd)
  b0       = rnorm(n=1,0,Chris.sd)
  b.age    = rnorm(n=1,0,Chris.sd)
  b.tbsa   = rnorm(n=1,0,Chris.sd)
  b.race   = rnorm(n=1,0,Chris.sd)
  b.inh_inj = rnorm(n=1,0,Chris.sd)
  chris.mul.inits[[i]] <- list(b0=b0,b.age=b.age,b.tbsa=b.tbsa,
                               b.race= b.race,b.inh_inj=b.inh_inj)
}

chris.mul.model = wide.mul.model
cat(chris.mul.model,file="D:/Study/Postgraduate/paper/chris_mul_model.txt")

#-- run jags
# Initialize the model
chris.mul.init.JAGS <- jags.model(file = "D:/Study/Postgraduate/paper/chris_mul_model.txt",
                                   data=chris.mul.data, inits=chris.mul.inits,
                                   n.chains=num.chains,
                                   n.adapt=adapt.length)

# specify parameters to monitor
chris.mul.params <- c("b0","b.age",'b.tbsa',"b.race",'b.inh_inj','p')

# Run JAGS, using coda.samples a wrapper for jags.samples
t1 <- proc.time()[3]
chris.mul.samples <- coda.samples(model=chris.mul.init.JAGS,
                                  variable.names=chris.mul.params,
                                  n.iter=burnin+chain.length)
t2 <- proc.time()[3]
logit.mul.chris.matrix = rbind(chris.mul.samples[[1]],chris.mul.samples[[2]],chris.mul.samples[[3]])



cat("elapsed time=",(t2-t1),"seconds \n")
# elapsed time= 37.1 seconds

#plot(burn.Chris.samples,ask=TRUE)
chris.mul.results <- summary(window(chris.mul.samples,start=(burnin+1)))
print(chris.mul.results$statistics)
#               Mean           SD       Naive SE Time-series SE
# b.age      1.517618614 0.1397131232 0.001140752874  0.00374174513
# b.inh_inj  0.962484842 0.2777029988 0.002267435490  0.00327875550
# b.race     0.139749509 0.2340240507 0.001910798372  0.00383362555
# b.tbsa     1.467349198 0.1370857924 0.001119300808  0.00267174251
# b0        -3.104548716 0.2067447287 0.001688063641  0.00545471183


# 4.3.5 N(0,pi2/(3*(p+1)))
# sd is 3*(4+1)/pi^2 in jags
p       <- 4
beta.sd <- pi/sqrt(3*(p+1))
beta.prec <- 1/beta.sd^2

#--- Data file ---
pi.mul.data <- list("death" = input.burn$death,
                       "age.scale" = input.burn$age.scale,
                       "tbsa.scale" = input.burn$tbsa.scale,
                       "race" = input.burn$race,
                       "inh_inj" = input.burn$inh_inj,
                       "beta.prec" = beta.prec,
                       n=1000)

#--- Initial Values for Parameters
set.seed(1)
pi.mul.inits <- list()
for(i in 1:num.chains) {
  # get "incompatable" initial values...?
  # b0 <- rnorm(n=1,mean=0,sd=beta.sd)
  # b1 <- rnorm(n=1,mean=0,sd=beta.sd)
  b0       = rnorm(n=1,0,1)
  b.age    = rnorm(n=1,0,1)
  b.tbsa   = rnorm(n=1,0,1)
  b.race   = rnorm(n=1,0,1)
  b.inh_inj = rnorm(n=1,0,1)
  pi.mul.inits[[i]] <- list(b0=b0,b.age=b.age,b.tbsa=b.tbsa,
                               b.race= b.race,b.inh_inj=b.inh_inj)
}
#-- run jags
# Initialize the model

pi.mul.model = wide.mul.model
cat(pi.mul.model,file="D:/Study/Postgraduate/paper/pi_mul_model.txt")

#-- run jags
# Initialize the model
pi.mul.init.JAGS <- jags.model(file = "D:/Study/Postgraduate/paper/pi_mul_model.txt",
                                   data=pi.mul.data, inits=pi.mul.inits, n.chains=num.chains,
                                   n.adapt=adapt.length)

# specify parameters to monitor
pi.mul.params <- c("b0","b.age",'b.tbsa',"b.race",'b.inh_inj','p')

# Run JAGS, using coda.samples a wrapper for jags.samples
pi.mul.samples <- coda.samples(model=pi.mul.init.JAGS,
                                  variable.names=pi.mul.params,
                                  n.iter=burnin+chain.length)
logit.mul.pi.matrix = rbind(pi.mul.samples[[1]],pi.mul.samples[[2]],pi.mul.samples[[3]])



pi.mul.results <- summary(window(pi.mul.samples,start=(burnin+1)))
print(pi.mul.results$statistics)
#                Mean           SD     Naive SE Time-series SE
# b.age      1.714869322 0.1633050645 1.333380e-03   5.824164e-03
# b.inh_inj  1.155125439 0.3002378353 2.451432e-03   3.759670e-03
# b.race     0.293750126 0.2604734947 2.126757e-03   4.714010e-03
# b.tbsa     1.566980442 0.1482912213 1.210793e-03   3.660925e-03
# b0        -3.477508102 0.2474308708 2.020265e-03   8.446120e-03

c_alpha5 = c(pi.mul.results$statistics[1,1]-1.96*pi.mul.results$statistics[1,2]/sqrt(1000), pi.mul.results$statistics[1,1]+1.96*pi.mul.results$statistics[1,2]/sqrt(1000))
c_beta5 = c(pi.mul.results$statistics[2,1]-1.96*pi.mul.results$statistics[2,2]/sqrt(1000), pi.mul.results$statistics[2,1]+1.96*pi.mul.results$statistics[2,2]/sqrt(1000))
c_alpha5 
c_beta5
# > c_alpha5 
# [1] 1.776251 1.797592
# > c_beta5
# [1] 0.9907973 1.0194405

#4.5 posterior distribution

# b.age
par(mfrow = c(2,2))
x = sample(1:18000,2500,replace = FALSE)
plot(x = x,y=(logit.mul.normal.matrix [,1]/(2.081))[x],pch = 20,xlab = "simulation point",ylab = "b.age",main = "Normal")
abline(h = 1,col = "blue",lwd = 3)
plot(x = x,y=(logit.mul.cauchy.matrix [,1]/(2.081))[x],pch = 20,xlab = "simulation point",ylab = "b.age",main = "Cauchy")
abline(h = 1,col = "blue",lwd = 3)
plot(x = x,y=(logit.mul.chris.matrix [,1]/(2.081))[x],pch = 20,xlab = "simulation point",ylab = "b.age",main = "Christensen")
abline(h = 1,col = "blue",lwd = 3)
plot(x = x,y=(logit.mul.pi.matrix [,1]/(2.081))[x],pch = 20,xlab = "simulation point",ylab = "b.age",main = "N(0,pi2/(3*(p+1)))")
abline(h = 1,col = "blue",lwd = 3)

# b.inh_inj
par(mfrow = c(2,2))
x = c(1:18000)
plot(x = x,y=(logit.mul.normal.matrix [,2]/(1.110))[x],pch = 20,xlab = "simulation point",ylab = "beta0",main = "Normal")
abline(h = 1,col = "blue",lwd = 3)
plot(x = x,y=(logit.mul.cauchy.matrix [,2]/(1.110))[x],pch = 20,xlab = "simulation point",ylab = "beta0",main = "Cauchy")
abline(h = 1,col = "blue",lwd = 3)
plot(x = x,y=(logit.mul.chris.matrix [,2]/(1.110))[x],pch = 20,xlab = "simulation point",ylab = "beta0",main = "Christensen")
abline(h = 1,col = "blue",lwd = 3)
plot(x = x,y=(logit.mul.pi.matrix [,2]/(1.110))[x],pch = 20,xlab = "simulation point",ylab = "beta0",main = "N(0,pi2/(3*(p+1)))")
abline(h = 1,col = "blue",lwd = 3)

# b.race
par(mfrow = c(2,2))
x = c(1:18000)
plot(x = x,y=(logit.mul.normal.matrix [,3]/(0.340))[x],pch = 20,xlab = "simulation point",ylab = "beta0",main = "Normal")
abline(h = 1,col = "blue",lwd = 3)
plot(x = x,y=(logit.mul.cauchy.matrix [,3]/(0.340))[x],pch = 20,xlab = "simulation point",ylab = "beta0",main = "Cauchy")
abline(h = 1,col = "blue",lwd = 3)
plot(x = x,y=(logit.mul.chris.matrix [,3]/(0.340))[x],pch = 20,xlab = "simulation point",ylab = "beta0",main = "Christensen")
abline(h = 1,col = "blue",lwd = 3)
plot(x = x,y=(logit.mul.pi.matrix [,3]/(0.340))[x],pch = 20,xlab = "simulation point",ylab = "beta0",main = "N(0,pi2/(3*(p+1)))")
abline(h = 1,col = "blue",lwd = 3)

# b.tbsa
par(mfrow = c(2,2))
x = c(1:18000)
plot(x = x,y=(logit.mul.normal.matrix [,4]/1.725)[x],pch = 20,xlab = "simulation point",ylab = "beta0",main = "Normal")
abline(h = 1,col = "blue",lwd = 3)
plot(x = x,y=(logit.mul.cauchy.matrix [,4]/1.725)[x],pch = 20,xlab = "simulation point",ylab = "beta0",main = "Cauchy")
abline(h = 1,col = "blue",lwd = 3)
plot(x = x,y=(logit.mul.chris.matrix [,4]/1.725)[x],pch = 20,xlab = "simulation point",ylab = "beta0",main = "Christensen")
abline(h = 1,col = "blue",lwd = 3)
plot(x = x,y=(logit.mul.pi.matrix [,4]/1.725)[x],pch = 20,xlab = "simulation point",ylab = "beta0",main = "N(0,pi2/(3*(p+1)))")
abline(h = 1,col = "blue",lwd = 3)

# b0
par(mfrow = c(2,2))
x = c(1:6000)
plot(x = x,y=(logit.mul.normal.matrix [,5]/(-3.485))[x],pch = 20,xlab = "simulation point",ylab = "beta0",main = "Normal")
abline(h = 1,col = "blue",lwd = 3)
plot(x = x,y=(logit.mul.cauchy.matrix [,5]/(-3.485))[x],pch = 20,xlab = "simulation point",ylab = "beta0",main = "Cauchy")
abline(h = 1,col = "blue",lwd = 3)
plot(x = x,y=(logit.mul.chris.matrix [,5]/(-3.485))[x],pch = 20,xlab = "simulation point",ylab = "beta0",main = "Christensen")
abline(h = 1,col = "blue",lwd = 3)
plot(x = x,y=(logit.mul.pi.matrix [,5]/(-3.485))[x],pch = 20,xlab = "simulation point",ylab = "beta0",main = "N(0,pi2/(3*(p+1)))")
abline(h = 1,col = "blue",lwd = 3)

# Draw the posterior distribution
# draw the prediction graph of the mle method
fitted.prob.burn = predict(burn.model,type = "response")
set.seed(1)
#random1 = sample(1:1000,size= 100,replace = FALSE)
#burn_age = aplore3::burn1000$age[random1]
#fitted.prob.burn = fitted.prob.burn[random1]
plot(sort(burn_age),fitted.prob.burn[order(burn_age)],type = 'l')
sample(1:1000,size= 100,replace = FALSE)


# MLE
# The model is death ~ age+tbsa+race+inh_inj
burn.model=glm(formula = survival~age.scale+tbsa.scale+race+inh_inj,family = binomial(link = 'logit'),data = df.flame)
summary(burn.model)
# Coefficients:
#                  Estimate Std. Error z value Pr(>|z|)    
# (Intercept)      -3.4850     0.2732 -12.757  < 2e-16 ***
#   age.scale       2.0812     0.2091   9.954  < 2e-16 ***
#   tbsa.scale      1.7247     0.1733   9.951  < 2e-16 ***
#   race.scale      0.3398     0.1629   2.086    0.037 *  
#   inh_inj.scale   1.1096     0.2559   4.337 1.45e-05 ***
alpha = -4.182
beta1 = 2.081
beta2 = 1.725
beta3 = 0.624
beta4 = 1.523
# The estimated values of MLE
par(mfrow = c(2,2))
L0 = alpha+beta1*input.burn$age.scale+beta2*input.burn$tbsa.scale+beta3*input.burn$race+beta4*input.burn$inh_inj
E0 = expit(L0)
d0 = data.frame(x = E0)
P0 = ggplot(d0,aes(x))+
  geom_histogram(aes(y=..density..),bins = 30)+ # scale histogram y
  geom_density(col = "blue")+
  ggtitle("MLE")+ xlab("Survival Rate")

# N(0,10e-6)
alpha = -4.221
beta1 = 2.103
beta2 = 1.756
beta3 = 0.631
beta4 = 1.533
E1 = expit(alpha+beta1*input.burn$age.scale+beta2*input.burn$tbsa.scale+beta3*input.burn$race+beta4*input.burn$inh_inj)
d1 = data.frame(x = E1)
P1 = ggplot(d1,aes(x))+
  geom_histogram(aes(y=..density..),bins = 30)+ # scale histogram y
  geom_density(col = "blue")+
  ggtitle("N(0,10e-6)")+ xlab("Survival Rate")


# Draw the comparing scatter plot of x and y
#normal.matrix = Bayesian.burn.normal$BUGSoutput$sims.matrix[,c(5,1,3,4,2)]
#compare.y.normal = c()
#compare.x  = c()
#for (i in 1:1000){
#  compare.x  = c(compare.x,rep(L0[i],1200))
#  compare.yi = normal.matrix[,1]+normal.matrix[,2]*df.flame$age.scale[i]+
#    normal.matrix[,3]*df.flame$tbsa.scale[i]+normal.matrix[,4]*df.flame$race[i]+normal.matrix[,5]*df.flame$inh_inj[i]
#  compare.y.normal = c(compare.y.normal,compare.yi)
#str(compare.x)
#str(compare.y.normal)
#plot(x = compare.x ,y = compare.y.normal,pch = 20)


# Cauchy
alpha = -4.127
beta1 = 2.062
beta2 = 1.737
beta3 = 0.568
beta4 = 1.457
E2 = expit(alpha+beta1*input.burn$age.scale+beta2*input.burn$tbsa.scale+beta3*input.burn$race+beta4*input.burn$inh_inj)
d2 = data.frame(x = E2)
P2 = ggplot(d2,aes(x))+
  geom_histogram(aes(y=..density..),bins = 30)+ # scale histogram y
  geom_density(col = "blue")+
  ggtitle("Cauchy Priors")+ xlab("Survival Rate")



# Christensen et al.
alpha = -3.081
beta1 = 1.741
beta2 = 1.534
beta3 = 0.256
beta4 = 0.941
E3 = expit(alpha+beta1*input.burn$age.scale+beta2*input.burn$tbsa.scale+beta3*input.burn$race+beta4*input.burn$inh_inj)
L3 = alpha+beta1*df.flame$age.scale+beta2*df.flame$tbsa.scale+beta3*df.flame$race+beta4*df.flame$inh_inj
d3 = data.frame(x = E3)
P3 = ggplot(d3,aes(x))+
  geom_histogram(aes(y=..density..),bins = 30)+ # scale histogram y
  geom_density(col = "blue")+
  ggtitle("Christensen Methods")+ xlab("Survival Rate")


# N(0,pi2/(3*(p+1)))
alpha = -3.472
beta1 = 1.714
beta2 = 1.569
beta3 = 0.299
beta4 = 1.149
E4 = expit(alpha+beta1*input.burn$age.scale+beta2*input.burn$tbsa.scale+beta3*input.burn$race+beta4*input.burn$inh_inj)
L4 = alpha+beta1*df.flame$age.scale+beta2*df.flame$tbsa.scale+beta3*df.flame$race+beta4*df.flame$inh_inj
d4 = data.frame(x = E4)
P4 = ggplot(d4,aes(x))+
  geom_histogram(aes(y=..density..),bins = 30)+ # scale histogram y
  geom_density(col = "blue")+
  ggtitle("N(0,pi^2/(3*(p+1)))")+ xlab("Survival Rate")

arrange(P0,P1,P2,P3,P4)

par(mfrow = c(2,3))

hist(E0,xlab = expression(pi(x)),main = "MLE")
hist(E1,xlab = expression(pi(x)),main = expression(N(0,10^6)))
hist(E2,xlab = expression(pi(x)),main = "Cauchy priors")
hist(E3,xlab = expression(pi(x)),main = "Christensen Method")
hist(E4,xlab = expression(pi(x)),main = expression(N(0,pi^2/(3*(p+1)))))

# Draw the posterior results compared with MLE
par(mfrow=c(2,2))
plot(E0,E1,pch = 20,xlab = "MLE",ylab = "N(0,10e-6)")
abline(a=0,b=1,col = 'blue')
plot(E0,E2,pch = 20,xlab = "MLE",ylab = "Cauchy")
abline(a=0,b=1,col = 'blue')
plot(E0,E3,pch = 20,xlab = "MLE",ylab = "(0,S^2/b)")
abline(a=0,b=1,col = 'blue')
plot(E0,E4,pch = 20,xlab = "MLE",ylab = "N(0,pi2/(3*(p+1)))")
abline(a=0,b=1,col = 'blue')

hist(x = df.flame$age,y = L0)

cauchy.matrix = Bayesian.burn.cauchy$BUGSoutput$sims.matrix[,c(5,1,4,3,2)]
posterior_point = matrix(NA,nrow = 1200,ncol = 1000)


# compare the posterior value with MLE results
dim(normal.matrix)  # [1]  999 1006
y = as.vector(t(normal.matrix))
length(y)

compare.y.normal = c()

# 4.1 Create list for the estimation of MLE

y <- as.vector(t(normal.matrix))

par(mfrow =c(2,2))
# 4.2.1 Create list for the estimation of normal priors
logit.mul.normal.matrix
logit.mul.cauchy.matrix
logit.mul.chris.matrix
logit.mul.pi.matrix 

dim(logit.mul.normal.matrix)

y <- as.vector(t(logit.mul.normal.matrix))
length(y)
x <- rep(E0,999)
set.seed(1)
subset.values <- sort(sample(1:999000,size=10000,replace=FALSE))
x.subset <- x[subset.values]
y.subset <- y[subset.values]
plot(x.subset,y.subset,pch = 20,main = "N(0,10^6)")
abline(0,1,col="red")

# 4.2.2 Create list for the estimation of normal priors in Christensen methods
Bayes.output.Chris <-Bayesian.burn.chris1$BUGSoutput$sims.matrix
dim(Bayes.output.Chris) # [1]  999 1006
posterior.survival.matrix.Chris <- Bayes.output.Chris[,7:1006]# 999 1000
post.mean.pred.Chris <- apply(posterior.survival.matrix.Chris,2,mean)

compare.y.Chris = c()
# Create list for the estimation of MLE
compare.x  = c()

y <- as.vector(t(posterior.survival.matrix.Chris))
length(y)
x <- rep(E0,999)
set.seed(1)
subset.values <- sort(sample(1:999000,size=10000,replace=FALSE))
x.subset <- x[subset.values]
y.subset <- y[subset.values]
plot(x.subset,y.subset,pch = 20,main = "N(0,S^2/b)")
abline(0,1,col="red")

# 4.2.3 Create list for the estimation of normal priors in pi method

Bayes.output.pi <-Bayesian.burn.pi$BUGSoutput$sims.matrix
dim(Bayes.output.pi) # [1]  999 1006
posterior.survival.matrix.pi <- Bayes.output.pi[,7:1006]# 999 1000
post.mean.pred.pi <- apply(posterior.survival.matrix.pi,2,mean)

compare.y.pi = c()
# Create list for the estimation of MLE
compare.x  = c()

y <- as.vector(t(posterior.survival.matrix.pi))
length(y)
x <- rep(E0,999)
set.seed(1)
subset.values <- sort(sample(1:999000,size=10000,replace=FALSE))
x.subset <- x[subset.values]
y.subset <- y[subset.values]
plot(x.subset,y.subset,pch = 20,main = "N(0,pi^2/3*(p+1))")
abline(0,1,col="red")

# 4.2.4 Create list for the estimation of normal priors in Cauchy method
Bayes.output.cauchy <-Bayesian.burn.cauchy$BUGSoutput$sims.matrix
dim(Bayes.output.cauchy) # [1]  999 1006
posterior.survival.matrix.cauchy <- Bayes.output.cauchy[,7:1006]# 999 1000
post.mean.pred.cauchy <- apply(posterior.survival.matrix.cauchy,2,mean)

compare.y.cauchy = c()
# Create list for the estimation of MLE
compare.x  = c()

y <- as.vector(t(posterior.survival.matrix.cauchy))
length(y)
x <- rep(E0,999)
set.seed(1)
subset.values <- sort(sample(1:999000,size=10000,replace=FALSE))
x.subset <- x[subset.values]
y.subset <- y[subset.values]
plot(x.subset,y.subset,pch = 20,main = "Cauchy priors")
abline(0,1,col="red")


# 4.5 Different Size
# Size = 1000
# Size = 500
# Size = 250
# Size = 100
# Size = 50
set.seed(1)
sample500 = sample(1:1000,500,replace = FALSE)
sample250 = sample(1:1000,250,replace = FALSE)
sample100 = sample(1:1000,100,replace = FALSE)
sample50 = sample(1:1000,50,replace = FALSE)

# Wide Normal prior
beta.sd = 10^3
beta.prec = 1/beta.sd^2
burn.normal500=wide.mul.model
cat(burn.normal500, file='D:/Study/Postgraduate/paper/burn.normal500.txt')
library(R2jags)
# create the input dataset

input.databurn.normal500 = list("death" = input.burn$death[sample500],
                             "age.scale" = input.burn$age.scale[sample500],
                             "tbsa.scale" = input.burn$tbsa.scale[sample500],
                             "race" = input.burn$race[sample500],
                             "inh_inj" = input.burn$inh_inj[sample500],
                             "n"=500,
                             beta.prec = beta.prec)

Bayesian.burn.normal500 <- jags.model(data=input.databurn.normal500,file='D:/Study/Postgraduate/paper/burn.normal500.txt',
                             inits = wide.mul.inits,
                             n.chains=num.chains, n.adapt=adapt.length)
# use the same parameters as wide.mul.params
normal500.mul.samples = coda.samples(model = Bayesian.burn.normal500,
                                     variable.names = wide.mul.params,
                                     n.iter = burnin+chain.length)
normal500.mul.matrix = rbind(normal500.mul.samples[[1]],normal500.mul.samples[[2]],normal500.mul.samples[[3]])
normal500.mul.results = summary(window(normal500.mul.samples,start=(burnin+1)))
print(normal500.mul.results$statistics) 
#             Mean           SD     Naive SE Time-series SE
# b.age      2.0566199019 0.2854707747 2.330859e-03   1.401770e-02
# b.inh_inj  1.1285729643 0.4691797623 3.830837e-03   6.832236e-03
# b.race     0.5866785068 0.4253054205 3.472604e-03   8.802580e-03
# b.tbsa     1.6044370146 0.2253732038 1.840165e-03   7.085706e-03
# b0        -4.1615911369 0.4759024715 3.885727e-03   2.320148e-02


# 250
# Wide Normal prior
beta.sd = 10^3
beta.prec = 1/beta.sd^2
burn.normal250=wide.mul.model
cat(burn.normal250, file='D:/Study/Postgraduate/paper/burn.normal250.txt')
library(R2jags)
# create the input dataset
input.databurn.normal250 = list("death" = input.burn$death[sample250],
                                "age.scale" = input.burn$age.scale[sample250],
                                "tbsa.scale" = input.burn$tbsa.scale[sample250],
                                "race" = input.burn$race[sample250],
                                "inh_inj" = input.burn$inh_inj[sample250],
                                "n"=250,
                                beta.prec = beta.prec)

Bayesian.burn.normal250 <- jags.model(data=input.databurn.normal250,file='D:/Study/Postgraduate/paper/burn.normal250.txt',
                                      inits = wide.mul.inits,
                                      n.chains=num.chains, n.adapt=adapt.length)
# use the same parameters as wide.mul.params
normal250.mul.samples = coda.samples(model = Bayesian.burn.normal250,
                                     variable.names = wide.mul.params,
                                     n.iter = burnin+chain.length)
normal250.mul.matrix = rbind(normal250.mul.samples[[1]],normal250.mul.samples[[2]],normal250.mul.samples[[3]])
normal250.mul.results = summary(window(normal250.mul.samples,start=(burnin+1)))
print(normal250.mul.results$statistics)
#                Mean           SD     Naive SE Time-series SE
# b.age      3.429077e+00 8.625827e-01 7.042958e-03   9.379440e-02
# b.inh_inj  1.835615e+00 1.225788e+00 1.000852e-02   2.604867e-02
# b.race     9.134729e-01 8.800601e-01 7.185660e-03   2.731707e-02
# b.tbsa     3.801662e+00 8.578766e-01 7.004533e-03   6.857052e-02
# b0        -6.550494e+00 1.507342e+00 1.230739e-02   1.476953e-01

# 100
# Wide Normal prior
beta.sd = 10^3
beta.prec = 1/beta.sd^2
burn.normal100=wide.mul.model
cat(burn.normal100, file='D:/Study/Postgraduate/paper/burn.normal100.txt')
library(R2jags)
# create the input dataset
input.databurn.normal100 = list("death" = input.burn$death[sample100],
                                "age.scale" = input.burn$age.scale[sample100],
                                "tbsa.scale" = input.burn$tbsa.scale[sample100],
                                "race" = input.burn$race[sample100],
                                "inh_inj" = input.burn$inh_inj[sample100],
                                "n"=100,
                                beta.prec = beta.prec)

Bayesian.burn.normal100 <- jags.model(data=input.databurn.normal100,file='D:/Study/Postgraduate/paper/burn.normal100.txt',
                                      inits = wide.mul.inits,
                                      n.chains=num.chains, n.adapt=adapt.length)
# use the same parameters as wide.mul.params
normal100.mul.samples = coda.samples(model = Bayesian.burn.normal100,
                                     variable.names = wide.mul.params,
                                     n.iter = burnin+chain.length)
normal100.mul.matrix = rbind(normal100.mul.samples[[1]],normal100.mul.samples[[2]],normal100.mul.samples[[3]])
normal100.mul.results = summary(window(normal100.mul.samples,start=(burnin+1)))
print(normal100.mul.results$statistics) 
#                Mean           SD     Naive SE Time-series SE
# b.age      2.2055521374 0.8987237169 7.338048e-03   5.127031e-02
# b.inh_inj  0.0681154567 2.2167764103 1.809990e-02   6.356331e-02
# b.race    -3.3541361172 2.7455133259 2.241702e-02   1.795406e-01
# b.tbsa     5.4067834791 1.6022198417 1.308207e-02   1.057605e-01
# b0        -3.7798804482 1.3588705020 1.109513e-02   7.574644e-02


# n=50
# Wide Normal prior
beta.sd = 10^3
beta.prec = 1/beta.sd^2
burn.normal50=wide.mul.model
cat(burn.normal50, file='D:/Study/Postgraduate/paper/burn.normal50.txt')
library(R2jags)
# create the input dataset
input.databurn.normal50 = list("death" = input.burn$death[sample50],
                                "age.scale" = input.burn$age.scale[sample50],
                                "tbsa.scale" = input.burn$tbsa.scale[sample50],
                                "race" = input.burn$race[sample50],
                                "inh_inj" = input.burn$inh_inj[sample50],
                                "n"=50,
                                beta.prec = beta.prec)

Bayesian.burn.normal50 <- jags.model(data=input.databurn.normal50,file='D:/Study/Postgraduate/paper/burn.normal50.txt',
                                      inits = wide.mul.inits,
                                      n.chains=num.chains, n.adapt=adapt.length)
# use the same parameters as wide.mul.params
normal50.mul.samples = coda.samples(model = Bayesian.burn.normal50,
                                     variable.names = wide.mul.params,
                                     n.iter = burnin+chain.length)
normal50.mul.matrix = rbind(normal50.mul.samples[[1]],normal50.mul.samples[[2]],normal50.mul.samples[[3]])
normal50.mul.results = summary(window(normal50.mul.samples,start=(burnin+1)))
print(normal50.mul.results$statistics) 
#                    Mean                SD           Naive SE     Time-series SE
# b.age      26.8259523029885187 11.50854194926967 0.0939668515304199 2.3615412992198821
# b.inh_inj  -0.3345402538278647  8.52130795690281 0.0695761881184337 1.0756891287652546
# b.race      3.7233195399728687  8.54839826529455 0.0697973795602151 1.4813411981262867
# b.tbsa     35.0509830081063782 14.46415047385610 0.1180992940792766 3.1671001130776557
# b0        -32.1052420045691065 15.51368017810806 0.1266686682303147 3.2898506744219556

# Christensen distribution
beta.sd = sqrt(1/3)
beta.prec = 1/beta.sd^2

burn.chris500 = wide.mul.model

burn.chris.inits <- list()
for(i in 1:num.chains) {
  b0       <- rnorm(n=1,mean=0,sd=beta.prec)
  b.age    <- rnorm(n=1,mean=0,sd=beta.prec)
  b.tbsa   <- rnorm(n=1,mean=0,sd=beta.prec)
  b.race   <- rnorm(n=1,mean=0,sd=beta.prec)
  b.inh_inj <- rnorm(n=1,mean=0,sd=beta.prec)
  burn.chris.inits[[i]] <- list(b0=b0,b.age=b.age,b.tbsa=b.tbsa,
                                b.race=b.race,b.inh_inj=b.inh_inj)
}


cat(burn.chris500, file='D:/Study/Postgraduate/paper/burn.chris500.txt')
library(R2jags)
# create the input dataset

input.databurn.chris500 = list("death" = input.burn$death[sample500],
                                "age.scale" = input.burn$age.scale[sample500],
                                "tbsa.scale" = input.burn$tbsa.scale[sample500],
                                "race" = input.burn$race[sample500],
                                "inh_inj" = input.burn$inh_inj[sample500],
                                "n"=500,
                                beta.prec = beta.prec)
Bayesian.burn.chris500 <- jags.model(file = "D:/Study/Postgraduate/paper/burn.chris500.txt",
                                data=input.databurn.chris500, inits=burn.chris.inits, n.chains=num.chains,
                                n.adapt=adapt.length)
chris500.mul.params = wide.mul.params
# use the same parameters as wide.mul.params
chris500.mul.samples = coda.samples(model = Bayesian.burn.chris500,
                                     variable.names = chris500.mul.params,
                                     n.iter = burnin+chain.length)
chris500.mul.matrix = rbind(chris500.mul.samples[[1]],chris500.mul.samples[[2]],chris500.mul.samples[[3]])
chris500.mul.results = summary(window(chris500.mul.samples,start=(burnin+1)))
print(chris500.mul.results$statistics) 
#                Mean          SD     Naive SE Time-series SE
# b.age      1.300340048 0.167848728 0.00137047913  0.00365805067
# b.inh_inj  0.604722829 0.334938180 0.00273475879  0.00359609015
# b.race    -0.015230531 0.296651420 0.00242214870  0.00441526303
# b.tbsa     1.231691119 0.165284380 0.00134954131  0.00274972665
# b0        -2.707990570 0.233727422 0.00190837641  0.00495824560



# n = 250
beta.sd = sqrt(1/3)
beta.prec = 1/beta.sd^2
burn.chris250=wide.mul.model
cat(burn.chris250, file='D:/Study/Postgraduate/paper/burn.chris250.txt')
library(R2jags)
# create the input dataset

input.databurn.chris250 = list("death" = input.burn$death[sample250],
                               "age.scale" = input.burn$age.scale[sample250],
                               "tbsa.scale" = input.burn$tbsa.scale[sample250],
                               "race" = input.burn$race[sample250],
                               "inh_inj" = input.burn$inh_inj[sample250],
                               "n"=250,
                               beta.prec = beta.prec)

Bayesian.burn.chris250 <- jags.model(data=input.databurn.chris250,file='D:/Study/Postgraduate/paper/burn.chris250.txt',
                                     inits = burn.chris.inits,
                                     n.chains=num.chains, n.adapt=adapt.length)
# use the same parameters as wide.mul.params
chris250.mul.samples = coda.samples(model = Bayesian.burn.chris250,
                                    variable.names = wide.mul.params,
                                    n.iter = burnin+chain.length)
chris250.mul.matrix = rbind(chris250.mul.samples[[1]],chris250.mul.samples[[2]],chris250.mul.samples[[3]])
chris250.mul.results = summary(window(chris250.mul.samples,start=(burnin+1)))
print(chris250.mul.results$statistics) 
#                Mean          SD      Naive SE Time-series SE
# b.age      0.995547975 0.219706045 0.00179389234  0.00418759732
# b.inh_inj  0.398082775 0.474963517 0.00387806088  0.00438518338
# b.race    -0.397703304 0.368286004 0.00300704263  0.00496580632
# b.tbsa     1.664726993 0.271627088 0.00221782589  0.00470751344
# b0        -2.235964976 0.289270027 0.00236187988  0.00526236812

# n = 100
beta.sd = sqrt(1/3)
beta.prec = 1/beta.sd^2
burn.chris100=wide.mul.model
cat(burn.chris100, file='D:/Study/Postgraduate/paper/burn.chris100.txt')
library(R2jags)
# create the input dataset

input.databurn.chris100 = list("death" = input.burn$death[sample100],
                               "age.scale" = input.burn$age.scale[sample100],
                               "tbsa.scale" = input.burn$tbsa.scale[sample100],
                               "race" = input.burn$race[sample100],
                               "inh_inj" = input.burn$inh_inj[sample100],
                               "n"=100,
                               beta.prec = beta.prec)

Bayesian.burn.chris100 <- jags.model(data=input.databurn.chris100,file='D:/Study/Postgraduate/paper/burn.chris100.txt',
                                     inits = burn.chris.inits,
                                     n.chains=num.chains, n.adapt=adapt.length)
# use the same parameters as wide.mul.params
chris100.mul.samples = coda.samples(model = Bayesian.burn.chris100,
                                    variable.names = wide.mul.params,
                                    n.iter = burnin+chain.length)
chris100.mul.matrix = rbind(chris100.mul.samples[[1]],chris100.mul.samples[[2]],chris100.mul.samples[[3]])
chris100.mul.results = summary(window(chris100.mul.samples,start=(burnin+1)))
print(chris100.mul.results$statistics) 
#               Mean         SD     Naive SE Time-series SE
# b.age      0.61068429 0.27578823 0.0022518015   0.0036673244
# b.inh_inj  0.02742337 0.51355305 0.0041931431   0.0045560871
# b.race    -0.74804328 0.44139016 0.0036039355   0.0050023165
# b.tbsa     1.58063571 0.37387438 0.0030526715   0.0047279686
# b0        -1.44099347 0.32511914 0.0026545867   0.0040988417

# n = 50
beta.sd = sqrt(1/3)
beta.prec = 1/beta.sd^2
burn.chris50=wide.mul.model
cat(burn.chris50, file='D:/Study/Postgraduate/paper/burn.chris50.txt')
library(R2jags)

# create the input dataset
input.databurn.chris50 = list("death" = input.burn$death[sample50],
                               "age.scale" = input.burn$age.scale[sample50],
                               "tbsa.scale" = input.burn$tbsa.scale[sample50],
                               "race" = input.burn$race[sample50],
                               "inh_inj" = input.burn$inh_inj[sample50],
                               "n"=50,
                               beta.prec = beta.prec)

Bayesian.burn.chris50 <- jags.model(data=input.databurn.chris50,file='D:/Study/Postgraduate/paper/burn.chris50.txt',
                                     inits = wide.mul.inits,
                                     n.chains=num.chains, n.adapt=adapt.length)
# use the same parameters as wide.mul.params
chris50.mul.samples = coda.samples(model = Bayesian.burn.chris50,
                                    variable.names = wide.mul.params,
                                    n.iter = burnin+chain.length)
chris50.mul.matrix = rbind(chris50.mul.samples[[1]],chris50.mul.samples[[2]],chris50.mul.samples[[3]])
chris50.mul.results = summary(window(chris50.mul.samples,start=(burnin+1)))
print(chris50.mul.results$statistics) 
#               Mean         SD     Naive SE Time-series SE
# b.age      0.43780673 0.36087398 0.0029465237   0.0040808816
# b.inh_inj -0.02388815 0.50313863 0.0041081097   0.0042504830
# b.race    -0.28930174 0.45525784 0.0037171647   0.0043125640
# b.tbsa     1.05548396 0.38014410 0.0031038636   0.0039815648
# b0        -1.39203804 0.36365578 0.0029692370   0.0037678514


# Pi distribution
beta.sd = sqrt((pi^2)/15)
beta.prec = 1/beta.sd^2
burn.pi500=wide.mul.model
cat(burn.pi500, file='D:/Study/Postgraduate/paper/burn.pi500.txt')
library(R2jags)
# create the input dataset
input.databurn.pi500 = list("death" = input.burn$death[sample500],
                               "age.scale" = input.burn$age.scale[sample500],
                               "tbsa.scale" = input.burn$tbsa.scale[sample500],
                               "race" = input.burn$race[sample500],
                               "inh_inj" = input.burn$inh_inj[sample500],
                               "n"=500,
                               beta.prec = beta.prec)

Bayesian.burn.pi500 <- jags.model(data=input.databurn.pi500,file='D:/Study/Postgraduate/paper/burn.pi500.txt',
                                     inits = wide.mul.inits,
                                     n.chains=num.chains, n.adapt=adapt.length)
# use the same parameters as wide.mul.params
pi500.mul.samples = coda.samples(model = Bayesian.burn.pi500,
                                    variable.names = wide.mul.params,
                                    n.iter = burnin+chain.length)
pi500.mul.matrix = rbind(pi500.mul.samples[[1]],pi500.mul.samples[[2]],pi500.mul.samples[[3]])
pi500.mul.results = summary(window(pi500.mul.samples,start=(burnin+1)))
print(pi500.mul.results$statistics) 
#                Mean          SD     Naive SE Time-series SE
# b.age      1.523628870 0.196098586 1.601138e-03   5.316004e-03
# b.inh_inj  0.765880970 0.377709758 3.083987e-03   4.466222e-03
# b.race     0.149584509 0.337302548 2.754064e-03   5.691874e-03
# b.tbsa     1.343829101 0.178376071 1.456435e-03   3.422582e-03
# b0        -3.131094960 0.293326826 2.395004e-03   8.406518e-03

#n = 250
beta.sd = sqrt((pi^2)/15)
beta.prec = 1/beta.sd^2
burn.pi250=wide.mul.model
cat(burn.pi250, file='D:/Study/Postgraduate/paper/burn.pi250.txt')
library(R2jags)
# create the input dataset
input.databurn.pi250 = list("death" = input.burn$death[sample250],
                            "age.scale" = input.burn$age.scale[sample250],
                            "tbsa.scale" = input.burn$tbsa.scale[sample250],
                            "race" = input.burn$race[sample250],
                            "inh_inj" = input.burn$inh_inj[sample250],
                            "n"=250,
                            beta.prec = beta.prec)

Bayesian.burn.pi250 <- jags.model(data=input.databurn.pi250,file='D:/Study/Postgraduate/paper/burn.pi250.txt',
                                  inits = wide.mul.inits,
                                  n.chains=num.chains, n.adapt=adapt.length)
# use the same parameters as wide.mul.params
pi250.mul.samples = coda.samples(model = Bayesian.burn.pi250,
                                 variable.names = wide.mul.params,
                                 n.iter = burnin+chain.length)
pi250.mul.matrix = rbind(pi250.mul.samples[[1]],pi250.mul.samples[[2]],pi250.mul.samples[[3]])
pi250.mul.results = summary(window(pi250.mul.samples,start=(burnin+1)))
print(pi250.mul.results$statistics) 
#                Mean          SD     Naive SE Time-series SE
# b.age      1.275305181 0.272051653 2.221292e-03   7.470244e-03
# b.inh_inj  0.591911573 0.606747338 4.954071e-03   6.205688e-03
# b.race    -0.310329717 0.442551671 3.613419e-03   7.028737e-03
# b.tbsa     1.953937507 0.333977536 2.726915e-03   7.148584e-03
# b0        -2.706858999 0.375473011 3.065724e-03   9.380116e-03

#n = 100
beta.sd = sqrt((pi^2)/15)
beta.prec = 1/beta.sd^2
burn.pi100=wide.mul.model
cat(burn.pi100, file='D:/Study/Postgraduate/paper/burn.pi100.txt')
library(R2jags)
# create the input dataset
input.databurn.pi100 = list("death" = input.burn$death[sample100],
                            "age.scale" = input.burn$age.scale[sample100],
                            "tbsa.scale" = input.burn$tbsa.scale[sample100],
                            "race" = input.burn$race[sample100],
                            "inh_inj" = input.burn$inh_inj[sample100],
                            "n"=100,
                            beta.prec = beta.prec)

Bayesian.burn.pi100 <- jags.model(data=input.databurn.pi100,file='D:/Study/Postgraduate/paper/burn.pi100.txt',
                                  inits = wide.mul.inits,
                                  n.chains=num.chains, n.adapt=adapt.length)
# use the same parameters as wide.mul.params
pi100.mul.samples = coda.samples(model = Bayesian.burn.pi100,
                                 variable.names = wide.mul.params,
                                 n.iter = burnin+chain.length)
pi100.mul.matrix = rbind(pi100.mul.samples[[1]],pi100.mul.samples[[2]],pi100.mul.samples[[3]])
pi100.mul.results = summary(window(pi100.mul.samples,start=(burnin+1)))
print(pi100.mul.results$statistics) 
#                Mean          SD     Naive SE Time-series SE
# b.age      0.801742033 0.337411608 2.754954e-03   0.0055003886
# b.inh_inj  0.028488007 0.682539854 5.572915e-03   0.0063033206
# b.race    -0.957006751 0.590876576 4.824487e-03   0.0085614428
# b.tbsa     2.042832516 0.478004579 3.902891e-03   0.0077221324
# b0        -1.720103790 0.411360548 3.358745e-03   0.0064025924

#n = 50
beta.sd = sqrt((pi^2)/15)
beta.prec = 1/beta.sd^2
burn.pi50=wide.mul.model
cat(burn.pi50, file='D:/Study/Postgraduate/paper/burn.pi50.txt')
library(R2jags)
# create the input dataset
input.databurn.pi50 = list("death" = input.burn$death[sample50],
                            "age.scale" = input.burn$age.scale[sample50],
                            "tbsa.scale" = input.burn$tbsa.scale[sample50],
                            "race" = input.burn$race[sample50],
                            "inh_inj" = input.burn$inh_inj[sample50],
                            "n"=50,
                            beta.prec = beta.prec)

Bayesian.burn.pi50 <- jags.model(data=input.databurn.pi50,file='D:/Study/Postgraduate/paper/burn.pi50.txt',
                                  inits = wide.mul.inits,
                                  n.chains=num.chains, n.adapt=adapt.length)
# use the same parameters as wide.mul.params
pi50.mul.samples = coda.samples(model = Bayesian.burn.pi50,
                                 variable.names = wide.mul.params,
                                 n.iter = burnin+chain.length)
pi50.mul.matrix = rbind(pi50.mul.samples[[1]],pi50.mul.samples[[2]],pi50.mul.samples[[3]])
pi50.mul.results = summary(window(pi50.mul.samples,start=(burnin+1)))
print(pi50.mul.results$statistics) 
#               Mean         SD     Naive SE Time-series SE
# b.age      0.65620154 0.45063937 0.0036794550   0.0061519574
# b.inh_inj -0.01267657 0.66985449 0.0054693390   0.0059275048
# b.race    -0.25063828 0.59680074 0.0048728577   0.0065583777
# b.tbsa     1.42035747 0.49862135 0.0040712263   0.0064639875
# b0        -1.73868625 0.47339072 0.0038652190   0.0056186038


# Cauchy distribution
#--- Data file ---
cauchy500.mul.data <- list(n=500,death=input.burn$death[sample500],age.scale=input.burn$age.scale[sample500],tbsa.scale=input.burn$tbsa.scale[sample500],
                        race= input.burn$race[sample500],inh_inj=input.burn$inh_inj[sample500])

set.seed(1)

input.databurn.cauchy500 = list("death" = input.burn$death[sample500],
                                "age.scale" = input.burn$age.scale[sample500],
                                "tbsa.scale" = input.burn$tbsa.scale[sample500],
                                "race" = input.burn$race[sample500],
                                "inh_inj" = input.burn$inh_inj[sample500],
                                "n"=500)
cauchy500.mul.inits <- list()
for(i in 1:num.chains) {
  # get "incompatable" initial values...?
  # b0 <- rnorm(n=1,mean=0,sd=beta.sd)
  # b1 <- rnorm(n=1,mean=0,sd=beta.sd)
  b0       = rnorm(n=1,0,1)
  b.age    = rnorm(n=1,0,1)
  b.tbsa   = rnorm(n=1,0,1)
  b.race   = rnorm(n=1,0,1)
  b.inh_inj = rnorm(n=1,0,1)
  cauchy500.mul.inits[[i]] <- list(b0=b0,b.age=b.age,b.tbsa=b.tbsa,
                                b.race= b.race,b.inh_inj=b.inh_inj)
}

cauchy500.mul.model = cauchy.mul.model

cat(cauchy500.mul.model,file="D:/Study/Postgraduate/paper/cauchy500_mul_model.txt")

#-- run jags
# Initialize the model
cauchy500.mul.init.JAGS <- jags.model( file = "D:/Study/Postgraduate/paper/cauchy500_mul_model.txt",
                                    data=cauchy500.mul.data, inits=cauchy500.mul.inits, n.chains=num.chains,
                                    n.adapt=adapt.length)

# specify parameters to monitor
cauchy500.mul.params <- c("b0","b.age",'b.tbsa',"b.race",'b.inh_inj','p')

# Run JAGS, using coda.samples a wrapper for jags.samples
cauchy500.mul.samples <- coda.samples(model=cauchy500.mul.init.JAGS,
                                   variable.names=cauchy500.mul.params,
                                   n.iter=burnin+chain.length)
logit.mul.cauchy500.matrix = rbind(cauchy500.mul.samples[[1]],cauchy500.mul.samples[[2]],cauchy500.mul.samples[[3]])

#plot(cauchy.mul.samples)
cauchy500.mul.results <- summary(window(cauchy500.mul.samples,start=(burnin+1)))
print(cauchy500.mul.results$statistics)
#                  Mean           SD     Naive SE Time-series SE
# b.age      1.9843813054 0.2700632418 0.002012932556  0.00670984377
# b.inh_inj  1.0422778285 0.4564696201 0.003402323667  0.00598260385
# b.race     0.5108664093 0.4041292297 0.003012201431  0.00644704095
# b.tbsa     1.5764257469 0.2187044945 0.001630127055  0.00369033285
# b0        -4.0111264983 0.4434487942 0.003305272161  0.01191818130
sd1 = sd(input.burn$age.scale[sample500])
sd2 = sd(input.burn$tbsa.scale[sample500])


#n = 250
cauchy250.mul.data <- list(n=250,death=input.burn$death[sample250],age.scale=input.burn$age.scale[sample250],tbsa.scale=input.burn$tbsa.scale[sample250],
                           race= input.burn$race[sample250],inh_inj=input.burn$inh_inj[sample250])

set.seed(1)
cauchy250.mul.inits<- list()
for(i in 1:num.chains) {
  # get "incompatable" initial values...?
  # b0 <- rnorm(n=1,mean=0,sd=beta.sd)
  # b1 <- rnorm(n=1,mean=0,sd=beta.sd)
  b0       = rnorm(n=1,0,1)
  b.age    = rnorm(n=1,0,1)
  b.tbsa   = rnorm(n=1,0,1)
  b.race   = rnorm(n=1,0,1)
  b.inh_inj = rnorm(n=1,0,1)
  cauchy250.mul.inits[[i]] <- list(b0=b0,b.age=b.age,b.tbsa=b.tbsa,
                                   b.race= b.race,b.inh_inj=b.inh_inj)
}

cauchy250.mul.model = cauchy.mul.model

cat(cauchy250.mul.model,file="D:/Study/Postgraduate/paper/cauchy250_mul_model.txt")

#-- run jags
# Initialize the model
cauchy250.mul.init.JAGS <- jags.model( file = "D:/Study/Postgraduate/paper/cauchy500_mul_model.txt",
                                       data=cauchy250.mul.data, inits=cauchy250.mul.inits, n.chains=num.chains,
                                       n.adapt=adapt.length)

# specify parameters to monitor
cauchy250.mul.params <- c("b0","b.age",'b.tbsa',"b.race",'b.inh_inj','p')

# Run JAGS, using coda.samples a wrapper for jags.samples
cauchy250.mul.samples <- coda.samples(model=cauchy250.mul.init.JAGS,
                                      variable.names=cauchy250.mul.params,
                                      n.iter=burnin+chain.length)
logit.mul.cauchy250.matrix = rbind(cauchy250.mul.samples[[1]],cauchy250.mul.samples[[2]],cauchy250.mul.samples[[3]])

#plot(cauchy.mul.samples)
cauchy250.mul.results <- summary(window(cauchy250.mul.samples,start=(burnin+1)))
print(cauchy250.mul.results$statistics)

#                 Mean           SD       Naive SE Time-series SE
# b.age      2.71353078997 0.6461507636 0.004816123437 0.023645296131
# b.inh_inj  1.25258383518 1.0071030142 0.007506502667 0.011679512669
# b.race     0.49281223976 0.7092903177 0.005286737887 0.013741024627
# b.tbsa     3.20854265346 0.6679948627 0.004978939739 0.018461826851
# b0        -5.25137043318 1.0801704402 0.008051115106 0.041102688615
sd1 = sd(input.burn$age.scale[sample250])
sd2 = sd(input.burn$tbsa.scale[sample250])

#n=100
cauchy100.mul.data <- list(n=100,death=input.burn$death[sample100],age.scale=input.burn$age.scale[sample100],tbsa.scale=input.burn$tbsa.scale[sample100],
                           race= input.burn$race[sample100],inh_inj=input.burn$inh_inj[sample100])

set.seed(1)
cauchy100.mul.inits<- list()
for(i in 1:num.chains) {
  # get "incompatable" initial values...?
  # b0 <- rnorm(n=1,mean=0,sd=beta.sd)
  # b1 <- rnorm(n=1,mean=0,sd=beta.sd)
  b0       = rnorm(n=1,0,1)
  b.age    = rnorm(n=1,0,1)
  b.tbsa   = rnorm(n=1,0,1)
  b.race   = rnorm(n=1,0,1)
  b.inh_inj = rnorm(n=1,0,1)
  cauchy100.mul.inits[[i]] <- list(b0=b0,b.age=b.age,b.tbsa=b.tbsa,
                                   b.race= b.race,b.inh_inj=b.inh_inj)
}

cauchy100.mul.model = cauchy.mul.model

cat(cauchy100.mul.model,file="D:/Study/Postgraduate/paper/cauchy100_mul_model.txt")

#-- run jags
# Initialize the model
cauchy100.mul.init.JAGS <- jags.model( file = "D:/Study/Postgraduate/paper/cauchy100_mul_model.txt",
                                       data=cauchy100.mul.data, inits=cauchy100.mul.inits, n.chains=num.chains,
                                       n.adapt=adapt.length)

# specify parameters to monitor
cauchy100.mul.params <- c("b0","b.age",'b.tbsa',"b.race",'b.inh_inj','p')

# Run JAGS, using coda.samples a wrapper for jags.samples
cauchy100.mul.samples <- coda.samples(model=cauchy100.mul.init.JAGS,
                                      variable.names=cauchy100.mul.params,
                                      n.iter=burnin+chain.length)
logit.mul.cauchy100.matrix = rbind(cauchy100.mul.samples[[1]],cauchy100.mul.samples[[2]],cauchy100.mul.samples[[3]])

#plot(cauchy.mul.samples)
cauchy100.mul.results <- summary(window(cauchy100.mul.samples,start=(burnin+1)))
print(cauchy100.mul.results$statistics)
#               Mean          SD       Naive SE Time-series SE
# b.age      1.6946278374 0.720854546 0.005372932552  0.01715377655
# b.inh_inj -0.0303732832 1.372477157 0.010229840736  0.01528515281
# b.race    -1.7271968678 1.626755742 0.012125121409  0.02059673214
# b.tbsa     4.2473148916 1.206583241 0.008993340490  0.01600816385
# b0        -3.0404413471 1.050112252 0.007827074598  0.02549678572

# n=50
cauchy50.mul.data <- list(n=50,death=input.burn$death[sample50],age.scale=input.burn$age.scale[sample50],tbsa.scale=input.burn$tbsa.scale[sample50],
                           race= input.burn$race[sample50],inh_inj=input.burn$inh_inj[sample50])

set.seed(1)
cauchy50.mul.inits<- list()
for(i in 1:num.chains) {
  # get "incompatable" initial values...?
  # b0 <- rnorm(n=1,mean=0,sd=beta.sd)
  # b1 <- rnorm(n=1,mean=0,sd=beta.sd)
  b0       = rnorm(n=1,0,1)
  b.age    = rnorm(n=1,0,1)
  b.tbsa   = rnorm(n=1,0,1)
  b.race   = rnorm(n=1,0,1)
  b.inh_inj = rnorm(n=1,0,1)
  cauchy50.mul.inits[[i]] <- list(b0=b0,b.age=b.age,b.tbsa=b.tbsa,
                                   b.race= b.race,b.inh_inj=b.inh_inj)
}

cauchy50.mul.model = cauchy.mul.model

cat(cauchy50.mul.model,file="D:/Study/Postgraduate/paper/cauchy50_mul_model.txt")


#-- run jags
# Initialize the model
cauchy50.mul.init.JAGS <- jags.model( file = "D:/Study/Postgraduate/paper/cauchy50_mul_model.txt",
                                       data=cauchy50.mul.data, inits=cauchy50.mul.inits, n.chains=num.chains,
                                       n.adapt=adapt.length)

# specify parameters to monitor
cauchy50.mul.params <- c("b0","b.age",'b.tbsa',"b.race",'b.inh_inj','p')

# Run JAGS, using coda.samples a wrapper for jags.samples
cauchy50.mul.samples <- coda.samples(model=cauchy50.mul.init.JAGS,
                                      variable.names=cauchy50.mul.params,
                                      n.iter=burnin+chain.length)
logit.mul.cauchy50.matrix = rbind(cauchy50.mul.samples[[1]],cauchy50.mul.samples[[2]],cauchy50.mul.samples[[3]])

#plot(cauchy.mul.samples)
cauchy50.mul.results <- summary(window(cauchy50.mul.samples,start=(burnin+1)))
print(cauchy50.mul.results$statistics)
#                Mean          SD      Naive SE Time-series SE
# b.age      4.1490079365 2.549660738 0.01900404910  0.13497011013
# b.inh_inj  0.0061260088 1.543942910 0.01150787100  0.02384578346
# b.race     0.1853735315 1.456494923 0.01085607219  0.02752522453
# b.tbsa     6.4907654030 3.601482571 0.02684386616  0.17897065527
# b0        -5.3545956360 2.683930094 0.02000483379  0.12924847434




sample1000.p1000 = sample(1:18000000,1000,replace = FALSE)
sample1000.p500 = sample(1:9000000,1000,replace = FALSE)
sample1000.p250 = sample(1:4500000,1000,replace = FALSE)
sample1000.p100 = sample(1:1800000,1000,replace = FALSE)
sample1000.p50 = sample(1:900000,1000,replace = FALSE)

sample1000.c = sample(1:18000,10000,replace = FALSE)

# Probability
logit.mul.normal.matrix.p = as.vector(logit.mul.normal.matrix[,6:1005])[sample1000.p1000]#18000  1005
normal500.mul.matrix.p = as.vector(normal500.mul.matrix[,6:505])[sample1000.p500]#18000   505
normal250.mul.matrix.p = as.vector(normal250.mul.matrix[,6:255])[sample1000.p250]#18000   255
normal100.mul.matrix.p = as.vector(normal100.mul.matrix[,6:105])[sample1000.p100]#18000   105
normal50.mul.matrix.p = as.vector(normal50.mul.matrix[,6:55])[sample1000.p50]# 18000    55

logit.mul.chris.matrix.p = as.vector(logit.mul.chris.matrix[,6:1005])[sample1000.p1000]# dimension is wrong
chris500.mul.matrix.p = as.vector(chris500.mul.matrix[,6:505])[sample1000.p500]#18000   505
chris250.mul.matrix.p = as.vector(chris250.mul.matrix[,6:255])[sample1000.p250]#18000   255
chris100.mul.matrix.p = as.vector(chris100.mul.matrix[,6:105])[sample1000.p100]#18000   105
chris50.mul.matrix.p = as.vector(chris50.mul.matrix[,6:55])[sample1000.p50]#18000    55

logit.mul.pi.matrix.p = as.vector(logit.mul.pi.matrix[,6:1005])[sample1000.p1000]#18000  1005
pi500.mul.matrix.p = as.vector(pi500.mul.matrix[,6:505])[sample1000.p500]#18000   505
pi250.mul.matrix.p = as.vector(pi250.mul.matrix[,6:255])[sample1000.p250]#18000   255
pi100.mul.matrix.p = as.vector(pi100.mul.matrix[,6:105])[sample1000.p100]#18000   105
pi50.mul.matrix.p = as.vector(pi50.mul.matrix[,6:55])[sample1000.p50]#18000    55

logit.mul.cauchy.matrix.p = as.vector(logit.mul.cauchy.matrix[,6:1005])[sample1000.p1000]# 18000  1005
logit.mul.cauchy500.matrix.p = as.vector(logit.mul.cauchy500.matrix[,6:505])[sample1000.p500] # 18000   505
logit.mul.cauchy250.matrix.p = as.vector(logit.mul.cauchy250.matrix[,6:255])[sample1000.p250] # 18000   255
logit.mul.cauchy100.matrix.p = as.vector(logit.mul.cauchy100.matrix[,6:105])[sample1000.p100] # 18000   105
logit.mul.cauchy50.matrix.p = as.vector(logit.mul.cauchy50.matrix[,6:55])[sample1000.p50] # 18000    55

# Coefficient
logit.mul.normal.matrix.c = logit.mul.normal.matrix[,1:5][sample1000.c,] #18000  1005
normal500.mul.matrix.c = normal500.mul.matrix[,1:5][sample1000.c,]#18000   505
normal250.mul.matrix.c = normal250.mul.matrix[,1:5][sample1000.c,]#18000   255
normal100.mul.matrix.c = normal100.mul.matrix[,1:5][sample1000.c,]#18000   105
normal50.mul.matrix.c = normal50.mul.matrix[,1:5][sample1000.c,]# 18000    55



logit.mul.chris.matrix.c = logit.mul.chris.matrix[,1:5][sample1000.c,]# dimension is wrong
chris500.mul.matrix.c = chris500.mul.matrix[,1:5][sample1000.c,]#18000   505
chris250.mul.matrix.c = chris250.mul.matrix[,1:5][sample1000.c,]#18000   255
chris100.mul.matrix.c = chris100.mul.matrix[,1:5][sample1000.c,]#18000   105
chris50.mul.matrix.c = chris50.mul.matrix[,1:5][sample1000.c,]#18000    55

logit.mul.pi.matrix.c = logit.mul.pi.matrix[,1:5][sample1000.c,]#18000  1005
pi500.mul.matrix.c = pi500.mul.matrix[,1:5][sample1000.c,]#18000   505
pi250.mul.matrix.c = pi250.mul.matrix[,1:5][sample1000.c,]#18000   255
pi100.mul.matrix.c = pi100.mul.matrix[,1:5][sample1000.c,]#18000   105
pi50.mul.matrix.c = pi50.mul.matrix[,1:5][sample1000.c,]#18000    55

logit.mul.cauchy.matrix.c = logit.mul.cauchy.matrix[,1:5][sample1000.c,]# 18000  1005
logit.mul.cauchy500.matrix.c = logit.mul.cauchy500.matrix[,1:5][sample1000.c,] # 18000   505
logit.mul.cauchy250.matrix.c = logit.mul.cauchy250.matrix[,1:5][sample1000.c,] # 18000   255
logit.mul.cauchy100.matrix.c = logit.mul.cauchy100.matrix[,1:5][sample1000.c,] # 18000   105
logit.mul.cauchy50.matrix.c = logit.mul.cauchy50.matrix[,1:5][sample1000.c,] # 18000    55

# We sample 1000 dataset from the total of 18000000
# Different sample sizes in the same method
# naive normal + b0
# plot(density(logit.mul.normal.matrix.p))
# lines(density(normal500.mul.matrix.p),col = c("#2166AC"))
# lines(density(normal250.mul.matrix.p),col = c("#4393C3"))
# lines(density(normal100.mul.matrix.p),col = c("#92C5DE"))
# lines(density(normal50.mul.matrix.p),col = c("#D1E5F0"))

boxplot.burn.normal = c(logit.mul.normal.matrix.p,
                normal500.mul.matrix.p,
                normal250.mul.matrix.p,
                normal100.mul.matrix.p,
                normal50.mul.matrix.p)
size = c(rep('1000',1000),rep('0500',1000),rep('0250',1000),rep('0100',1000),rep('0050',1000))
burn_size = data.frame(size = size, posterior = boxplot.burn.normal)
p.normal = ggplot(burn_size,aes(x = size,y = posterior))+
  xlab('Normal distribution')+
  geom_boxplot()

boxplot.burn.chris = c(logit.mul.chris.matrix.p,
                       chris500.mul.matrix.p,
                       chris250.mul.matrix.p,
                       chris100.mul.matrix.p,
                       chris50.mul.matrix.p)
size = c(rep('1000',1000),rep('0500',1000),rep('0250',1000),rep('0100',1000),rep('0050',1000))
burn_size = data.frame(size = size, posterior = boxplot.burn.chris)
p.chris = ggplot(burn_size,aes(x = size,y = posterior))+
  xlab('N(0,S^2/b)')+
  geom_boxplot()
p.chris

boxplot.burn.pi = c(logit.mul.pi.matrix.p,
                    pi500.mul.matrix.p,
                    pi250.mul.matrix.p,
                    pi100.mul.matrix.p,
                    pi50.mul.matrix.p)
size = c(rep('1000',1000),rep('0500',1000),rep('0250',1000),rep('0100',1000),rep('0050',1000))
burn_size = data.frame(size = size, posterior = boxplot.burn.pi)
p.pi = ggplot(burn_size,aes(x = size,y = posterior))+
  xlab('N(0,pi^2/3*(p+1))')+
  geom_boxplot()
p.pi



boxplot.burn.cauchy = c(logit.mul.cauchy.matrix.p,
                        logit.mul.cauchy500.matrix.p,
                        logit.mul.cauchy250.matrix.p,
                        logit.mul.cauchy100.matrix.p,
                        logit.mul.cauchy50.matrix.p)
size = c(rep('1000',1000),rep('0500',1000),rep('0250',1000),rep('0100',1000),rep('0050',1000))
burn_size = data.frame(size = size, posterior = boxplot.burn.cauchy)
p.cauchy = ggplot(burn_size,aes(x = size,y = posterior))+
  xlab('Cauchy')+
  geom_boxplot()
p.cauchy
arrange(p.normal,p.chris,p.pi,p.cauchy)

# Different methods in the same sample size
# n = 50
# beta0 n=50
# The model is death ~ age+tbsa+race+inh_inj
burn.model50=glm(formula = death[sample50] ~ age.scale[sample50]+tbsa.scale[sample50]+race[sample50]+inh_inj[sample50],family = binomial(link = 'logit'),data = input.burn)
summary(burn.model50)
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)            -474.3    94819.5  -0.005    0.996
# age.scale[sample50]     627.8    76028.2   0.008    0.993
# tbsa.scale[sample50]    837.4    99672.5   0.008    0.993
# race[sample50]         -214.1    82547.6  -0.003    0.998
# inh_inj[sample50]      -274.0    91384.9  -0.003    0.998
str(burn.model50)

burn.model50$coefficients[1]

# create the tag variable
x.tag = c('0-MLE',rep('1-N(0,10^6)',10000),rep('2-Cauchy',10000),rep('3-N(0,S^2/b)',10000),rep('4-N(0,pi^2/3*(p+1))',10000))


# n = 50
# beta0 n=50
burn.size50.beta0 = c(burn.model50$coefficients[1],
                      normal50.mul.matrix.c[,'b0'],
                      logit.mul.cauchy50.matrix.c[,'b0'],
                      chris50.mul.matrix.c[,'b0'],
                      pi50.mul.matrix.c[,'b0'])
data.burn.size50.beta0 = data.frame(x.tag = x.tag, posterior = burn.size50.beta0)

p.size50.beta0 = ggplot(data.burn.size50.beta0,aes(x.tag,posterior))+geom_boxplot()+
  xlab("The prior types")+ylab("beta0")+ggtitle("size = 50")

# beta1 n=50
burn.size50.beta1 = c(burn.model50$coefficients[1],
                      normal50.mul.matrix.c[,'b.age'],
                      logit.mul.cauchy50.matrix.c[,'b.age'],
                      chris50.mul.matrix.c[,'b.age'],
                      pi50.mul.matrix.c[,'b.age'])
data.burn.size50.beta1 = data.frame(x.tag = x.tag, posterior = burn.size50.beta1)
p.size50.beta1 = ggplot(data.burn.size50.beta1,aes(x.tag,posterior))+geom_boxplot()+
  xlab("The prior types")+ylab("b.age(scale)")+ggtitle("size = 50")

# beta2 n=50
burn.size50.beta2 = c(burn.model50$coefficients[2],
                      normal50.mul.matrix.c[,'b.tbsa'],
                      logit.mul.cauchy50.matrix.c[,'b.tbsa'],
                      chris50.mul.matrix.c[,'b.tbsa'],
                      pi50.mul.matrix.c[,'b.tbsa'])
data.burn.size50.beta2 = data.frame(x.tag = x.tag, posterior = burn.size50.beta2)
p.size50.beta2 = ggplot(data.burn.size50.beta2,aes(x.tag,posterior))+geom_boxplot()+
  xlab("The prior types")+ylab("b.tbsa(scale)")+ggtitle("size = 50")

# beta3 n=50
burn.size50.beta3 = c(burn.model50$coefficients[3],
                      normal50.mul.matrix.c[,'b.race'],
                      logit.mul.cauchy50.matrix.c[,'b.race'],
                      chris50.mul.matrix.c[,'b.race'],
                      pi50.mul.matrix.c[,'b.race'])
data.burn.size50.beta3 = data.frame(x.tag = x.tag, posterior = burn.size50.beta3)
p.size50.beta3 = ggplot(data.burn.size50.beta3,aes(x.tag,posterior))+geom_boxplot()+
  xlab("The prior types")+ylab("b.race")+ggtitle("size = 50")

# beta4 n=50
burn.size50.beta4 = c(burn.model50$coefficients[4],
                      normal50.mul.matrix.c[,'b.inh_inj'],
                      logit.mul.cauchy50.matrix.c[,'b.inh_inj'],
                      chris50.mul.matrix.c[,'b.inh_inj'],
                      pi50.mul.matrix.c[,'b.inh_inj'])
data.burn.size50.beta4 = data.frame(x.tag = x.tag, posterior = burn.size50.beta4)
p.size50.beta4 = ggplot(data.burn.size50.beta4,aes(x.tag,posterior))+geom_boxplot()+
  xlab("The prior types")+ylab("b.inh_inj")+ggtitle("size = 50")

arrange(p.size50.beta0,p.size50.beta1,p.size50.beta2,p.size50.beta3,p.size50.beta4)

# n = 100
burn.model100=glm(formula = death[sample100] ~ age.scale[sample100]+tbsa.scale[sample100]+race[sample100]+inh_inj[sample100],family = binomial(link = 'logit'),data = input.burn)

# beta0 n=100
burn.size100.beta0 = c(burn.model100$coefficients[1],
                       normal100.mul.matrix.c[,'b0'],
                       logit.mul.cauchy100.matrix.c[,'b0'],
                       chris100.mul.matrix.c[,'b0'],
                       pi100.mul.matrix.c[,'b0'])
data.burn.size100.beta0 = data.frame(x.tag = x.tag, posterior = burn.size100.beta0)

p.size100.beta0 = ggplot(data.burn.size100.beta0,aes(x.tag,posterior))+geom_boxplot()+
  xlab("The prior types")+ylab("b0")+ggtitle("size = 100")

# beta1 n=100
burn.size100.beta1 = c(burn.model100$coefficients[2],
                       normal100.mul.matrix.c[,'b.age'],
                       logit.mul.cauchy100.matrix.c[,'b.age'],
                       chris100.mul.matrix.c[,'b.age'],
                       pi100.mul.matrix.c[,'b.age'])
data.burn.size100.beta1 = data.frame(x.tag = x.tag, posterior = burn.size100.beta1)
p.size100.beta1 = ggplot(data.burn.size100.beta1,aes(x.tag,posterior))+geom_boxplot()+
  xlab("The prior types")+ylab("b.age")+ggtitle("size = 100")

# b.tbsa n=100
burn.size100.beta2 = c(burn.model100$coefficients[3],
                       normal100.mul.matrix.c[,'b.tbsa'],
                       logit.mul.cauchy100.matrix.c[,'b.tbsa'],
                       chris100.mul.matrix.c[,'b.tbsa'],
                       pi100.mul.matrix.c[,'b.tbsa'])
data.burn.size100.beta2 = data.frame(x.tag = x.tag, posterior = burn.size100.beta2)
p.size100.beta2 = ggplot(data.burn.size100.beta2,aes(x.tag,posterior))+geom_boxplot()+
  xlab("The prior types")+ylab("b.tbsa")+ggtitle("size = 100")

# beta3 n=100
burn.size100.beta3 = c(burn.model100$coefficients[4],
                       normal100.mul.matrix.c[,'b.race'],
                       logit.mul.cauchy100.matrix.c[,'b.race'],
                       chris100.mul.matrix.c[,'b.race'],
                       pi100.mul.matrix.c[,'b.race'])
data.burn.size100.beta3 = data.frame(x.tag = x.tag, posterior = burn.size100.beta3)
p.size100.beta3 = ggplot(data.burn.size100.beta3,aes(x.tag,posterior))+geom_boxplot()+
  xlab("The prior types")+ylab("b.race")+ggtitle("size = 100")

# beta4 n=100
burn.size100.beta4 = c(burn.model100$coefficients[5],
                       normal100.mul.matrix.c[,'b.inh_inj'],
                       logit.mul.cauchy100.matrix.c[,'b.inh_inj'],
                       chris100.mul.matrix.c[,'b.inh_inj'],
                       pi100.mul.matrix.c[,'b.inh_inj'])
data.burn.size100.beta4 = data.frame(x.tag = x.tag, posterior = burn.size100.beta4)
p.size100.beta4 = ggplot(data.burn.size100.beta4,aes(x.tag,posterior))+geom_boxplot()+
  xlab("The prior types")+ylab("b.inh_inj")+ggtitle("size = 100")

arrange(p.size100.beta0,p.size100.beta1,p.size100.beta2,p.size100.beta3,p.size100.beta4)

# n = 250
burn.model250=glm(formula = death[sample250] ~ age.scale[sample250]+tbsa.scale[sample250]+race[sample250]+inh_inj[sample250],family = binomial(link = 'logit'),data = input.burn)

# beta0 n=250
burn.size250.beta0 = c(burn.model250$coefficients[1],
                       normal250.mul.matrix.c[,'b0'],
                       logit.mul.cauchy250.matrix.c[,'b0'],
                       chris250.mul.matrix.c[,'b0'],
                       pi250.mul.matrix.c[,'b0'])
data.burn.size250.beta0 = data.frame(x.tag = x.tag, posterior = burn.size250.beta0)

p.size250.beta0 = ggplot(data.burn.size250.beta0,aes(x.tag,posterior))+geom_boxplot()+
  xlab("The prior types")+ylab("b0")+ggtitle("size = 250")

# beta1 n=250
burn.size250.beta1 = c(burn.model250$coefficients[2],
                       normal250.mul.matrix.c[,'b.age'],
                       logit.mul.cauchy250.matrix.c[,'b.age'],
                       chris250.mul.matrix.c[,'b.age'],
                       pi250.mul.matrix.c[,'b.age'])
data.burn.size250.beta1 = data.frame(x.tag = x.tag, posterior = burn.size250.beta1)
p.size250.beta1 = ggplot(data.burn.size250.beta1,aes(x.tag,posterior))+geom_boxplot()+
  xlab("The prior types")+ylab("b.age")+ggtitle("size = 250")

# beta2 n=250
burn.size250.beta2 = c(burn.model250$coefficients[3],
                       normal250.mul.matrix.c[,'b.tbsa'],
                       logit.mul.cauchy250.matrix.c[,'b.tbsa'],
                       chris250.mul.matrix.c[,'b.tbsa'],
                       pi250.mul.matrix.c[,'b.tbsa'])
data.burn.size250.beta2 = data.frame(x.tag = x.tag, posterior = burn.size250.beta2)
p.size250.beta2 = ggplot(data.burn.size250.beta2,aes(x.tag,posterior))+geom_boxplot()+
  xlab("The prior types")+ylab("b.tbsa")+ggtitle("size = 250")

# beta3 n=250
burn.size250.beta3 = c(burn.model250$coefficients[4],
                       normal250.mul.matrix.c[,'b.race'],
                       logit.mul.cauchy250.matrix.c[,'b.race'],
                       chris250.mul.matrix.c[,'b.race'],
                       pi250.mul.matrix.c[,'b.race'])
data.burn.size250.beta3 = data.frame(x.tag = x.tag, posterior = burn.size250.beta3)
p.size250.beta3 = ggplot(data.burn.size250.beta3,aes(x.tag,posterior))+geom_boxplot()+
  xlab("The prior types")+ylab("b.race")+ggtitle("size = 250")

# beta4 n=250
burn.size250.beta4 = c(burn.model250$coefficients[5],
                       normal250.mul.matrix.c[,'b.inh_inj'],
                       logit.mul.cauchy250.matrix.c[,'b.inh_inj'],
                       chris250.mul.matrix.c[,'b.inh_inj'],
                       pi250.mul.matrix.c[,'b.inh_inj'])
data.burn.size250.beta4 = data.frame(x.tag = x.tag, posterior = burn.size250.beta4)
p.size250.beta4 = ggplot(data.burn.size250.beta4,aes(x.tag,posterior))+geom_boxplot()+
  xlab("The prior types")+ylab("b.inh_inj")+ggtitle("size = 250")

arrange(p.size250.beta0,p.size250.beta1,p.size250.beta2,p.size250.beta3,p.size250.beta4)


# n = 500
burn.model500=glm(formula = death[sample500] ~ age.scale[sample500]+tbsa.scale[sample500]+race[sample500]+inh_inj[sample500],family = binomial(link = 'logit'),data = input.burn)
summary(burn.model500)
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)            -4.0505     0.4585  -8.835  < 2e-16 ***
#   age.scale[sample500]    1.9927     0.2755   7.233 4.71e-13 ***
#   tbsa.scale[sample500]   1.5470     0.2178   7.104 1.21e-12 ***
#   race[sample500]         0.5799     0.4180   1.387   0.1654    
# inh_inj[sample500]      1.1109     0.4580   2.425   0.0153 *  

# beta0 n=500
burn.size500.beta0 = c(burn.model500$coefficients[1],
                       normal500.mul.matrix.c[,'b0'],
                       logit.mul.cauchy500.matrix.c[,'b0'],
                       chris500.mul.matrix.c[,'b0'],
                       pi500.mul.matrix.c[,'b0'])
data.burn.size500.beta0 = data.frame(x.tag = x.tag, posterior = burn.size500.beta0)

p.size500.beta0 = ggplot(data.burn.size500.beta0,aes(x.tag,posterior))+geom_boxplot()+
  xlab("The prior types")+ylab("b0")+ggtitle("size = 500")

# beta1 n=500
burn.size500.beta1 = c(burn.model500$coefficients[2],
                       normal500.mul.matrix.c[,'b.age'],
                       logit.mul.cauchy500.matrix.c[,'b.age'],
                       chris500.mul.matrix.c[,'b.age'],
                       pi500.mul.matrix.c[,'b.age'])
data.burn.size500.beta1 = data.frame(x.tag = x.tag, posterior = burn.size500.beta1)
p.size500.beta1 = ggplot(data.burn.size500.beta1,aes(x.tag,posterior))+geom_boxplot()+
  xlab("The prior types")+ylab("b.age")+ggtitle("size = 500")

# beta2 n=500
burn.size500.beta2 = c(burn.model500$coefficients[3],
                       normal500.mul.matrix.c[,'b.tbsa'],
                       logit.mul.cauchy500.matrix.c[,'b.tbsa'],
                       chris500.mul.matrix.c[,'b.tbsa'],
                       pi500.mul.matrix.c[,'b.tbsa'])
data.burn.size500.beta2 = data.frame(x.tag = x.tag, posterior = burn.size500.beta2)
p.size500.beta2 = ggplot(data.burn.size500.beta2,aes(x.tag,posterior))+geom_boxplot()+
  xlab("The prior types")+ylab("b.tbsa")+ggtitle("size = 500")
9
# beta3 n=500
burn.size500.beta3 = c(burn.model500$coefficients[4],
                       normal500.mul.matrix.c[,'b.tbsa'],
                       logit.mul.cauchy500.matrix.c[,'b.race'],
                       chris500.mul.matrix.c[,'b.race'],
                       pi500.mul.matrix.c[,'b.race'])
data.burn.size500.beta3 = data.frame(x.tag = x.tag, posterior = burn.size500.beta3)
p.size500.beta3 = ggplot(data.burn.size500.beta3,aes(x.tag,posterior))+geom_boxplot()+
  xlab("The prior types")+ylab("b.race")+ggtitle("size = 500")

# beta4 n=500
burn.size500.beta4 = c(burn.model500$coefficients[5],
                       normal500.mul.matrix.c[,'b.inh_inj'],
                       logit.mul.cauchy500.matrix.c[,'b.inh_inj'],
                       chris500.mul.matrix.c[,'b.inh_inj'],
                       pi500.mul.matrix.c[,'b.inh_inj'])
data.burn.size500.beta4 = data.frame(x.tag = x.tag, posterior = burn.size500.beta4)
p.size500.beta4 = ggplot(data.burn.size500.beta4,aes(x.tag,posterior))+geom_boxplot()+
  xlab("The prior types")+ylab("b.inh_inj")+ggtitle("size = 500")

arrange(p.size500.beta0,p.size500.beta1,p.size500.beta2,p.size500.beta3,p.size500.beta4)

# Importance sampling: to estimate the mean and variance of the posterior

# Chapter 5 Exploration
# c-log-log
#define the cloglog function in R

expexp.pdf <- function(x) {exp(x)*exp(-exp(x))}

integrate(expexp.pdf, lower = -10, upper = 10)

beta.seq <- seq(-10,10,length=200)
plot(beta.seq,expexp.pdf(beta.seq),type="l",main="Cloglog pdf")

#--- Importance Sampling to calculate E[beta]
set.seed(1)
num.sims <- 100000
beta.star <- rnorm(num.sims,mean=0,sd=1)
MC.Estim.E.Beta <- mean(beta.star*cloglog.pdf(beta.star)/
                          dnorm(beta.star,mean=0,sd=1))
cat("Estimate of E[Beta]=",round(MC.Estim.E.Beta,3),"\n")

#--estimate the variance by estimating E[Beta^2]
MC.Estim.E.Beta.sq <- mean(beta.star^2*cloglog.pdf(beta.star)/
                             dnorm(beta.star,mean=0,sd=1))
MC.Estim.Var.Beta <- MC.Estim.E.Beta.sq-MC.Estim.E.Beta^2
cat("Estimate of Var[Beta]=",round(MC.Estim.Var.Beta,3),"\n")
sd.b <- sqrt(MC.Estim.Var.Beta)
# > sd.b
# [1] 1.252077

#--- Use Sampling Importance Resampling to generate a sample from cloglog dist'n
weights <- cloglog.pdf(beta.star)/dnorm(beta.star,mean=0,sd=10)
b.star.corrected <- sample(x=beta.star,size=num.sims,replace=TRUE,prob=weights)

#------ plots of theoretical pdf, and empirical histogram ---
par(mfrow=c(2,1),oma=c(0,0,3,0))
#Theoretical
plot(beta.seq,cloglog.pdf(beta.seq),type="l",xlab = "probability",main="Cloglog probability density function")
abline(v=MC.Estim.E.Beta,col="blue",lty = 2)
abline(v=c(-2*sd.b,2*sd.b),col="red",lty = 2)

# Empirical
hist(b.star.corrected,xlab="Simulated betas",ylab="")
abline(v=MC.Estim.E.Beta,col="blue",lty = 2)
abline(v=c(-2*sd.b,2*sd.b),col="red",lty = 2)

# Empirical density plot
plot(density(b.star.corrected),xlab="Simulated betas",ylab = "frequency",main="Beta distribution resampled according to importance sampling")
abline(v=MC.Estim.E.Beta,col="blue",lty = 2)
abline(v=c(-2*sd.b,2*sd.b),col="red",lty = 2)

# plug the random variable from into the c-log-log function
# b.star.corrected: the extreme value distribution, and these will be values on the 
# real number line, thus negative and positive values
# real number data???survival probabilities
random.original = 1-exp(-exp(b.star.corrected))

hist(random.original)
# get the uniform distribution

# when I plug those values into the 

# 5.2 frequentist analysis
# the same as the analysis before

# N(0,sigma^2)
par(mfrow=c(1,2))
normal.sigma = function(sigma){
  induced.prior.normal = matrix(data = NA,nrow = 100, ncol = 1000)
  induced.prior.random = matrix(data = NA,nrow = 100, ncol = 1000)
  set.seed(1)
  alpha = rnorm(1000,0,sigma)
  beta = rnorm(1000,0,sigma)
  x = rnorm(100,0,1)
  for(i in 1:100){
    induced.prior.normal[i,] = cloglog(alpha+beta*age.scale[i])
    induced.prior.random[i,] = cloglog(alpha+beta*x[i])
  }
  hist(induced.prior.normal,xlab = expression(pi(x)), ylab=("frequency"),main = paste("original data",'sigma=',sigma))
  hist(induced.prior.random,xlab = expression(pi(x)), ylab=("frequency"),main = paste("random normal variables",'sigma=',sigma))
}

par(mfrow = c(5,2))
normal.sigma(1.5)
normal.sigma(25)
normal.sigma(0.75)
normal.sigma(2.133)
normal.sigma(10^6)
normal.sigma(25)

#N(0,S^2/b)
# the age.scale has been standardized, so S = 1
normal.b = function(b){
  induced.prior.chris = matrix(data = NA,nrow = 100, ncol = 1000)
  induced.prior.random = matrix(data = NA,nrow = 100, ncol = 1000)
  set.seed(1)
  alpha = rnorm(1000,0,1^2/b)
  beta = rnorm(1000,0,1^2/b)
  x = rnorm(100,0,1)
  for(i in 1:100){
    induced.prior.chris[i,] = cloglog(alpha+beta*age.scale[i])
    induced.prior.random[i,] = cloglog(alpha+beta*x[i])
  }
  hist(induced.prior.chris,xlab = expression(pi(x)), ylab=("frequency"),main = paste("original data",'b =',b))
  hist(induced.prior.random,xlab = expression(pi(x)), ylab=("frequency"),main = paste("random normal variables",'b =',b))
}
par(mfrow = c(5,2))
normal.b(1)
normal.b(0.9)
normal.b(0.8)
normal.b(0.7)
normal.b(0.6)



#N(digamma(1),pi^2/(3*(p+1)))
par(mfrow = c(1,2))
p = 1
set.seed(1)
induced.prior.pi = matrix(data = NA,nrow = 100, ncol = 1000)
alpha = rnorm(1000,digamma(1),sqrt(pi^2/(3*(p+1))))
beta = rnorm(1000,0,sqrt(pi^2/(3*(p+1))))
x = rnorm(100,0,1)
for(i in 1:100){
  induced.prior.pi[i,] = expexp(alpha+beta*age.scale[i])
  induced.prior.random[i,] = expexp(alpha+beta*x[i])
}
hist(induced.prior.pi,xlab = expression(pi(x)), ylab=("frequency"),main = "pi prior distribution original dataset")
hist(induced.prior.random,xlab = expression(pi(x)), ylab=("frequency"),main = "pi prior distribution random variable N(0,1)")

# cauchy priors
cauchy.scale = function(scale){
  induced.prior.cauchy = matrix(data = NA,nrow = 100, ncol = 1000)
  induced.prior.random = matrix(data = NA,nrow = 100, ncol = 1000)
  set.seed(1)
  alpha = rcauchy(1000,0,scale)
  beta = rcauchy(1000,0,scale)
  x = rnorm(100,0,1)
  for(i in 1:100){
    induced.prior.cauchy[i,] = cloglog(alpha+beta*age.scale[i])
    induced.prior.random[i,] = cloglog(alpha+beta*x[i])
  }
  hist(induced.prior.cauchy,xlab = expression(pi(x)), ylab=("frequency"),main = paste("original data",'scale=',scale))
  hist(induced.prior.random,xlab = expression(pi(x)), ylab=("frequency"),main = ("random normal variables"))
}
par(mfrow = c(4,2))
scale = c(1.5,1.25,1,0.5,0.1)
cauchy.scale(0.5)


# CMPs
par(mfrow = c(4,2))
cmp.induced.clog = function(alpha1,beta1){
  pi40 = rbeta(1000,alpha1,beta1)
  pi60 = rbeta(1000,alpha1,beta1)
  x40 = (40-mean(age.scale))/sd(age.scale)
  x60 = (60-mean(age.scale))/sd(age.scale)
  m1 = exp(-exp(pi40))
  m2 = exp(-exp(pi60))
  beta = (m1-m2)/(x60-x40)
  alpha = 1-x40*beta-m1
  x = rnorm(100,0,1)
  induced.prior.cmp = matrix(data = NA,nrow = 100, ncol = 1000)
  induced.prior.random = matrix(data = NA,nrow = 100, ncol = 1000)
  
  
  for(i in 1:100){
    induced.prior.cmp[i,] = expexp(alpha+beta*age.scale[i])
    induced.prior.random[i,] = expexp(alpha+beta*x[i])
  }
  
  hist(induced.prior.cmp,xlab = expression(pi(x)), ylab=("frequency"),main = paste('beta(',alpha1,',',beta1,')',"CMPs method original dataset"))
  hist(induced.prior.random,xlab = expression(pi(x)), ylab=("frequency"),main = paste('beta(',alpha1,',',beta1,')',"CMPs method random variable"))
}

cmp.induced.clog(5,5)
cmp.induced.clog(1,1)
cmp.induced.clog(0.5,0.5)
cmp.induced.clog(0.25,0.25)
cmp.induced.clog(0.05,0.05)



par(mfrow = c(2,3))
# normal distribution
sigma = 25
induced.prior.normal = matrix(data = NA,nrow = 100, ncol = 1000)
set.seed(1)
alpha = rnorm(1000,0,sigma)
beta = rnorm(1000,0,sigma)
x = rnorm(100,0,1)
for(i in 1:100){
  induced.prior.normal[i,] = cloglog(alpha+beta*age.scale[i])
}
hist(induced.prior.normal,xlab = expression(pi(x)), ylab=("frequency"),main = paste("N(0,sigma^2)",'sigma=',sigma))

# christensen et al.
b = 0.8
induced.prior.chris = matrix(data = NA,nrow = 100, ncol = 1000)
set.seed(1)
alpha = rnorm(1000,0,1^2/b)
beta = rnorm(1000,0,1^2/b)
x = rnorm(100,0,1)
for(i in 1:100){
  induced.prior.chris[i,] = cloglog(alpha+beta*age.scale[i])
}
hist(induced.prior.chris,xlab = expression(pi(x)), ylab=("frequency"),main = paste("N(0,S^2/b)",'b =',b))

# pi
p = 1
set.seed(1)
induced.prior.pi = matrix(data = NA,nrow = 100, ncol = 1000)
alpha = rnorm(1000,digamma(1),sqrt(pi^2/(3*(p+1))))
beta = rnorm(1000,0,sqrt(pi^2/(3*(p+1))))
x = rnorm(100,0,0.5)
for(i in 1:100){
  induced.prior.pi[i,] = expexp(alpha+beta*age.scale[i])
}
hist(induced.prior.pi,xlab = expression(pi(x)), ylab=("frequency"),main = "pi^2 prior distribution Original dataset")


#Cauchy
scale= 1
induced.prior.cauchy = matrix(data = NA,nrow = 100, ncol = 1000)
set.seed(1)
alpha = rcauchy(1000,0,scale)
beta = rcauchy(1000,0,scale)
x = rnorm(100,0,1)
for(i in 1:100){
  induced.prior.cauchy[i,] = cloglog(alpha+beta*age.scale[i])
}
hist(induced.prior.cauchy,xlab = expression(pi(x)), ylab=("frequency"),main = paste("Cauchy distribution",'scale=',scale))

#CMPs
alpha1 = 1
  beta1 = 1
  pi40 = rbeta(1000,alpha1,beta1)
pi60 = rbeta(1000,alpha1,beta1)
x40 = (40-mean(age.scale))/sd(age.scale)
x60 = (60-mean(age.scale))/sd(age.scale)
m1 = exp(-exp(pi40))
m2 = exp(-exp(pi60))
beta = (m1-m2)/(x60-x40)
alpha = 1-x40*beta-m1
x = rnorm(100,0,1)
induced.prior.cmp = matrix(data = NA,nrow = 100, ncol = 1000)
induced.prior.random = matrix(data = NA,nrow = 100, ncol = 1000)


for(i in 1:100){
  induced.prior.cmp[i,] = expexp(alpha+beta*age.scale[i])
  induced.prior.random[i,] = expexp(alpha+beta*x[i])
}

hist(induced.prior.cmp,xlab = expression(pi(x)), ylab=("frequency"),main = paste('beta(',alpha1,',',beta1,')',"CMPs method original dataset"))



# 5.3 coefficients
model.cloglog = glm(y ~ x,data = input.chdage,family = binomial(link = cloglog))
# Coefficients:
#              Estimate Std. Error z value Pr(>|z|)    
# (Intercept)    -0.7350     0.1864  -3.943 8.05e-05 ***
#   x             0.9281     0.1910   4.858 1.18e-06 ***

summary(model.cloglog)
c_alpha2 = c( -0.7350 -1.96*0.1864/sqrt(100),-0.7350+1.96*0.1864/sqrt(100))
c_beta2 = c(0.9281-1.96*0.1910/sqrt(100), 0.9281+1.96*0.1910/sqrt(100))
# > c_alpha2
# [1] -0.7715344 -0.6984656
# > c_beta2
# [1] 0.890664 0.965536



# Naive normal
#-- Normal Prior Case ---------------------
beta.sd <- 25
beta.prec <- 1/beta.sd^2

#--- Data file ---
cloglog.normal.data <- list(chd=chd,age.scale=age.scale,beta.prec=beta.prec)

#--- Initial Values for Parameters
set.seed(1)
cloglog.normal.inits <- list()
for(i in 1:num.chains) {
  b0 <- rnorm(n=1,mean=0,sd=1)
  b1 <- rnorm(n=1,mean=0,sd=1)
  cloglog.normal.inits[[i]] <- list(b0=b0,b1=b1)
}


cloglog_naivenormal <- "model {
 # priors
 b0 ~ dnorm(0,beta.prec)
 b1 ~ dnorm(0,beta.prec)

 # obs'n model
 for(i in 1:100) {
  cloglog(p[i]) <- b0+b1*age.scale[i]
  chd[i] ~ dbern(p[i])
 }
}
"
cat(cloglog_naivenormal,file="D:/Study/Postgraduate/paper/cloglog_naivenormal.txt")

#-- run jags
# Initialize the model
cloglog.init.JAGS <- jags.model( file = "D:/Study/Postgraduate/paper/cloglog_naivenormal.txt",
                              data=cloglog.normal.data, 
                              inits=cloglog.normal.inits, 
                              n.chains=num.chains,
                              n.adapt=adapt.length)
# specify parameters to monitor
wide.params <- c("b0","b1",'p')

# Run JAGS, using coda.samples a wrapper for jags.samples
cloglog.normal.samples <- coda.samples(model=cloglog.init.JAGS,
                             variable.names=wide.params,
                             n.iter=burnin+chain.length)
#create the matrix

cloglog.normal.matrix = rbind(cloglog.normal.samples[[1]],cloglog.normal.samples[[2]],cloglog.normal.samples[[3]])

#plot(cloglog.normal.samples)
cloglog.jags.normal.results <- summary(window(cloglog.normal.samples,start=(burnin+1)))
print(cloglog.jags.normal.results$statistics)
#            Mean         SD     Naive SE Time-series SE
# b0     -0.76314399 0.18877696 0.0014070604   0.0023149575
# b1      0.94329413 0.18741785 0.0013969302   0.0022753903

c(cloglog.jags.normal.results$statistics[1,1]-1.96*cloglog.jags.normal.results$statistics[1,2]/sqrt(100),cloglog.jags.normal.results$statistics[1,1]+1.96*cloglog.jags.normal.results$statistics[1,2]/sqrt(100))
c(cloglog.jags.normal.results$statistics[2,1]-1.96*cloglog.jags.normal.results$statistics[2,2]/sqrt(100),cloglog.jags.normal.results$statistics[2,1]+1.96*cloglog.jags.normal.results$statistics[2,2]/sqrt(100))

# [1] -0.8001443 -0.7261437
# [1] 0.9065602 0.9800280



# Christensen
b      <- 0.8   # much like pi prior
beta.sd <- 1^2/b
beta.prec <- 1/beta.sd

#--- Data file ---
cloglog.Chris.data <- list(n=n,chd=chd,age.scale=age.scale,
                   beta.prec=beta.prec)

#--- Initial Values for Parameters
set.seed(1)
cloglog.Chris.inits <- list()
for(i in 1:num.chains) {
  b0 <- rnorm(n=1,mean=0,sd=beta.sd)
  b1 <- rnorm(n=1,mean=0,sd=beta.sd)
  cloglog.Chris.inits[[i]] <- list(b0=b0,b1=b1)
}


cloglog.Chris.model <- "model {
 # priors
 b0 ~ dnorm(0,beta.prec)
 b1 ~ dnorm(0,beta.prec)

 # obs'n model
 for(i in 1:n) {
  cloglog(p[i]) <- b0+b1*age.scale[i]
  chd[i] ~ dbern(p[i])
 }
}
"
cat(cloglog.Chris.model,file="D:/Study/Postgraduate/paper/cloglog.Chris_model.txt")

#-- run jags
# Initialize the model
cloglog.Chris.init.JAGS <- jags.model(file = "D:/Study/Postgraduate/paper/cloglog.Chris_model.txt",
                              data=cloglog.Chris.data, inits=cloglog.Chris.inits, n.chains=num.chains,
                              n.adapt=adapt.length)

# specify parameters to monitor
Chris.params <- c("b0","b1",'p')

# Run JAGS, using coda.samples a wrapper for jags.samples
cloglog.Chris.samples <- coda.samples(model=cloglog.Chris.init.JAGS,
                              variable.names=Chris.params,
                              n.iter=burnin+chain.length)
cloglog.Chris.matrix = rbind(cloglog.Chris.samples[[1]],cloglog.Chris.samples[[2]],cloglog.Chris.samples[[3]])

plot(cloglog.Chris.samples)
cloglog.Chris.results <- summary(window(cloglog.Chris.samples,start=(burnin+1)))
print(cloglog.Chris.results$statistics)
#             Mean         SD     Naive SE Time-series SE
# b0     -0.73321671 0.18257117 0.0013608051   0.0021875164
# b1      0.90923025 0.18327827 0.0013660755   0.0021792315
c(cloglog.Chris.results$statistics[1,1]-cloglog.Chris.results$statistics[1,2]/sqrt(100),cloglog.Chris.results$statistics[1,1]+cloglog.Chris.results$statistics[1,2]/sqrt(100))
c(cloglog.Chris.results$statistics[2,1]-cloglog.Chris.results$statistics[2,2]/sqrt(100),cloglog.Chris.results$statistics[2,1]+cloglog.Chris.results$statistics[2,2]/sqrt(100))
# [1] -0.7514738 -0.7149596
# [1] 0.8909024 0.9275581

#create the matrix
cloglog.chris.matrix=rbind(cloglog.Chris.samples[[1]],cloglog.Chris.samples[[2]],cloglog.Chris.samples[[3]])

# Cauchy(0,1)
#--- Data file ---
cloglog.Cauchy.data <- list(chd=chd,age.scale=age.scale)

#--- Initial Values for Parameters
set.seed(1)
cloglog.Cauchy.inits <- list()
for(i in 1:num.chains) {
  b0 <- rnorm(n=1,mean=0,sd=1)
  b1 <- rnorm(n=1,mean=0,sd=1)
  cloglog.Cauchy.inits[[i]] <- list(b0=b0,b1=b1)
}


cloglog.Cauchy.model <- "
model{
  b0~dt(0,1,1)
  b1~dt(0,1,1)
  
  for (i in 1:100) {
      chd[i] ~ dbern(p[i])
      cloglog(p[i]) = b0+b1*age.scale[i]
  }
}"
cat(cloglog.Cauchy.model,file="D:/Study/Postgraduate/paper/cloglog.Cauchy_model.txt")

#-- run jags
# Initialize the model
cloglog.Cauchy.init.JAGS <- jags.model(file = "D:/Study/Postgraduate/paper/cloglog.Cauchy_model.txt",
                                      data=cloglog.Cauchy.data, inits=cloglog.Cauchy.inits, n.chains=num.chains,
                                      n.adapt=adapt.length)

# specify parameters to monitor
Cauchy.params <- c("b0","b1",'p')

# Run JAGS, using coda.samples a wrapper for jags.samples
cloglog.Cauchy.samples <- coda.samples(model=cloglog.Cauchy.init.JAGS,
                                      variable.names=Cauchy.params,
                                      n.iter=burnin+chain.length)
cloglog.Cauchy.matrix=rbind(cloglog.Cauchy.samples[[1]],cloglog.Cauchy.samples[[2]],cloglog.Cauchy.samples[[3]])


plot(cloglog.Cauchy.samples)
cloglog.Cauchy.results <- summary(window(cloglog.Cauchy.samples,start=(burnin+1)))
print(cloglog.Cauchy.results$statistics)
#            Mean         SD     Naive SE Time-series SE
# b0     -0.71247863 0.18230079 0.0013587899   0.0021930673
# b1      0.88980370 0.18529870 0.0013811350   0.0022081104
c(-0.712-0.182/sqrt(100),-0.712+0.182/sqrt(100))
c(0.890-0.185/sqrt(100),0.890+0.185/sqrt(100))
# > c(-0.712-0.182/sqrt(100),-0.712+0.182/sqrt(100))
# [1] -0.7302 -0.6938
# > c(0.890-0.185/sqrt(100),0.890+0.185/sqrt(100))
# [1] 0.8715 0.9085


#create the matrix
#create the matrix


# 5.3.4 Use the normal(0,pi^2/(3*(p+1))) prior
#------------------- pi Prior Case ---------------------
p       <- 1
beta.sd <- pi^2/sqrt(6*(p+1))
beta.prec <- 1/beta.sd
n = 100
#--- Data file ---
cloglog.pi.data <- list(n=n,chd=chd,age.scale=age.scale,beta.prec=beta.prec,mu = -0.579)

#--- Initial Values for Parameters
set.seed(1)
cloglog.pi.inits <- list()
for(i in 1:num.chains) {
  b0 <- rnorm(n=1,mean=digamma(1),sd=beta.sd)
  b1 <- rnorm(n=1,mean=0,sd=beta.sd)
  cloglog.pi.inits[[i]] <- list(b0=b0,b1=b1)
}

cloglog.pi.model <- "model {
 # priors
 b0 ~ dnorm(mu,beta.prec)
 b1 ~ dnorm(0,beta.prec)

 # obs'n model
 for(i in 1:n) {
  cloglog(p[i]) <- b0+b1*age.scale[i]
  chd[i] ~ dbern(p[i])
 }
}
"
cat(cloglog.pi.model,file="D:/Study/Postgraduate/paper/cloglog_pi_model.txt")

#-- run jags
# Initialize the model
cloglog.pi.init.JAGS <- jags.model(file = "D:/Study/Postgraduate/paper/cloglog_pi_model.txt",
                            data=cloglog.pi.data, inits=cloglog.pi.inits, n.chains=num.chains,
                            n.adapt=adapt.length)

# specify parameters to monitor
cloglog.pi.params <- c("b0","b1",'p')

# Run JAGS, using coda.samples a wrapper for jags.samples
cloglog.pi.samples <- coda.samples(model=cloglog.pi.init.JAGS,
                           variable.names=cloglog.pi.params,
                           n.iter=burnin+chain.length)
cloglog.pi.matrix=rbind(cloglog.pi.samples[[1]],cloglog.pi.samples[[2]],cloglog.pi.samples[[3]])

#create the matrix
b0 = c(pi.samples[[1]][,1],pi.samples[[2]][,1],pi.samples[[3]][,1])
b1 = c(pi.samples[[1]][,2],pi.samples[[2]][,2],pi.samples[[3]][,2])
cloglog.pi.matrix = matrix(NA,nrow = 18000,ncol = 2)
cloglog.pi.matrix[,1]=b0
cloglog.pi.matrix[,2]=b1

plot(cloglog.pi.samples)
cloglog.pi.results <- summary(window(cloglog.pi.samples,start=(burnin+1)))
print(cloglog.pi.results$statistics)
#         Mean         SD     Naive SE Time-series SE
# b0     -0.75049408 0.18706727 0.0013943171   0.0022628098
# b1      0.92215922 0.18449505 0.0013751449   0.0022011768
c(-0.750-1.96*0.187/sqrt(100),-0.750+1.96*0.187/sqrt(100))
c(0.922-1.96*0.184/sqrt(100),0.922+1.96*0.184/sqrt(100))
# > c(-0.750-1.96*0.187/sqrt(100),-0.750+1.96*0.187/sqrt(100))
# [1] -0.786652 -0.713348
# > c(0.922-1.96*0.184/sqrt(100),0.922+1.96*0.184/sqrt(100))
# [1] 0.885936 0.958064


#CMPs
x40 = (40-mean(aplore3::chdage$age))/sd(aplore3::chdage$age)
x60 = (60-mean(aplore3::chdage$age))/sd(aplore3::chdage$age)
modelString_cloglogcmp="
model{
  pi40 ~ dbeta(1,1)
  pi60 ~ dbeta(1,1)
  m1 = log(-log(pi40))
  m2 = log(-log(pi60))
  beta = (m2-m1)/(x60-x40)
  alpha = m2-beta*x60
  for (i in 1:100){
      y[i] ~ dbern(p[i])
      cloglog(p[i]) = alpha+beta*x[i]
  }
}"

writeLines(modelString_cloglogcmp, con='D:/Study/Postgraduate/paper/modelString_cloglogcmp.bug')
# create the input dataset
input.chdage.cloglogcmp = list('x' = age.scale,'y' = chd,'x40' = x40,'x60' = x60) 

chd.jags.cloglogcmp <- jags(data=input.chdage.cloglogcmp,model.file='D:/Study/Postgraduate/paper/modelString_cloglogcmp.bug',
                     param=c('alpha','beta','p'),
                     n.chains=3, n.iter=19000,n.burnin = 1000, n.thin=3)
print(chd.jags.cloglogcmp)
#3 chains, each with 19000 iterations (first 1000 discarded), n.thin = 3
# n.sims = 18000 iterations saved
#           mu.vect sd.vect    2.5%     25%     50%     75%   97.5%  Rhat n.eff
# alpha     -0.744   0.185  -1.126  -0.863  -0.737  -0.616  -0.401 1.001 18000
# beta       0.904   0.185   0.556   0.778   0.900   1.029   1.275 1.001 18000

c(-0.744-0.185/sqrt(100),-0.744+0.185/sqrt(100))
c(0.904-0.185/sqrt(100),0.904+0.185/sqrt(100))
# > c(-0.744-0.185/sqrt(100),-0.744+0.185/sqrt(100))
# [1] -0.7625 -0.7255
# > c(0.904-0.185/sqrt(100),0.904+0.185/sqrt(100))
# [1] 0.8855 0.9225

# Create the matrix
cloglog.cmp.matrix = chd.jags.cloglogcmp$BUGSoutput$sims.matrix[,4:103]

# 5.4 posterior distribution
# draw the prediction graph of the mle method
library(reshape2)

summary(model.cloglog)
par(mfrow = c(1,1))
fitted.prob.prevalence = predict(chd.glm,type = "response")
plot(sort(age.scale),
     expexp(-0.763+0.944*age.scale),type = 'l',ylab = "")

age.order = (aplore3::chdage$age)

#normal
fitted.bayes1 = expexp(-0.762+0.944*age.scale)
#cauchy
fitted.bayes2 = expexp(-0.712+0.890*age.scale)
#Chris
fitted.bayes3 = expexp(-0.733+0.909*age.scale)
#pi
fitted.bayes4 = expexp(-0.750+0.922*age.scale)
#cmps
fitted.bayes5 = expexp(-0.746+0.906*age.scale)

# N(0,10^6)
lines(sort(age.scale),fitted.bayes1[order(age.scale)],lty = 2, col = 'red')
abline(h = expexp(-0.762+0.944*x60),col = 'red')
# Cauchy
lines(sort(age.scale),fitted.bayes2[order(age.scale)],lty = 2, col = 'blue')
abline(h = expexp(-0.712+0.890*x60),col = 'blue',lty = 2)
# chris
lines(sort(age.scale),fitted.bayes3[order(age.scale)],lty = 2, col = 'orange')
abline(h = expexp(-0.733+0.909*x60),col = 'orange',lty = 2)
# N(0,pi^2/3)
lines(sort(age.scale),fitted.bayes4[order(age.scale)],lty = 2, col = 'pink')
abline(h = expexp(-0.750+0.922*x60),col = 'pink',lty = 2)
# CMPs
lines(sort(age.scale),fitted.bayes5[order(age.scale)],lty = 2, col = 'purple')
abline(h = expexp(-0.746+0.906*x60),col = 'purple',lty = 2)
abline(v = x60,lty = 2)

axis(1,at=20:69,las = 3)
legend("topleft",legend = c("MLE","N(0,10^6)","Cauchy","N(0,S^2/b","N(0,pi^2/3)","CMPs"),
       col = c('black','red','blue','green','pink','purple'),lty = c(1,2,2,2,2,2),cex = 1)

par(mfrow=c(3,2))
plot(x = fitted.prob.prevalence,y = fitted.bayes1,xlab = "CHD rate(MLE)",ylab = "CHD rate",pch = 20,main = expression(N(0,1.2^2)))
abline(a = 0 ,b = 1, col  = "blue")
plot(x = fitted.prob.living,y = fitted.bayes2,xlab = "CHD rate(MLE)",ylab = "CHD rate",pch = 20,main = "Cauchy(0,0.5)" )
abline(a = 0 ,b = 1, col  = "blue")
plot(x = fitted.prob.living,y = fitted.bayes3,xlab = "CHD rate(MLE)",ylab = "CHD rate",pch = 20,main = "N(0,S^2/b)" )
abline(a = 0 ,b = 1, col  = "blue")
plot(x = fitted.prob.living,y = fitted.bayes4,xlab = "CHD rate(MLE)",ylab = "CHD rate",pch = 20,main = "N(0,pi^2/3)" )
abline(a = 0 ,b = 1, col  = "blue")
plot(x = fitted.prob.living,y = fitted.bayes5,xlab = "CHD rate(MLE)",ylab = "CHD rate",pch = 20,main = "CMPs" )
abline(a = 0 ,b = 1, col  = "blue")

# show the fluctuation of b

data.mle.x = c()
data.normal.y = c()
data.cauchy.y = c()
data.chris.y = c()
data.pi.y = c()
data.cmps.y = c()

data.mle.x = rep(fitted.prob.prevalence,3000)

posterior.survival.matrix.normal <- cloglog.normal.matrix[,3:102]
data.normal.y <- as.vector(t(posterior.survival.matrix.normal))

posterior.survival.matrix.cauchy<- cloglog.cauchy.matrix[,3:102]
data.cauchy.y <- as.vector(t(posterior.survival.matrix.cauchy))

posterior.survival.matrix.chris<- cloglog.chris.matrix[,3:102]
data.chris.y <- as.vector(t(posterior.survival.matrix.chris))

posterior.survival.matrix.pi<- cloglog.pi.matrix[,3:102]
data.pi.y <- as.vector(t(posterior.survival.matrix.pi))

Bayes.output.cmps <-chd.jags.cmp$BUGSoutput$sims.matrix
posterior.survival.matrix.cmps <- Bayes.output.cmps[,4:103]
data.cmps.y <- as.vector(t(posterior.survival.matrix.cmps))


# show the deviation of the Bayesian analysis

subset.values <- sort(sample(1:1800000,size=10000,replace=FALSE))
data.mle.x.subset = data.mle.x[subset.values]
data.normal.y.subset = data.normal.y[subset.values]
data.cauchy.y.subset = data.cauchy.y[subset.values]
data.chris.y.subset = data.chris.y[subset.values]
data.pi.y.subset = data.pi.y[subset.values]
data.cmps.y.subset = data.cmps.y[subset.values]

par(mfrow = c(2,3))
plot(data.mle.x.subset,data.normal.y.subset,pch = 20,xlab = "Estimate CHD rate (Frequencist)",ylab = "Estimate CHD rate (Bayesian)",main = "N(0,10^6)")
abline(a = 0 ,b = 1, col  = "blue")
plot(data.mle.x.subset,data.cauchy.y.subset,pch = 20,xlab = "Estimate CHD rate (Frequencist)",ylab = "Estimate CHD rate (Bayesian)",main = "Cauchy(0,0.8)")
abline(a = 0 ,b = 1, col  = "blue")
plot(data.mle.x.subset,data.chris.y.subset,pch = 20,xlab = "Estimate CHD rate (Frequencist)",ylab = "Estimate CHD rate (Bayesian)",main = "N(0,S^2/b),b = 0.8")
abline(a = 0 ,b = 1, col  = "blue")
plot(data.mle.x.subset,data.pi.y.subset,pch = 20,xlab = "Estimate CHD rate (Frequencist)",ylab = "Estimate CHD rate (Bayesian)",main = "pi priors")
abline(a = 0 ,b = 1, col  = "blue")
plot(data.mle.x.subset,data.cmps.y.subset,pch = 20,xlab = "Estimate CHD rate (Frequencist)",ylab = "Estimate CHD rate (Bayesian)",main = "CMPs methods")
abline(a = 0 ,b = 1, col  = "blue")

