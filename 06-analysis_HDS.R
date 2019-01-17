library(INLA)
library(inlabru)

# Practice ---------------

data(Seeds)

df = data.frame(y = Seeds$r, Ntrials = Seeds$n, Seeds[, 3:5])

family1 = "binomial"
control.family1 = list(control.link=list(model="logit"))
# number of trials is df$Ntrials

hyper1 = list(theta = list(prior="pc.prec", param=c(1,0.01)))
formula1 = y ~ x1 + x2 + f(plate, model="iid", hyper=hyper1)

q1 = INLA:::inla.rw1(n = 5)

crossprod(q1)

INLA:::inla.rw2(n = 5)

# Mesh example

s <- 3 ### this factor will only changes C, not G
pts <- rbind(c(1,1), c(2,1),
             c(2.6, 1), c(0.7,1.7), 4:5/3, c(2,1.7))*s
n <- nrow(pts)
mesh0 <- inla.mesh.2d(pts[1:2,], max.edge=3*s,
                      offset=1*s, n=6, cutoff=s*1/2)
mesh <- inla.mesh.2d(rbind(c(3.3,1)*s, c(2.4,2)*s,
                           mesh0$loc[-c(3:4),1:2]),
                     max.edge=3*s, offset=1e-5, cutoff=s*1/2, n=100)
(m <- mesh$n)
dmesh <- inla.mesh.dual(mesh)
fem <- inla.mesh.fem(mesh, order=1)
A <- inla.spde.make.A(mesh, pts)

# Geostats example

data("SPDEtoy")

str(SPDEtoy)

pl.dom <- cbind(c(0, 1, 1, 0.7, 0), c(0, 0, 0.7, 1, 1))

mesh5 <- inla.mesh.2d(pl.dom, max.e = c(0.092, 0.2))


spde5 <- inla.spde2.pcmatern(
  mesh = mesh5,
  alpha = 2,
  ### mesh and smoothness parameter
  prior.range = c(0.3, 0.5),
  ### P(practic.range<0.3)=0.5
  prior.sigma = c(1, 0.01)
) ### P(sigma>1)=0.01

coords <- as.matrix(SPDEtoy[,1:2])
A5 <- inla.spde.make.A(mesh5, loc=coords)

stk5 <- inla.stack(
  data = list(resp = SPDEtoy$y),
  A = list(A5, 1),
  effects = list(i = 1:spde5$n.spde,
                 m = rep(1, nrow(SPDEtoy))),
  tag = 'est'
)

res5 <- inla(
  resp ~ 0 + m + f(i, model = spde5),
  data = inla.stack.data(stk5),
  control.predictor = list(A = inla.stack.A(stk5))
)

res5.field <- inla.spde2.result(res5, 'i', spde5, do.transf=TRUE)


# Prediction during estimation
pts3 <- rbind(c(.1,.1), c(.5,.55), c(.7,.9))

dim(A5pts3 <- inla.spde.make.A(mesh5, loc=pts3))

(jj3 <- which(colSums(A5pts3)>0))

round(A5pts3[, jj3],3)

stk5p.rf <- inla.stack(
  data = list(resp = NA),
  A = list(A5pts3),
  effects = list(i = 1:spde5$n.spde),
  tag = 'prd5r'
)

stk5.jp <- inla.stack(stk5, stk5p.rf)

res5p <- inla(
  resp ~ 0 + m + f(i, model = spde5),
  data = inla.stack.data(stk5.jp),
  control.predictor = list(A = inla.stack.A(stk5.jp), compute = TRUE)
)

# Prediction after estimation

drop(A5pts3%*%res5$summary.random$i$mean)

inla.mesh.project(inla.mesh.projector(mesh5, loc = pts3),
                  res5$summary.random$i$mean)

# Project onto mesh

# project on grid with domain [0,1] x [0,1]
pgrid0 <- inla.mesh.projector(mesh5, xlim=0:1, ylim=0:1, dims=c(101,101))

# project posterior mean and standard deviation
prd0.m <- inla.mesh.project(pgrid0, res5$summary.ran$i$mean)
prd0.s <- inla.mesh.project(pgrid0, res5$summary.ran$i$s)


# Predict response
stk5.presp <- inla.stack(
  data = list(resp = NA),
  A = list(A5pts3, 1),
  effects = list(i = 1:spde5$n.spde, m = rep(1, 3)),
  tag = 'prd5.resp'
)

stk5.full <- inla.stack(stk5, stk5.presp)

r5presp <- inla(
  resp ~ 0 + m + f(i, model = spde5),
  data = inla.stack.data(stk5.full),
  control.predictor = list(A = inla.stack.A(stk5.full), compute = TRUE)
)

(indd3r <- inla.stack.index(stk5.full, 'prd5.resp')$data)

round(r5presp$summary.fitted.values[indd3r,], 3)

# Marginal distribution of point response
marg3r <- r5presp$marginals.fitted.values[indd3r]

# Response on a grid, disabling marginal distributions.
stkgrid <- inla.stack(
  data = list(resp = NA),
  A = list(pgrid0$proj$A, 1),
  effects = list(i = 1:spde5$n.spde,
                 m = rep(1, 101 * 101)),
  tag = 'prd.gr'
)

stk.all <- inla.stack(stk5, stkgrid)
res5g <- inla(
  resp ~ 0 + m + f(i, model = spde5),
  data = inla.stack.data(stk.all),
  control.predictor = list(A = inla.stack.A(stk.all),
                           compute = TRUE),
  quantiles = NULL,
  control.results = list(
    return.marginals.random = FALSE,
    return.marginals.predictor = FALSE
  )
)
res5g$cpu

# Indexes of prediction grid
igr <- inla.stack.index(stk.all, 'prd.gr')$data


