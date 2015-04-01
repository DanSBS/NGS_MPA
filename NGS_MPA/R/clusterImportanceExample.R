# Testing out some concepts
# 
# Daniel V. Samarov
###############################################################################

require(MASS)
c1 <- mvrnorm(200, c(-1,0), matrix(c(.1,.25,.25,.9),2,2))
c2 <- mvrnorm(200, c(1, 0), matrix(c(.1,.25,.25,.9),2,2))
c3 <- mvrnorm(200, c(-1,-3), matrix(c(.1,.25,.25,.9),2,2))
Xs <- rbind(c1,c2,c3)

