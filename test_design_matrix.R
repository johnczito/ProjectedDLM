source("helpers/_helpers.R")

y = c(1, 4, 3, 5, 2, 6, 4)
TT = length(y)
p = 3
intercept = TRUE

X_manual = cbind(rep(1, TT - p),
                 c(3, 5, 2, 6),
                 c(4, 3, 5, 2),
                 c(1, 4, 3, 5))

X_helper_1 = arp_design_matrix(y, p, intercept)

X_helper_2 = varp_design_matrix(matrix(y), p, intercept)

X_manual == X_helper_1
X_manual == X_helper_2

arp_design_matrix(y, 1, FALSE) == matrix(y[1:(TT - 1)])

n = 2
Y = matrix(rnorm(n * TT), TT, n)

varp_design_matrix(Y, 1, FALSE) == Y[1:(TT - 1), ]