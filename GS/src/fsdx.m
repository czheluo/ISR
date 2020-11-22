function fx=fsdx(x)
fx=.4488*exp(103.62*x-13.298*x.^2)./(1+.3659*exp(106.66*x))-.02;