function e=g(x)

e=x+exp(-x.^2/2)./( normcdf(x).*sqrt(2*pi));