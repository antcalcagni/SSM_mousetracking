function [c] = expscale_besseli(x)

c = exp(-x)*besseli(0,x);

end