function [circmean] = circ_mean(x)

sinr = sum(sin(x));
cosr = sum(cos(x));
circmean = atan2(sinr, cosr);

end