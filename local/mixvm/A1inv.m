function [A1] = A1inv(R)
    
if R < 0.53
  A1 = 2*R + R^3 + 5*R^5/6;
elseif R>=0.53 && R<0.85
  A1 = -.4 + 1.39*R + 0.43/(1-R);
else
  A1 = 1/(R^3 - 4*R^2 + 3*R);
end
   
end