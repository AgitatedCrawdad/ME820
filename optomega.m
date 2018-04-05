x = 0:0.25:6;
y = 0:0.25:4;

p = length(x);
q = length(y);
% p = 45;
% q=45;

t = cos(pi/p)+cos(pi/q);

omega = roots([t^2 -16 16]);