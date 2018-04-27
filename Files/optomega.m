x = 0:0.025:6;
y = 0:0.025:4;

p = length(x);
q = length(y);
% p = 45;
% q=45;

t = cos(pi/p)+cos(pi/q);

omega = roots([t^2 -16 16]);