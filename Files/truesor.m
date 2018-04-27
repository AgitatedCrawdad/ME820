% x = fminsearch(@fsor,1);
% y = fminsearch(@fslor,1);
% z = fminsearch(@fadi,0.1);
% z = fminsearch(@fadibig,1.33)
z = fminbnd(@fadibig,1.2,1.33);
% z = fminbnd(@fslorbig,1.3,1.35);
% z = fminbnd(@fslorbig,1.97,1.999);