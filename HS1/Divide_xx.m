function [xx_ind, xx_fixed] = Divide_xx(xx)
% Separating into independent Variables to be be optimised and fixed
% variables to be kept at baseline value. Only repair variables and x(4) kept as fixed here.  
xx_ind_index = [1 2 5:20];
xx_fixed_index = [3 4 21:23];
xx_ind = xx(:,xx_ind_index);
xx_fixed = xx(:,xx_fixed_index);
return