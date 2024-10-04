function [xx_ind, xx_fixed] = Divide_xx(xx)
% Separating into independent Variables to be be optimised and fixed
% variables to be kept at baseline value. Only repair variables and x(4) kept as fixed here.  
xx_ind_index = [1:7 9:18];
xx_fixed_index = [8];
xx_ind = xx(:,xx_ind_index);
xx_fixed = xx(:,xx_fixed_index);
return