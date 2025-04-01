function [xx_ind, xx_fixed] = Divide_xx(xx)
% Separating into independent Variables to be be optimised and fixed
% variables to be kept at baseline value. Only repair variables kept as fixed here.  
xx_ind_index = [1:5];
xx_fixed_index = [6:8];
xx_ind = xx(:,xx_ind_index);
xx_fixed = xx(:,xx_fixed_index);
return