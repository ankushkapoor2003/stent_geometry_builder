function [xx] = Recreate_xx(xx_ind, xx_fixed)
xx = [xx_ind(:,1:2) xx_fixed(:,1:2) xx_ind(:,3:18) xx_fixed(:,3:5)];
return