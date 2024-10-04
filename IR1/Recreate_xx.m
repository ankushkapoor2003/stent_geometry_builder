function [xx] = Recreate_xx(xx_ind, xx_fixed)
xx = [xx_ind(:,1:7) xx_fixed(:,1) xx_ind(:,8:17)];
return