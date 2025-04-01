function [xrep_new] = Repair(x_ind, x_rep)

    % Evaluating bounds on repair and dependent variables
    xrep_ll = 0.5*x_rep;
    xrep_ul = 1.5*x_rep;
    cr_yscale_limit = [0.5 2];
    cr_yscales_limit = [0.5 2.5]; 
 

    % The objective function for repair aims to enforce constraints and minimise the difference between connector height and height of AB Strut 
    % Objective function,constraint evaluation, and fmincon
    fun = @(x_rep)obj(x_rep, x_ind);
    mycon = @(x_rep)cons(x_rep, x_ind, cr_yscale_limit, cr_yscales_limit);
    x0 = x_rep;
    options = optimset('Algorithm','sqp','Display','off');
    [x,fval,exitflag,~] = fmincon(fun,x0,[],[],[],[],xrep_ll,xrep_ul,mycon,options);
    xrep_new = x;
    if exitflag == -2
        error("No feasible repair point found")
    elseif exitflag == 0
        error("Exceeded number of function evals. solution not found")
    end
    exitflag;
end

function score = obj(x_rep, x_ind)
% Renaming variables
Ry = x_ind(1); % Y axis radius of eliptical crown
Rc = x_ind(2); % Crimped radius of the stent
Sr = x_ind(3); % Radius of the circular strut
SL = x_ind(4);% Stent length
Nn = x_ind(5);% % Total number of unit struts in the main helix of the stent - including the last 3/4 strut
alpha = x_rep(1); % Helix angle - repair variable.
SLt = x_rep(2);% Length of the transistion part of the stent
beta = x_rep(3);

% Parameters for main helix 
Nu = 5; % Number of unit strut components in a unit strut
Nr = 8; % Number of complete unit strut components in a full 360 degree circle
Nt = Nr - Nu;% Number of complete strut components between bottom and top intersecting crowns
welded_crown_shift_factor = 0.2;% Shift in the welded crown as a proportion of Sr
Nt2 = 2; 
Nt3 = 3;
Nt1 = Nr - Nt2 - Nt3 + 1; 

% Calculating dependent variables 
Rys = 1.1*Sr + Ry;% Radius of the intersecting crowns
SSy = ((SL - SLt)*Nr + Rys*(0 - 2*Nu*Nn + 4))/((Nr + 1)*(Nu*Nn - 2) - (Nr)*(Nu*Nn - 1)); % Height of short strut
TSy = (2*Rys + (Nr + 1)*SSy)/Nr;
Rx = 0.5*((2*pi*Rc) - ((TSy - SSy)*Nr/tand(alpha)));% X radius of elliptical crown
TSx = 4*Rx - ((TSy - SSy)/tand(alpha));% Width of tall strut
AAX = (4*Rx - TSx)*Nu;% X distance covered by a strut
AAY = (TSy - SSy)*Nu;% y distance covered by a strut
cr_yscale = Ry/Rx;
cr_yscales = Rys/Rx;
Sy0 = 0.5*(SLt) -2*AAY - 3*Rys - (4*Rx - TSx)*(Nt1 + Nt2)*tand(beta) - TSy;
%Rys_shifted = Rys - welded_crown_shift_factor*Sr;
%cr_yscales_shifted = Rys_shifted/Rx;

score = abs(Sy0 - SSy);
end


function [c, ceq] = cons(x_rep, x_ind, cr_yscale_limit,cr_yscales_limit)

% Renaming variables
% Independent variables
Ry = x_ind(1); % Y axis radius of eliptical crown
Rc = x_ind(2); % Crimped radius of the stent
Sr = x_ind(3); % Radius of the circular strut
SL = x_ind(4);% Stent length
Nn = x_ind(5);% % Total number of unit struts in the main helix of the stent - including the last 3/4 strut

% Repair variables
alpha = x_rep(1); % Helix angle - repair variable.
SLt = x_rep(2);% Length of the transistion part of the stent
beta = x_rep(3);

% Parameters for main helix 
Nu = 5; % Number of unit strut components in a unit strut
Nr = 8; % Number of complete unit strut components in a full 360 degree circle
Nt = Nr - Nu;% Number of complete strut components between bottom and top intersecting crowns
welded_crown_shift_factor = 0.2;% Shift in the welded crown as a proportion of Sr
% Parameters - For the top strut
Nt2 = 2; 
Nt3 = 3;
Nt1 = Nr - Nt2 - Nt3 + 1; 

% Calculating dependent variables 
Rys = 1.1*Sr + Ry;% Radius of the intersecting crowns
SSy = ((SL - SLt)*Nr + Rys*(0 - 2*Nu*Nn + 4))/((Nr + 1)*(Nu*Nn - 2) - (Nr)*(Nu*Nn - 1)); % Height of short strut
TSy = (2*Rys + (Nr + 1)*SSy)/Nr;
Rx = 0.5*((2*pi*Rc) - ((TSy - SSy)*Nr/tand(alpha)));% X radius of elliptical crown
TSx = 4*Rx - ((TSy - SSy)/tand(alpha));% Width of tall strut
AAX = (4*Rx - TSx)*Nu;% X distance covered by a strut
AAY = (TSy - SSy)*Nu;% y distance covered by a strut
cr_yscale = Ry/Rx;
cr_yscales = Rys/Rx;

% Top and Transition related dependent variables
Sy0 = 0.5*(SLt) -2*AAY - 3*Rys - (4*Rx - TSx)*(Nt1 + Nt2)*tand(beta) - TSy;

% Welded crown shift based dependent variables
%Rys_shifted = Rys - welded_crown_shift_factor*Sr;
%cr_yscales_shifted = Rys_shifted/Rx;

xdep(1) = Rys;
xdep(2) = SSy;
xdep(3) = TSy;
xdep(4) = Rx;
xdep(5) = TSx;
xdep(6) = AAX;
xdep(7) = AAY;
xdep(8) = cr_yscale;
xdep(9) = cr_yscales;
xdep(10) = Sy0;
%xdep(11) = Rys_shifted;
%xdep(12) = cr_yscales_shifted;

% Specifying constraints
c = -1.*(xdep);
c(11) = xdep(2) - xdep(3);
c(12) = -1*xdep(8) + cr_yscale_limit(1);
c(13) = xdep(8) - cr_yscale_limit(2);
c(14) = -1*xdep(9) + cr_yscales_limit(1);
c(15) = xdep(9) - cr_yscales_limit(2);
c(16) = 0.20*xdep(2) - xdep(10);
c(17) = xdep(5) - 0.5*xdep(4);

ceq = [];
end
