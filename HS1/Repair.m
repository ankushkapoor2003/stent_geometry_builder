function [xrep_new] = Repair(x_ind, x_rep)

    % Evaluating dependent variables for received design
    x_dep = CalculateDependentVariables(x_rep, x_ind); 

    % Evaluating bounds on repair and dependent variables
    xrep_ll = 0.7*x_rep;
    if x_ind(3) == 5
        xrep_ll = 0.5*x_rep;
    end
    xrep_ul = 1.3*x_rep;
    cr_scale_limit(1) = 0.5;
    cry_scale_limit(2) = 1.5;
    min_radius_limit = 0.04;

    % The objective function for repair aims to enforce constraints and minimise the difference between connector height and height of AB Strut 
    % Objective function,constraint evaluation, and fmincon
    fun = @(x_rep)obj(x_rep, x_ind);
    mycon = @(x_rep)cons(x_rep, x_ind, cry_scale_limit, min_radius_limit);
    x0 = x_rep;
    options = optimset('Algorithm','sqp','Display','off');
    [x,fval,exitflag,~] = fmincon(fun,x0,[],[],[],[],xrep_ll,xrep_ul,mycon,options);
    xrep_new = x;
    if exitflag == -2
        error("No feasible repair point found")
    end
end

function score = obj(x_rep, x_ind)
% This function computes the IGD with target
% Renaming variables
NS = x_ind(1); % Number of unit struts in a helix x(1)
crimped_radius = x_ind(2);% x(2)
alpha = x_rep(1); %Helix angle x(4)
k = x_ind(3);% Approximate number of units in a  % x(5)
SL = x_ind(4);% Stent length x(6)
L1 = x_ind(5); %End ROw height (7)
wc = x_ind(6); %x(8) Width of the S Connector
cp_width_ab = x_rep(2);
cp_width_da = x_rep(3); % 
cp_height_da = x_rep(4);% 

% Alternate notations used in formulae
R = crimped_radius;
WDAN = cp_width_da;
WAB = cp_width_ab;
HDAN = cp_height_da;

% Calculating dependent variables 
NS2 = (2*NS+1)/2;
HH = SL - 2*L1; % xdep(1)
crx = (1/(4*k))*(2*pi()*R - (k - 1)*WDAN -2*wc - 2*k*WAB); % xdep(2)
AAX = 2*WAB + 4*crx + WDAN; % xdep(3)
HAB = -1*NS2*AAX*tand(alpha) + (SL - 2*L1); % xdep(4)
cry = 0.25*(2*HAB -2*HDAN -AAX*tand(alpha));% xdep(5)
cr_yscale = cry/crx; % xdep(6)
AAY = tand(alpha)*AAX; %xdep(7)
hc = (k*AAY/2);% xdep(8)
score = abs(HAB - hc);
end

function [c, ceq] = cons(x_rep, x_ind, yscale_limit,min_radius_limit)
% Renaming variables
NS = x_ind(1); % Number of unit struts in a helix x(1)
crimped_radius = x_ind(2);% x(2)
alpha = x_rep(1); %Helix angle x(4)
k = x_ind(3);% Approximate number of units in a  % x(5)
SL = x_ind(4);% Stent length x(6)
L1 = x_ind(5); %End ROw height (7)
wc = x_ind(6); %x(8) Width of the S Connector
cp_width_ab = x_rep(2);
cp_width_da = x_rep(3); % 
cp_height_da = x_rep(4);% 

% Alternate notations used in formulae
R = crimped_radius;
WDAN = cp_width_da;
WAB = cp_width_ab;
HDAN = cp_height_da;

% Calculating dependent variables 
NS2 = (2*NS+1)/2;
HH = SL - 2*L1; % xdep(1)
crx = (1/(4*k))*(2*pi()*R - (k - 1)*WDAN -2*wc - 2*k*WAB); % xdep(2)
AAX = 2*WAB + 4*crx + WDAN; % xdep(3)
HAB = -1*NS2*AAX*tand(alpha) + (SL - 2*L1); % xdep(4)
cry = 0.25*(2*HAB -2*HDAN -AAX*tand(alpha));% xdep(5)
cr_yscale = cry/crx; % xdep(6)
AAY = tand(alpha)*AAX; %xdep(7)
hc = (k*AAY/2);% xdep(8)
xdep(1) = HH;
xdep(2) = crx;
xdep(3) = AAX;
xdep(4) = HAB;
xdep(5) = cry;
xdep(6) = cr_yscale;
xdep(7) = AAY;
xdep(8) = hc;
% Specifying constraints
c = -1.*(xdep);
c(2) = c(2) + min_radius_limit;
c(5) = c(5) + min_radius_limit;
c(9) = -1*xdep(6) + yscale_limit(1);
c(10) = xdep(6) - yscale_limit(2);
ceq = [];
end


function xdep = CalculateDependentVariables(x_rep, x_ind)
% Renaming variables
NS = x_ind(1); % Number of unit struts in a helix x(1)
crimped_radius = x_ind(2);% x(2)
alpha = x_rep(1); %Helix angle x(4)
k = x_ind(3);% Approximate number of units in a  % x(5)
SL = x_ind(4);% Stent length x(6)
L1 = x_ind(5); %End ROw height (7)
wc = x_ind(6); %x(8) Width of the S Connector
cp_width_ab = x_rep(2);
cp_width_da = x_rep(3); % 
cp_height_da = x_rep(4);% 

% Alternate notations used in formulae
R = crimped_radius;
WDAN = cp_width_da;
WAB = cp_width_ab;
HDAN = cp_height_da;

% Calculating dependent variables 
NS2 = (2*NS + 1)/2;
HH = SL - 2*L1; % xdep(1)
crx = (1/(4*k))*(2*pi()*R - (k - 1)*WDAN -2*wc - 2*k*WAB); % xdep(2)
AAX = 2*WAB + 4*crx + WDAN; % xdep(3)
HAB = -1*NS2*AAX*tand(alpha) + (SL - 2*L1); % xdep(4)
cry = 0.25*(2*HAB -2*HDAN -AAX*tand(alpha));% xdep(5)
cr_yscale = cry/crx; % xdep(6)
AAY = tand(alpha)*AAX; %xdep(7)
hc = (k*AAY/2);% xdep(8)
xdep(1) = HH;
xdep(2) = crx;
xdep(3) = AAX;
xdep(4) = HAB;
xdep(5) = cry;
xdep(6) = cr_yscale;
xdep(7) = AAY;
xdep(8) = hc;
end
