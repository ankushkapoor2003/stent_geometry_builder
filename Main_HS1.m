clc;
close all;
clear all;

%% This code generates the IR1 stent based in the design variables x values
% Design variable x values to be entred below:
% Base values of x. 
x(1) = 19;% Number of unit struts in a helix
x(2) = 0.46;% Crimped radius
x(3) = 18;% Stent length SL x(5)
x(4) = 0.95; %End Row height (6)
x(5) = 0.1; % Width of the S Connector
x(6) = 1.63339;% ww_ab(1)
x(7) = 0.39003;% ww_ab(2) 
x(8) = 0.75959;% ww_ab(3)
x(9) = 0.81258;% ww_ab(4)
x(10) = 1.31384;% ww_da(1)
x(11) = 0.04872;% ww_da(2)
x(12) = 0.00164;% ww_da(3)
x(13) = 0.07540;% ww_connector(1)
x(14) = 0.61571;% ww_connector(2)
x(15) = 4.16285;% ww_connector(3)
x(16) = 1.79921;% ww_connector(4)
x(17) = 0.06;% strut width
x(18) = 0.06;% strut_thickness

% Add hs1 folder in the path and call the function for stent generation
hs1_folder = fullfile(pwd, 'HS1');
addpath(hs1_folder);
Main(x)
rmpath(hs1_folder);