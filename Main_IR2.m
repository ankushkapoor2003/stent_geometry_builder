clc;
close all;
clear all;

%% This code generates the IR2 stent based in the design variables x values
% Design variable x values to be entred below:
x(1) = 18; % Stent length (SL)
x(2) = 0.46; % Crimped radius of the stent (Rc)
x(3) = 0.9; % Strut length (Excluding the crown part) (SRL)
x(4) = 6; % Number of unit struts in each ring (Nu)
x(5) = 0.08; % Strut thickness (including the connector thickness) (St)
x(6) = 0.08; % Strut width (including the connector width) (Sw)
x(7) = 9; % Number of rings in a stent (Nr)
x(8) = 0.12; % Connector shape variable 1 (100% variationa allowed p1
x(9) = 0.5; % Connector shape variable 2 (100%) variation allowed here as well p2

% Add ir2 folder in the path and call the function for stent generation
ir2_folder = fullfile(pwd, 'IR2');
addpath(ir2_folder);
Main(x)
rmpath(ir2_folder);