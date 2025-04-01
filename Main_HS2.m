clc;
close all;
clear all;

%% This code generates the HS2 stent based in the design variables x values
% Design variable x values to be entred below:
x(1) = 0.13; % Y axis radius of eliptical crown
x(2) = 0.4570; % Crimped radius of the stent
x(3) = 0.04; % Radius of the circular strut
x(4) = 15.2;% Stent length
x(5) = 16;% % Total number of unit struts in the main helix of the stent - including the last 3/4 strut

% Add hs2 folder in the path and call the function for stent generation
hs2_folder = fullfile(pwd, 'HS2');
addpath(hs2_folder);
Main(x)
rmpath(hs2_folder);