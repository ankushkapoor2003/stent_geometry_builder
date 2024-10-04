clc;
close all;
clear all;

%% This code generates the IR1 stent based in the design variables x values
% Design variable x values to be entred below:
x(1) = 0.07; % Row 1 NURBS Control Points Width 
x(2) = 0.8; % Row 1 NURBS Control Points Height
x(3) = 0.61317;%ww(1)
x(4) = 0.09022;%ww(2)
x(5) = 0.86459;%ww(3)
x(6) = 0.36024;%ww(4)
x(7) = 0.458; %Stent Crimped radius
x(8) = 0.9; % Row 2 Tall Strut height 
x(9) = 0.7; %Row 2 short strut height 
x(10) = 21/17; %Row 2 Short Strut top crown radius / Row 2 other crown radius
x(11) = 0.095; % Connector lateral crown radius
x(12) = 0.16; % Connector lateral length
x(13) = 3; % Total Number of repittions of SS-TS-C unit
x(14) = 18; % Stent Length
x(15) = 13; % Number of rows in the stent
x(16) = 1;% Stent width factor (1) will be multiplied with base value before sending to the geometry creater
x(17) = 1;% Stent thicknesss factor (1) will be multiplied with base value before sending to the geometry creater

% Add ir1 folder in the path and call the function for stent generation
ir1_folder = fullfile(pwd, 'IR1');
addpath(ir1_folder);
Main(x)
rmpath(ir1_folder);