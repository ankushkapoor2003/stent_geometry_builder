function Main(x_c)
% This function takes in the x values and saves the HS1 stent geometry
% Saving the current wd path
original_dir = pwd;

% Setting up the results directory
result_dir_name = 'HS1_Results';
result_dir_path = fullfile(original_dir, result_dir_name);

% Generating the results directory if not already created
if ~exist(result_dir_path, 'dir')
    mkdir(result_dir_path);
end

% Fixed parameters for resolution
distance_thresh = 1E-5; % 1 micron distance threshold between adjacent points 
distance_thresh_conn = 1E-5;
distance_thresh_er = 1E-4;
resolution_threshold = [distance_thresh distance_thresh_conn distance_thresh_er];

% Base values of x. 
x_b(1) = 19;% Number of unit struts in a helix
x_b(2) = 0.46;% Crimped radius
x_b(3) = 38; % Helix angle alpha x(3) - Repair Variable
x_b(4) = 3;% Approximate number of units in a rotation, k % x(4)
x_b(5) = 18;% Stent length SL x(5)
x_b(6) = 0.95; %End Row height (6)
x_b(7) = 0.1; % Width of the S Connector
x_b(8) = 1.63339;% ww_ab(1)
x_b(9) = 0.39003;% ww_ab(2) 
x_b(10) = 0.75959;% ww_ab(3)
x_b(11) = 0.81258;% ww_ab(4)
x_b(12) = 1.31384;% ww_da(1)
x_b(13) = 0.04872;% ww_da(2)
x_b(14) = 0.00164;% ww_da(3)
x_b(15) = 0.07540;% ww_connector(1)
x_b(16) = 0.61571;% ww_connector(2)
x_b(17) = 4.16285;% ww_connector(3)
x_b(18) = 1.79921;% ww_connector(4)
x_b(19) = 0.06;% strut width
x_b(20) = 0.06;% strut_thickness
x_b(21) = 0.3;% 0.3 was default cp_width_ab WAB - Repair Variable
x_b(22) = 0.125; % x(22) cp_width_da WDAN - Repair Variable
x_b(23) = 0.675;% x(23) cp_height_da HDAN - Repair Variable

% Separating into base Variables to be be optimised and fixed
[~, x_fixed] = Divide_xx(x_b);

% Setting up the input for HS1 stent generation
x_ind = x_c;
x_ind = round(x_ind,5);
x_ind(:,1) = round(x_ind(:,1));
x = Recreate_xx(x_ind, x_fixed);

% Rounding the integer variables before sending to geometry creator
x(1) = round(x(1)); x(4) = round(x(4));

% Repair and raise error if repair fails
try
    xrep_old = x([3,21:23]);
    x_rep = Repair(x([1:2,4:7]), xrep_old);
    x([3,21:23]) = x_rep;
catch
    error('Error:HS1', 'Repair failed for HS1 stent. Kindly try with different variables');
end

% Trying to generate stent with increasingly relaxed threshold, then
% generating the stl, stent points and evaluating distance between target
% and candidate to report distance score. 
try
    try
        HS1_Generation(x, resolution_threshold, result_dir_path);
    catch
        py.sld_interface.exit_sw();
        resolution_threshold(1) = resolution_threshold(1)*10;
        try
            HS1_Generation(x, resolution_threshold, result_dir_path);
        catch
            py.sld_interface.exit_sw();
            resolution_threshold(1) = resolution_threshold(1)*10;
            HS1_Generation(x, resolution_threshold, result_dir_path);
        end
    end
catch
    error('Error:HS2', "HS2 stent geometry construction failed. Kindly try with different design variables")
end
return