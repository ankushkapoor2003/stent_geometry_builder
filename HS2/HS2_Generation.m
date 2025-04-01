function [helix_sw_name] = HS2_generation(x_c, result_location)


%% Generate x for candidate generation 
% Base values
x_b(1) = 0.13; % Y axis radius of eliptical crown
x_b(2) = 0.4570; % Crimped radius of the stent
x_b(3) = 0.04; % Radius of the circular strut
x_b(4) = 15.2;% Stent length
x_b(5) = 16;% % Total number of unit struts in the main helix of the stent - including the last 3/4 strut
x_b(6) = 23.5; % Helix angle - repair variable.
x_b(7) = 6.7;% Length of the transistion part of the stent - Repair variable
x_b(8) = 11.75; % Beta angle of the bottom of the first row

% Separating into base Variables to be be optimised and fixed
[~, x_fixed] = Divide_xx(x_b);

% Setting up the input for HSs stent generation
x_ind = x_c;
x_ind = round(x_ind,5);
x_ind(:,5) = round(x_ind(:,5));
x = Recreate_xx(x_ind, x_fixed);

% Perform Repair
xrep_old = x([6:8]);
x_rep = Repair(x([1:5]), xrep_old);
x([6:8]) = x_rep;

%% Start Solidworks
py.sld_interface.start();

%% Start stent building
Ry = x(1);% Y axis radius of eliptical crown
Rc = x(2); % Crimped radius of the stent 
Sr = x(3); % Radius of the circular strut
SL = x(4); % Stent length
Nn = x(5); % Total number of unit struts in the main helix of the stent - including the 
alpha = x(6);% Helix angle To be converted to repair variable.
SLt = x(7);% Length of the transistion part of the stent
beta = x(8);

% Parameters - Change in any of these parameters should correspond with the
% change in the same parameters in all functions of Repair - IMPORTANT
Nu = 5; % Number of unit strut components in a unit strut
Nr = 8; % Number of complete unit strut components in a full 360 degree circle
Nt2 = 2; %
Nt3 = 3;%
welded_crown_shift_factor = 0.2;% Shift in the welded crown as a proportion of Sr

% Solidworks related parameters 
distance_thresh_main = 0.03; % Primary for main helix
distance_thresh_transition = 0.03; % secondary for non main helix components
distance_thresh_row2 = 0.03;
distance_thresh_row1 = 0.03;

% Dependent variables
Rys = 1.1*Sr + Ry;% Radius of the intersecting crowns
Rys_shifted = Rys - welded_crown_shift_factor*Sr;
SSy = ((SL - SLt)*Nr + Rys*(0 - 2*Nu*Nn + 4))/((Nr + 1)*(Nu*Nn - 2) - (Nr)*(Nu*Nn - 1)); % Height of short strut
TSy = (2*Rys + (Nr + 1)*SSy)/Nr;
Rx = 0.5*((2*pi*Rc) - ((TSy - SSy)*Nr/tand(alpha)));% X radius of elliptical crown
TSx = 4*Rx - ((TSy - SSy)/tand(alpha));% Width of tall strut
Nt = Nr - Nu;% Number of complete strut components between bottom and top intersecting crowns
AAX = (4*Rx - TSx)*Nu;% X distance covered by a strut
AAY = (TSy - SSy)*Nu;% y distance covered by a strut
cr_yscales_shifted = Rys_shifted/Rx;

% Top and Transition related dependent variables
Nt1 = Nr - Nt2 - Nt3 + 1; 
Sy0 = 0.5*(SLt) -2*AAY - 3*Rys - (4*Rx - TSx)*(Nt1 + Nt2)*tand(beta) - TSy;
cr_yscale = Ry/Rx;
cr_yscales = Rys/Rx;
main_helix_shift_initial_x = 2*pi()*Rc + (Nt1 + Nt2)*(4*Rx - TSx) + Rx - TSx + Rx + 2*AAX;% 2AAY due to transisiton helix, which is extra to top helix 
main_helix_shift_initial_y = -0.5*SLt - Rx*cr_yscales; 

% Stent generation parameters
cut_off_distance_for_struts = 0.02; %Can be varied to see if model works better
np_strut = 500;
np_crown = 200;

% Generating the tall strut
slope_correction_cutoff_value_TS = 0.5*TSy - cut_off_distance_for_struts;
ts_strut = tall_strut_generator(TSx, TSy, np_strut);
ts_slopes = strut_slopes(ts_strut, np_strut);

% Code 4 tall strut
[ts_corrected_strut, ts_corrected_slope] = strut_slope_correction(ts_strut, ts_slopes, slope_correction_cutoff_value_TS);

% Generatig the short strut
% Code 2: Short Strut
short_strut = short_strut_generator(SSy, np_strut);

% generating the parsec points for top crowns 
z_te = ts_corrected_strut(end,1) + Rx - 0.5*TSx;
x_te = -(ts_corrected_strut(end,2) - 0.5*TSy - Rx);
z_te_slope = -1/ts_corrected_slope(end,1);
parsec_coeffs_top_A = parsec(Rx, x_te, z_te, z_te_slope);
crown_points_parsec_top = generate_crown_parsec(parsec_coeffs_top_A, "top", np_crown/4,x_te,Rx, cr_yscale);

% generationg the semicircular crown points for top
crown_points_circular = crown_top_circular(Rx, 4*np_crown);

% Code - 3: Crown top
crown_top = crown_aggregator_V3(crown_points_parsec_top, crown_points_parsec_top, crown_points_circular,'standard', cr_yscale);

% generating top welded crown
crown_points_parsec_top_w = generate_crown_parsec(parsec_coeffs_top_A, "top", np_crown/4,x_te,Rx, cr_yscales);

% Code - 13: Crown top welded
crown_top_w = crown_aggregator_V4(crown_points_parsec_top_w, crown_points_parsec_top_w, crown_points_circular,'shifted', cr_yscales_shifted, cr_yscales, Rx);

% generating the parsec points for bottom crown
z_te_b = -ts_corrected_strut(1,1) + Rx - 0.5*TSx;
x_te_b = (ts_corrected_strut(1,2) + 0.5*TSy + Rx);
z_te_slope_b = -1/ts_corrected_slope(1,1);
parsec_coeffs_bottom = parsec(Rx, x_te_b, z_te_b, z_te_slope_b);
crown_points_parsec_bottom = generate_crown_parsec(parsec_coeffs_bottom, "bottom", np_crown/4,x_te_b, Rx, cr_yscale);
crown_points_circular_bottom = crown_top_circular(Rx, 4*np_crown);
crown_bottom_rev = crown_aggregator_V3(crown_points_parsec_bottom, crown_points_parsec_bottom, crown_points_circular_bottom,'standard',cr_yscale);

% Code 1 Crown Bottom
crown_bottom = [-1.*crown_bottom_rev(:,1) -1.*crown_bottom_rev(:,2)];

% Generating the welded crown bottom
crown_points_parsec_bottom_w = generate_crown_parsec(parsec_coeffs_bottom, "bottom", np_crown/4,x_te_b, Rx, cr_yscales);
crown_bottom_rev_w = crown_aggregator_V4(crown_points_parsec_bottom_w, crown_points_parsec_bottom_w, crown_points_circular_bottom,'shifted', cr_yscales_shifted, cr_yscales, Rx);

% Code 11 Crown Bottom
crown_bottom_w = [-1.*crown_bottom_rev_w(:,1) -1.*crown_bottom_rev_w(:,2)];

% Code 34 TS with Parsec 
Parsec_top_shifted(:,1) = crown_points_parsec_top(:,2) - Rx + 0.5*TSx;
Parsec_top_shifted(:,2) = crown_points_parsec_top(:,1) + Rx*cr_yscale + 0.5*TSy;
Parsec_bottom_shifted = -1*Parsec_top_shifted;
ts_full_strut = [Parsec_bottom_shifted;ts_corrected_strut; flipud(Parsec_top_shifted)];

% Code 31 for crown bottom without parsec
crown_bottom_nopar_rev = crown_aggregator_V3([], [], crown_points_circular_bottom,'standard',cr_yscale);
crown_bottom_nopar = [-1.*crown_bottom_nopar_rev(:,1) -1.*crown_bottom_nopar_rev(:,2)];

% Code 41 for welded crown bottom without parsec
crown_bottom_nopar_rev_w = crown_aggregator_V3([], [], crown_points_circular_bottom,'standard',cr_yscales);
crown_bottom_nopar_w = [-1.*crown_bottom_nopar_rev_w(:,1) -1.*crown_bottom_nopar_rev_w(:,2)];

% Code 33 for crown top without parsec
crown_top_nopar = crown_aggregator_V3([], [], crown_points_circular,'standard',cr_yscale);

% Code 43 for welded crown bottom without parsec
crown_top_nopar_w = crown_aggregator_V3([], [], crown_points_circular,'standard',cr_yscales);

% Generating unit strut
unit_config = [11 2 3 4 1 2 3 4 1 2 3 4 1 2 13 4 1 2 3 4];
unit_feature_pos = feature_pos_unit(unit_config, Rx, cr_yscale, cr_yscales, TSx, TSy, SSy);
%unit_feature_pos_reverse = feature_pos_unit_reverse(unit_config, Rx, cr_yscale, cr_yscales, TSx, TSy, SSy);

% Generate unit config for transition strut
config_transition = [11 2 3 4 1 2 3 4 1 2 13 4 1 2 3 4 1 2 3 4 11 2 3 4 1 2 13 4 1 2 3 4 1 2 13 4 1 2 3 4];% Will have to be automated later based on the total number of Nr and Nu. For now, fixing it to generate a MVP
feature_pos_transition = feature_pos_unit(config_transition, Rx, cr_yscale, cr_yscales, TSx, TSy, SSy);

% Generate unit config for row_2
config_top2 = [43 34 41 32 33 34 31 32 33 34 31 32 33 34 31 32 43 34 41 32 33 34 31 32 43 34];% Will have to be automated later based on the total number of Nr and Nu. For now, fixing it to generate a MVP
[feature_pos_top2, TSy_r2, SSy_r2] = feature_pos_unit_top2(config_top2, Rx, cr_yscale, cr_yscales, TSx, TSy, SSy, alpha, beta);

% Generate unit config for row_2
config_top1 = [41 32 43 34 31 32 43 34 31 32 43 34 31 32 43 34 41 32 43 34 31 32 43 34 41 32 43 34 31 32 43 34 31 32];% Will have to be automated later based on the total number of Nr and Nu. For now, fixing it to generate a MVP
tantheta = (Sy0 + 2*Rx*cr_yscales)/(2*pi*Rc - (Nt1 + Nt2)*(4*Rx - TSx) - 2*Rx);
[feature_pos_top1, TSy_r1, SSy_r1] = feature_pos_unit_top1(config_top1, Rx, cr_yscale, cr_yscales, TSx, TSy, SSy, alpha, beta, tantheta, Sy0);

%% Building helical stent 
% Generate main helix of the stent
[helix_main, helix_sw_name, main_helix_connection_point] = helix_builder(ts_corrected_strut, short_strut, crown_top, crown_bottom, crown_top_w, crown_bottom_w, cr_yscale, unit_feature_pos, unit_config, AAX, AAY, Nn, main_helix_shift_initial_x, main_helix_shift_initial_y, Rc, Sr, distance_thresh_main);
[transition_strut, transition_strut_rev, helix_sw_name, unit_transition_connection_point] = transition_builder(ts_corrected_strut, short_strut, crown_top, crown_bottom, crown_top_w, crown_bottom_w, cr_yscale, feature_pos_transition, config_transition, AAX, AAY, Nn, main_helix_shift_initial_x, main_helix_shift_initial_y, TSy, SSy, TSx, Rx, cr_yscales, alpha, SL, Sy0, beta, Nt1, Nt2, Rc, Sr, distance_thresh_transition, main_helix_connection_point, helix_sw_name);
[row_2, row_2_rev, helix_sw_name, row2_connection_point] =row_2_builder(ts_full_strut, short_strut, crown_top_nopar, crown_bottom_nopar, crown_top_nopar_w, crown_bottom_nopar_w, cr_yscale, feature_pos_top2, config_top2, AAX, AAY, Nn, main_helix_shift_initial_x, main_helix_shift_initial_y, TSy, SSy, TSx, Rx, cr_yscales, alpha, beta, TSy_r2, SSy_r2, Nt1, Nt2, Sy0, Rc, SL, Sr, distance_thresh_row2, unit_transition_connection_point, helix_sw_name);
[row_1, row_1_rev, helix_sw_name] =row_1_builder(ts_full_strut, short_strut, crown_top_nopar, crown_bottom_nopar, crown_top_nopar_w, crown_bottom_nopar_w, cr_yscale, feature_pos_top1, config_top1, AAX, AAY, Nn, main_helix_shift_initial_x, main_helix_shift_initial_y, TSy, SSy, TSx, Rx, cr_yscales, alpha, beta, TSy_r1, SSy_r1, Nt1, Nt2, Sy0, Rc, SL, Sr, distance_thresh_row1, row2_connection_point, helix_sw_name);
%% Rotate to orient stent about -y axis so that it matches with the other
% stents 
%helix_sw_name = py.sld_interface.reorient_ro(helix_sw_name);% Can be commented if the user would not like to reorient the generated stent

%% Saving stent solidworks file, parasolid file and Mesh STL file 
% Saving 3D stent files
file_name = "stent.SLDPRT";
ps_file_name = "stent.x_t";
stl_file_name = "stent.STL";
fullpath_sw = fullfile(result_location,file_name);
py.sld_interface.save_sw_file(fullpath_sw);
fullpath_ps = fullfile(result_location,ps_file_name);
fullpath_stl = fullfile(result_location,stl_file_name);
py.sld_interface.exit_sw();
py.sld_interface.start();
py.sld_interface.sw_to_ps_file(fullpath_sw, fullpath_ps);
py.sld_interface.ps_to_stl_single_helix(fullpath_ps,fullpath_stl);
py.sld_interface.exit_sw();
% Stent built. 
end

%% Helper functions to generate geometries
function [row_1, row_1_rev, helix_sw_name] = row_1_builder(ts_full_strut, short_strut, crown_top, crown_bottom, crown_top_w, crown_bottom_w, cr_yscale, unit_feature_pos, unit_config, AAX, AAY, Nn, main_helix_shift_initial_x, main_helix_shift_initial_y, TSy, SSy, TSx, Rx, cr_yscales, alpha, beta, TSy_r1, SSy_r1, Nt1, Nt2, Sy0, Rc, SL, Sr, distance_thresh, row2_connection_point, helix_sw_name)
helix_shift_x = 0;
helix_shift_y = -1*Rx*cr_yscales - Sy0 - Rx*cr_yscales;
helix = [];
last_unit_flag = 0;
[row_1, row_1_rev, row1_sw_name, row1_rev_sw_name] = generate_row1_unit_strut(ts_full_strut, short_strut, crown_top, crown_bottom, crown_top_w, crown_bottom_w, cr_yscale, unit_feature_pos, unit_config, helix_shift_x, helix_shift_y, Nn, AAX, AAY, TSy, SSy, TSx, Rx, cr_yscales, alpha, beta, TSy_r1, SSy_r1, Rc, Nt1, Nt2, SL, Sy0, Sr, distance_thresh, row2_connection_point, helix_sw_name);
helix_sw_name = py.sld_interface.row_combine(helix_sw_name, row1_sw_name);% I created this function to combine rows in IR1 and reusing it. That is why the name of the function is row_combine. It just combines solid bodies
helix_sw_name = py.sld_interface.row_combine(helix_sw_name, row1_rev_sw_name);% I created this function to combine rows in IR1 and reusing it. That is why the name of the function is row_combine. It just combines solid bodies
end

function [row_2, row_2_rev, helix_sw_name, row2_connection_point] = row_2_builder(ts_full_strut, short_strut, crown_top, crown_bottom, crown_top_w, crown_bottom_w, cr_yscale, unit_feature_pos, unit_config, AAX, AAY, Nn, main_helix_shift_initial_x, main_helix_shift_initial_y, TSy, SSy, TSx, Rx, cr_yscales, alpha, beta, TSy_r2, SSy_r2, Nt1, Nt2, Sy0, Rc, SL, Sr, distance_thresh, unit_transition_connection_point, helix_sw_name);
helix_shift_x = 2*pi*Rc;
helix_shift_y = -2*Rx*cr_yscales - Sy0;
helix = [];
last_unit_flag = 0;
[row_2, row_2_rev, row2_sw_name, row2_rev_sw_name, row2_connection_point] = generate_row2_unit_strut(ts_full_strut, short_strut, crown_top, crown_bottom, crown_top_w, crown_bottom_w, cr_yscale, unit_feature_pos, unit_config, helix_shift_x, helix_shift_y, Nn, AAX, AAY, TSy, SSy, TSx, Rx, cr_yscales, alpha, beta, TSy_r2, SSy_r2, Nt1, Nt2, Rc, SL, Sy0, Sr, distance_thresh, unit_transition_connection_point,  helix_sw_name);
helix_sw_name = py.sld_interface.row_combine(helix_sw_name, row2_sw_name);% I created this function to combine rows in IR1 and reusing it. That is why the name of the function is row_combine. It just combines solid bodies
helix_sw_name = py.sld_interface.row_combine(helix_sw_name, row2_rev_sw_name);% I created this function to combine rows in IR1 and reusing it. That is why the name of the function is row_combine. It just combines solid bodies
end

function [transition_strut, transition_strut_rev, helix_sw_name, unit_transition_connection_point] = transition_builder(ts_corrected_strut, short_strut, crown_top, crown_bottom, crown_top_w, crown_bottom_w, cr_yscale, unit_feature_pos, unit_config, AAX, AAY, Nn, main_helix_shift_initial_x, main_helix_shift_initial_y, TSy, SSy, TSx, Rx, cr_yscales, alpha, SL, Sy0, beta, Nt1, Nt2, Rc, Sr, distance_thresh, main_helix_connection_point, helix_sw_name)
helix_shift_x = main_helix_shift_initial_x - 2*AAX;
helix_shift_y = main_helix_shift_initial_y + 2*AAY;
helix = [];
last_unit_flag = 0;
[transition_strut, transition_strut_rev, transition_strut_sw_name, transition_strut_rev_sw_name, unit_transition_connection_point]= generate_transition_unit_strut(ts_corrected_strut, short_strut, crown_top, crown_bottom, crown_top_w, crown_bottom_w, cr_yscale, unit_feature_pos, unit_config, helix_shift_x, helix_shift_y, Nn, AAX, AAY, TSy, SSy, TSx, Rx, cr_yscales, alpha, SL, Sy0, beta, Nt1, Nt2, Rc, Sr, distance_thresh, main_helix_connection_point, helix_sw_name);
helix_sw_name = py.sld_interface.row_combine(helix_sw_name, transition_strut_sw_name);% I created this function to combine rows in IR1 and reusing it. That is why the name of the function is row_combine. It just combines solid bodies
helix_sw_name = py.sld_interface.row_combine(helix_sw_name, transition_strut_rev_sw_name);% I created this function to combine rows in IR1 and reusing it. That is why the name of the function is row_combine. It just combines solid bodies
end

function [helix_main, helix_sw_name, main_helix_connection_point] = helix_builder(ts_corrected_strut, short_strut, crown_top, crown_bottom, crown_top_w, crown_bottom_w, cr_yscale, unit_feature_pos, unit_config, AAX, AAY, Nn, main_helix_shift_initial_x, main_helix_shift_initial_y, Rc, Sr, distance_thresh)
helix_shift_x = main_helix_shift_initial_x;
helix_shift_y = main_helix_shift_initial_y;
rel_z_shift = -1*AAY;
rel_theta_shift = 1*(AAX/Rc);
helix = [];
last_unit_flag = 0;
py.sld_interface.create_z_axis();% Create cylinder axis
for i=1:Nn
    if i == Nn
        last_unit_flag = 1;
    end
    if i == 1 || i == Nn
        sw_unit_strut_flag = 1;
    else 
        sw_unit_strut_flag = 0;
    end
    [unit_strut, unit_strut_sw_name, unit_strut_connection_point] = generate_unit_strut(ts_corrected_strut, short_strut, crown_top, crown_bottom, crown_top_w, crown_bottom_w, cr_yscale, unit_feature_pos, unit_config, helix_shift_x, helix_shift_y, last_unit_flag, Rc, Sr, distance_thresh, sw_unit_strut_flag);
    helix_shift_x = helix_shift_x + AAX;
    helix_shift_y = helix_shift_y - AAY;
    helix = [helix; unit_strut];
    if i == 1
        helix_sw_name = unit_strut_sw_name;
        current_unit_strut_name = helix_sw_name;
        main_helix_connection_point = unit_strut_connection_point;
    elseif i == Nn
        helix_sw_name = py.sld_interface.row_combine(helix_sw_name, current_unit_strut_name);% I created this function to combine rows in IR1 and reusing it. That is why the name of the function is row_combine. It just combines solid bodies
        helix_sw_name = py.sld_interface.row_combine(helix_sw_name, unit_strut_sw_name);% I created this function to combine rows in IR1 and reusing it. That is why the name of the function is row_combine. It just combines solid bodies
    else
        new_unit_strut_name = py.sld_interface.copy_move_vro(current_unit_strut_name, "Axis1", rel_z_shift/1000,rel_theta_shift);
        if i == 2
            current_unit_strut_name = new_unit_strut_name;
            continue
        else
            helix_sw_name = py.sld_interface.row_combine(helix_sw_name, current_unit_strut_name);% I created this function to combine rows in IR1 and reusing it. That is why the name of the function is row_combine. It just combines solid bodies
            current_unit_strut_name = new_unit_strut_name;
        end
    end

end
helix_main = helix;
end


function [row1, row1_rev, row1_sw_name, row1_rev_sw_name] = generate_row1_unit_strut(ts_full_strut, short_strut, crown_top, crown_bottom, crown_top_w, crown_bottom_w, cr_yscale, unit_feature_pos, unit_config, helix_shift_x, helix_shift_y, Nn, AAX, AAY, TSy, SSy, TSx, Rx, cr_yscales, alpha, beta, TSy_r1, SSy_r1, Rc, Nt1, Nt2, SL, Sy0, Sr, distance_thresh, row2_connection_point, helix_sw_name)
len = size(unit_config);
featurex = [];
featurey = [];
k = 1;
for j = 1:len(2)
    if unit_config(j) == 31
        strut_x = crown_bottom(:,2) + unit_feature_pos(1,j) + helix_shift_x;
        strut_y = crown_bottom(:,1) + unit_feature_pos(2,j) + helix_shift_y;
        featurex = [featurex; strut_x];
        featurey = [featurey; strut_y];
        strut = [strut_x strut_y];
        strut = flipud(strut);
        prev_strut_end = strut_common_point;
        strut = [prev_strut_end;strut];
        strut_close = removing_close_neibours(strut, distance_thresh);
        strut_common_point = strut_close(end - 1,:);
        sw_data_points = spline_data_points_3d(strut_close, Rc);
        spline_name = py.sld_interface.create_spline_3d(sw_data_points, prev_spline_name);
        prev_spline_name = spline_name;
        new_sweep = py.sld_interface.create_circular_sweep_profile(sw_data_points(1:6), spline_name, Sr);
        unit_strut_sw = py.sld_interface.row_combine(unit_strut_sw, new_sweep);% I created this function to combine rows in IR1 and reusing it. That is why the name of the function is row_combine. It just combines solid bodies        

    elseif unit_config(j) == 41
        strut_x = crown_bottom_w(:,2) + unit_feature_pos(1,j) + helix_shift_x ;
        strut_y = crown_bottom_w(:,1) + unit_feature_pos(2,j) + helix_shift_y;
        featurex = [featurex; strut_x];
        featurey = [featurey; strut_y];
        strut = [strut_x strut_y];
        strut = flipud(strut);

        if j ~= 1
            prev_strut_end = strut_common_point;
            strut = [prev_strut_end;strut];
        else 
            prev_spline_name = "START";
        end

        strut_close = removing_close_neibours(strut, distance_thresh);
        strut_common_point = strut_close(end - 1,:);
        sw_data_points = spline_data_points_3d(strut_close, Rc);
        spline_name = py.sld_interface.create_spline_3d(sw_data_points, prev_spline_name);
        prev_spline_name = spline_name;
        new_sweep = py.sld_interface.create_circular_sweep_profile(sw_data_points(1:6), spline_name, Sr);

        if j == 1
            unit_strut_sw = new_sweep;
            continue
        end

        unit_strut_sw = py.sld_interface.row_combine(unit_strut_sw, new_sweep);% I created this function to combine rows in IR1 and reusing it. That is why the name of the function is row_combine. It just combines solid bodies

    elseif unit_config(j) == 33
        strut_x = crown_top(:,2) + unit_feature_pos(1,j) + helix_shift_x;
        strut_y = crown_top(:,1) + unit_feature_pos(2,j) + helix_shift_y;
        featurex = [featurex; strut_x];
        featurey = [featurey; strut_y];
        strut = [strut_x strut_y];
        prev_strut_end = strut_common_point;
        strut = [prev_strut_end;strut];
        strut_close = removing_close_neibours(strut, distance_thresh);
        strut_common_point = strut_close(end - 1,:);
        sw_data_points = spline_data_points_3d(strut_close, Rc);
        spline_name = py.sld_interface.create_spline_3d(sw_data_points, prev_spline_name);
        prev_spline_name = spline_name;
        new_sweep = py.sld_interface.create_circular_sweep_profile(sw_data_points(1:6), spline_name, Sr);
        unit_strut_sw = py.sld_interface.row_combine(unit_strut_sw, new_sweep);% I created this function to combine rows in IR1 and reusing it. That is why the name of the function is row_combine. It just combines solid bodies        

    elseif unit_config(j) == 43
        strut_x = crown_top_w(:,2) + unit_feature_pos(1,j) + helix_shift_x;
        strut_y = crown_top_w(:,1) + unit_feature_pos(2,j) + helix_shift_y;
        featurex = [featurex; strut_x];
        featurey = [featurey; strut_y];
        strut = [strut_x strut_y];
        prev_strut_end = strut_common_point;
        strut = [prev_strut_end;strut];
        strut_close = removing_close_neibours(strut, distance_thresh);
        strut_common_point = strut_close(end - 1,:);
        sw_data_points = spline_data_points_3d(strut_close, Rc);
        spline_name = py.sld_interface.create_spline_3d(sw_data_points, prev_spline_name);
        prev_spline_name = spline_name;
        new_sweep = py.sld_interface.create_circular_sweep_profile(sw_data_points(1:6), spline_name, Sr);
        unit_strut_sw = py.sld_interface.row_combine(unit_strut_sw, new_sweep);% I created this function to combine rows in IR1 and reusing it. That is why the name of the function is row_combine. It just combines solid bodies        

    elseif unit_config(j) == 32
        strut_x = short_strut(:,1) + unit_feature_pos(1,j) + helix_shift_x;
        strut_y = short_strut(:,2).*(SSy_r1(k)/SSy) + unit_feature_pos(2,j) + helix_shift_y;
        
        if j == len(2) &&  ~isempty(row2_connection_point)               
            strut_x = [strut_x;row2_connection_point(1,1)];
            strut_y = [strut_y;row2_connection_point(1,2)];
        end 

        featurex = [featurex; strut_x];
        featurey = [featurey; strut_y];
        strut = [strut_x strut_y];
        prev_strut_end = strut_common_point;
        strut = [prev_strut_end;strut];
        strut_close = removing_close_neibours(strut, distance_thresh);
        strut_common_point = strut_close(end - 1,:);
        sw_data_points = spline_data_points_3d(strut_close, Rc);
        spline_name = py.sld_interface.create_spline_3d(sw_data_points, prev_spline_name);
        prev_spline_name = spline_name;
        new_sweep = py.sld_interface.create_circular_sweep_profile(sw_data_points(1:6), spline_name, Sr);
        unit_strut_sw = py.sld_interface.row_combine(unit_strut_sw, new_sweep);% I created this function to combine rows in IR1 and reusing it. That is why the name of the function is row_combine. It just combines solid bodies        
        
        if j == len(2)
            reverse_ref_pt_1 = sw_data_points(end-2:end);
            ref_X_pt_1 = strut_close(end,1);
        end

    elseif unit_config(j) == 34
        strut_x = ts_full_strut(:,1) + unit_feature_pos(1,j) + helix_shift_x ;
        strut_y = ts_full_strut(:,2).*(TSy_r1(k)/TSy) + unit_feature_pos(2,j) + helix_shift_y;
        featurex = [featurex; strut_x];
        featurey = [featurey; strut_y];
        k = k + 1;
        strut = [strut_x strut_y];
        strut = flipud(strut);
        prev_strut_end = strut_common_point;
        strut = [prev_strut_end;strut];
        strut_close = removing_close_neibours(strut, distance_thresh);
        strut_common_point = strut_close(end - 1,:);
        sw_data_points = spline_data_points_3d(strut_close, Rc);
        spline_name = py.sld_interface.create_spline_3d(sw_data_points, prev_spline_name);
        prev_spline_name = spline_name;
        new_sweep = py.sld_interface.create_circular_sweep_profile(sw_data_points(1:6), spline_name, Sr);
        unit_strut_sw = py.sld_interface.row_combine(unit_strut_sw, new_sweep);% I created this function to combine rows in IR1 and reusing it. That is why the name of the function is row_combine. It just combines solid bodies        
    end
end

unit_strut = [featurex featurey];
row1 = unit_strut;
row1_sw_name = unit_strut_sw;

% Creating method to generate the reverse struts
z_axis_name = "Axis1";
strut_end_reverse_x = -1*(strut_close(end,1) - helix_shift_x) + 2*pi*Rc - Rx + 2*pi()*Rc + (Nt1 + Nt2)*(4*Rx - TSx) + Rx - TSx + Rx + 2*AAX + Nn*AAX - Rx - (4*Rx - TSx) + 2*AAX + (Nt1 + Nt2)*(4*Rx - TSx) + Rx + Rx;
strut_end_reverse_y = -1*(strut_close(end,2) - helix_shift_y) - SL + 2*cr_yscales*Rx + Sy0;
theta_shift = (strut_end_reverse_x - strut_close(end,1))/Rc;
reverse_ref_pt_2 = spline_data_points_3d([strut_end_reverse_x strut_end_reverse_y], Rc);
p1_name = py.sld_interface.generate_ref_pts(reverse_ref_pt_1);
p2_name = py.sld_interface.generate_ref_pts([0 0 reverse_ref_pt_1(end)]);
rot_axis_name = py.sld_interface.generate_axis_between_pts(p1_name, p2_name);
reverse_strut_sw_name =  py.sld_interface.generate_reverse_pattern(row1_sw_name, rot_axis_name, z_axis_name, reverse_ref_pt_2(end) - reverse_ref_pt_1(end), theta_shift);
p3_name = py.sld_interface.generate_ref_pts(reverse_ref_pt_2);
featurex_rev = [];
featurey_rev = [];
k = 1;
for j = 1:len(2)
    if unit_config(j) == 31
        strut_x = crown_bottom(:,2) + unit_feature_pos(1,j) + helix_shift_x;
        strut_y = crown_bottom(:,1) + unit_feature_pos(2,j) + helix_shift_y;

        % Generating the reverse strut
        strut_x_rev = -1*(strut_x - helix_shift_x) + 2*pi*Rc - Rx + 2*pi()*Rc + (Nt1 + Nt2)*(4*Rx - TSx) + Rx - TSx + Rx + 2*AAX + Nn*AAX - Rx - (4*Rx - TSx) + 2*AAX + (Nt1 + Nt2)*(4*Rx - TSx) + Rx + Rx;
        strut_y_rev = -1*(strut_y - helix_shift_y) - SL + 2*cr_yscales*Rx + Sy0;
        featurex_rev = [featurex_rev; strut_x_rev];
        featurey_rev = [featurey_rev; strut_y_rev];
        
    elseif unit_config(j) == 41
        strut_x = crown_bottom_w(:,2) + unit_feature_pos(1,j) + helix_shift_x ;
        strut_y = crown_bottom_w(:,1) + unit_feature_pos(2,j) + helix_shift_y;

        % Generating the reverse strut
        strut_x_rev = -1*(strut_x - helix_shift_x) + 2*pi*Rc - Rx + 2*pi()*Rc + (Nt1 + Nt2)*(4*Rx - TSx) + Rx - TSx + Rx + 2*AAX + Nn*AAX - Rx - (4*Rx - TSx) + 2*AAX + (Nt1 + Nt2)*(4*Rx - TSx) + Rx + Rx;
        strut_y_rev = -1*(strut_y - helix_shift_y) - SL + 2*cr_yscales*Rx + Sy0;
        featurex_rev = [featurex_rev; strut_x_rev];
        featurey_rev = [featurey_rev; strut_y_rev];

    elseif unit_config(j) == 33
        strut_x = crown_top(:,2) + unit_feature_pos(1,j) + helix_shift_x;
        strut_y = crown_top(:,1) + unit_feature_pos(2,j) + helix_shift_y;

        % Generating the reverse strut
        strut_x_rev = -1*(strut_x - helix_shift_x) + 2*pi*Rc - Rx + 2*pi()*Rc + (Nt1 + Nt2)*(4*Rx - TSx) + Rx - TSx + Rx + 2*AAX + Nn*AAX - Rx - (4*Rx - TSx) + 2*AAX + (Nt1 + Nt2)*(4*Rx - TSx) + Rx + Rx;
        strut_y_rev = -1*(strut_y - helix_shift_y) - SL + 2*cr_yscales*Rx + Sy0;
        featurex_rev = [featurex_rev; strut_x_rev];
        featurey_rev = [featurey_rev; strut_y_rev];

    elseif unit_config(j) == 43
        strut_x = crown_top_w(:,2) + unit_feature_pos(1,j) + helix_shift_x;
        strut_y = crown_top_w(:,1) + unit_feature_pos(2,j) + helix_shift_y;

        % Generating the reverse strut
        strut_x_rev = -1*(strut_x - helix_shift_x) + 2*pi*Rc - Rx + 2*pi()*Rc + (Nt1 + Nt2)*(4*Rx - TSx) + Rx - TSx + Rx + 2*AAX + Nn*AAX - Rx - (4*Rx - TSx) + 2*AAX + (Nt1 + Nt2)*(4*Rx - TSx) + Rx + Rx;
        strut_y_rev = -1*(strut_y - helix_shift_y) - SL + 2*cr_yscales*Rx + Sy0;
        featurex_rev = [featurex_rev; strut_x_rev];
        featurey_rev = [featurey_rev; strut_y_rev];

    elseif unit_config(j) == 32
        strut_x = short_strut(:,1) + unit_feature_pos(1,j) + helix_shift_x;
        strut_y = short_strut(:,2).*(SSy_r1(k)/SSy) + unit_feature_pos(2,j) + helix_shift_y;
        if j == len(2) && ~isempty(row2_connection_point)
            strut_x = [strut_x;row2_connection_point(1,1)];
            strut_y = [strut_y;row2_connection_point(1,2)];
        end

        % Generating the reverse strut
        strut_x_rev = -1*(strut_x - helix_shift_x) + 2*pi*Rc - Rx + 2*pi()*Rc + (Nt1 + Nt2)*(4*Rx - TSx) + Rx - TSx + Rx + 2*AAX + Nn*AAX - Rx - (4*Rx - TSx) + 2*AAX + (Nt1 + Nt2)*(4*Rx - TSx) + Rx + Rx;
        strut_y_rev = -1*(strut_y - helix_shift_y) - SL + 2*cr_yscales*Rx + Sy0;
        featurex_rev = [featurex_rev; strut_x_rev];
        featurey_rev = [featurey_rev; strut_y_rev];    

    elseif unit_config(j) == 34
        strut_x = ts_full_strut(:,1) + unit_feature_pos(1,j) + helix_shift_x ;
        strut_y = ts_full_strut(:,2).*(TSy_r1(k)/TSy) + unit_feature_pos(2,j) + helix_shift_y;
        k = k + 1;

        % Generating the reverse strut
        strut_x_rev = -1*(strut_x - helix_shift_x) + 2*pi*Rc - Rx + 2*pi()*Rc + (Nt1 + Nt2)*(4*Rx - TSx) + Rx - TSx + Rx + 2*AAX + Nn*AAX - Rx - (4*Rx - TSx) + 2*AAX + (Nt1 + Nt2)*(4*Rx - TSx) + Rx + Rx;
        strut_y_rev = -1*(strut_y - helix_shift_y) - SL + 2*cr_yscales*Rx + Sy0;
        featurex_rev = [featurex_rev; strut_x_rev];
        featurey_rev = [featurey_rev; strut_y_rev];
    end
end

unit_strut = [featurex_rev featurey_rev];
row1_rev = unit_strut;
row1_rev_sw_name = reverse_strut_sw_name;
end

function [row2, row2_rev, row2_sw_name, row2_rev_sw_name, row2_connection_point] = generate_row2_unit_strut(ts_full_strut, short_strut, crown_top, crown_bottom, crown_top_w, crown_bottom_w, cr_yscale, unit_feature_pos, unit_config, helix_shift_x, helix_shift_y, Nn, AAX, AAY, TSy, SSy, TSx, Rx, cr_yscales, alpha, beta, TSy_r2, SSy_r2, Nt1, Nt2, Rc, SL, Sy0, Sr, distance_thresh, unit_transition_connection_point, helix_sw_name)
len = size(unit_config);
featurex = [];
featurey = [];
featurex_rev = [];
featurey_rev = [];
k = 1;
for j = 1:len(2)
    if unit_config(j) == 31
        strut_x = crown_bottom(:,2) + unit_feature_pos(1,j) + helix_shift_x;
        strut_y = crown_bottom(:,1) + unit_feature_pos(2,j) + helix_shift_y;
        featurex = [featurex; strut_x];
        featurey = [featurey; strut_y];
        strut = [strut_x strut_y];
        strut = flipud(strut);
        prev_strut_end = strut_common_point;
        strut = [prev_strut_end;strut];
        strut_close = removing_close_neibours(strut, distance_thresh);
        strut_common_point = strut_close(end - 1,:);
        sw_data_points = spline_data_points_3d(strut_close, Rc);
        spline_name = py.sld_interface.create_spline_3d(sw_data_points, prev_spline_name);
        prev_spline_name = spline_name;
        new_sweep = py.sld_interface.create_circular_sweep_profile(sw_data_points(1:6), spline_name, Sr);
        unit_strut_sw = py.sld_interface.row_combine(unit_strut_sw, new_sweep);% I created this function to combine rows in IR1 and reusing it. That is why the name of the function is row_combine. It just combines solid bodies        

    elseif unit_config(j) == 41
        strut_x = crown_bottom_w(:,2) + unit_feature_pos(1,j) + helix_shift_x;
        strut_y = crown_bottom_w(:,1) + unit_feature_pos(2,j) + helix_shift_y;
        featurex = [featurex; strut_x];
        featurey = [featurey; strut_y];
        strut = [strut_x strut_y];
        strut = flipud(strut);
        prev_strut_end = strut_common_point;
        strut = [prev_strut_end;strut];
        strut_close = removing_close_neibours(strut, distance_thresh);
        strut_common_point = strut_close(end - 1,:);
        sw_data_points = spline_data_points_3d(strut_close, Rc);
        spline_name = py.sld_interface.create_spline_3d(sw_data_points, prev_spline_name);
        prev_spline_name = spline_name;
        new_sweep = py.sld_interface.create_circular_sweep_profile(sw_data_points(1:6), spline_name, Sr);
        unit_strut_sw = py.sld_interface.row_combine(unit_strut_sw, new_sweep);% I created this function to combine rows in IR1 and reusing it. That is why the name of the function is row_combine. It just combines solid bodies        

    elseif unit_config(j) == 33
        strut_x = crown_top(:,2) + unit_feature_pos(1,j) + helix_shift_x;
        strut_y = crown_top(:,1) + unit_feature_pos(2,j) + helix_shift_y;
        featurex = [featurex; strut_x];
        featurey = [featurey; strut_y];
        strut = [strut_x strut_y];
        prev_strut_end = strut_common_point;
        strut = [prev_strut_end;strut];
        strut_close = removing_close_neibours(strut, distance_thresh);
        strut_common_point = strut_close(end - 1,:);
        sw_data_points = spline_data_points_3d(strut_close, Rc);
        spline_name = py.sld_interface.create_spline_3d(sw_data_points, prev_spline_name);
        prev_spline_name = spline_name;
        new_sweep = py.sld_interface.create_circular_sweep_profile(sw_data_points(1:6), spline_name, Sr);
        unit_strut_sw = py.sld_interface.row_combine(unit_strut_sw, new_sweep);% I created this function to combine rows in IR1 and reusing it. That is why the name of the function is row_combine. It just combines solid bodies        

    elseif unit_config(j) == 43
        strut_x = crown_top_w(:,2) + unit_feature_pos(1,j) + helix_shift_x;
        strut_y = crown_top_w(:,1) + unit_feature_pos(2,j) + helix_shift_y;
        featurex = [featurex; strut_x];
        featurey = [featurey; strut_y];
        strut = [strut_x strut_y];

        if j ~= 1
            prev_strut_end = strut_common_point;
            strut = [prev_strut_end;strut];
        
        else 
            prev_spline_name = "START";
        end

        strut_close = removing_close_neibours(strut, distance_thresh);
        strut_common_point = strut_close(end - 1,:);
        sw_data_points = spline_data_points_3d(strut_close, Rc);
        spline_name = py.sld_interface.create_spline_3d(sw_data_points, prev_spline_name);
        prev_spline_name = spline_name;
        new_sweep = py.sld_interface.create_circular_sweep_profile(sw_data_points(1:6), spline_name, Sr);
        
        if j == 1
            unit_strut_sw = new_sweep;
            row2_connection_point = strut_close(4,:);
            continue
        end

        unit_strut_sw = py.sld_interface.row_combine(unit_strut_sw, new_sweep);% I created this function to combine rows in IR1 and reusing it. That is why the name of the function is row_combine. It just combines solid bodies        

    elseif unit_config(j) == 32
        strut_x = short_strut(:,1) + unit_feature_pos(1,j) + helix_shift_x;
        strut_y = short_strut(:,2).*(SSy_r2(k)/SSy) + unit_feature_pos(2,j) + helix_shift_y;
        featurex = [featurex; strut_x];
        featurey = [featurey; strut_y];
        k = k + 1;
        strut = [strut_x strut_y];
        prev_strut_end = strut_common_point;
        strut = [prev_strut_end;strut];
        strut_close = removing_close_neibours(strut, distance_thresh);
        strut_common_point = strut_close(end - 1,:);
        sw_data_points = spline_data_points_3d(strut_close, Rc);
        spline_name = py.sld_interface.create_spline_3d(sw_data_points, prev_spline_name);
        prev_spline_name = spline_name;
        new_sweep = py.sld_interface.create_circular_sweep_profile(sw_data_points(1:6), spline_name, Sr);
        unit_strut_sw = py.sld_interface.row_combine(unit_strut_sw, new_sweep);% I created this function to combine rows in IR1 and reusing it. That is why the name of the function is row_combine. It just combines solid bodies        

    elseif unit_config(j) == 34
        strut_x = ts_full_strut(:,1) + unit_feature_pos(1,j) + helix_shift_x ;
        strut_y = ts_full_strut(:,2).*(TSy_r2(k)/TSy) + unit_feature_pos(2,j) + helix_shift_y;
        featurex = [featurex; strut_x];
        featurey = [featurey; strut_y];

        if j == len(2) && ~isempty(unit_transition_connection_point)
            strut_x = [unit_transition_connection_point(1,1); strut_x];
            strut_y = [unit_transition_connection_point(1,2); strut_y];
        end

        strut = [strut_x strut_y];
        strut = flipud(strut);
        prev_strut_end = strut_common_point;
        strut = [prev_strut_end;strut];
        strut_close = removing_close_neibours(strut, distance_thresh);
        strut_common_point = strut_close(end - 1,:);
        sw_data_points = spline_data_points_3d(strut_close, Rc);
        spline_name = py.sld_interface.create_spline_3d(sw_data_points, prev_spline_name);
        prev_spline_name = spline_name;
        new_sweep = py.sld_interface.create_circular_sweep_profile(sw_data_points(1:6), spline_name, Sr);
        unit_strut_sw = py.sld_interface.row_combine(unit_strut_sw, new_sweep);% I created this function to combine rows in IR1 and reusing it. That is why the name of the function is row_combine. It just combines solid bodies        
  
        if j == len(2)
            reverse_ref_pt_1 = sw_data_points(end-2:end);
            ref_X_pt_1 = strut_close(end,1);
        end
    end
end

k = 1;
unit_strut = [featurex featurey];
row2 = unit_strut;
row2_sw_name = unit_strut_sw;

% Creating method to generate the reverse struts
z_axis_name = "Axis1";
strut_end_reverse_x = -1*(strut_close(end,1) - helix_shift_x) + (Nt1 + Nt2)*(4*Rx-TSx) + Rx + 2*pi()*Rc + (Nt1 + Nt2)*(4*Rx - TSx) + Rx - TSx + Rx + 2*AAX + Nn*AAX - Rx - (4*Rx - TSx) + 2*AAX;
strut_end_reverse_y = -1*(strut_close(end,2) - helix_shift_y) -SL + Sy0 + 2*cr_yscales*Rx;
theta_shift = (strut_end_reverse_x - strut_close(end,1))/Rc;
reverse_ref_pt_2 = spline_data_points_3d([strut_end_reverse_x strut_end_reverse_y], Rc);
p1_name = py.sld_interface.generate_ref_pts(reverse_ref_pt_1);
p2_name = py.sld_interface.generate_ref_pts([0 0 reverse_ref_pt_1(end)]);
rot_axis_name = py.sld_interface.generate_axis_between_pts(p1_name, p2_name);
reverse_strut_sw_name =  py.sld_interface.generate_reverse_pattern(row2_sw_name, rot_axis_name, z_axis_name, reverse_ref_pt_2(end) - reverse_ref_pt_1(end), theta_shift);
p3_name = py.sld_interface.generate_ref_pts(reverse_ref_pt_2);

for j = 1:len(2)
    if unit_config(j) == 31
        strut_x = crown_bottom(:,2) + unit_feature_pos(1,j) + helix_shift_x;
        strut_y = crown_bottom(:,1) + unit_feature_pos(2,j) + helix_shift_y;

        % Generating the reverse strut
        strut_x_rev = -1*(strut_x - helix_shift_x) + (Nt1 + Nt2)*(4*Rx-TSx) + Rx + 2*pi()*Rc + (Nt1 + Nt2)*(4*Rx - TSx) + Rx - TSx + Rx + 2*AAX + Nn*AAX - Rx - (4*Rx - TSx) + 2*AAX;
        strut_y_rev = -1*(strut_y - helix_shift_y) -SL + Sy0 + 2*cr_yscales*Rx;
        featurex_rev = [featurex_rev; strut_x_rev];
        featurey_rev = [featurey_rev; strut_x_rev];         

    elseif unit_config(j) == 41
        strut_x = crown_bottom_w(:,2) + unit_feature_pos(1,j) + helix_shift_x;
        strut_y = crown_bottom_w(:,1) + unit_feature_pos(2,j) + helix_shift_y;

        % Generating the reverse strut
        strut_x_rev = -1*(strut_x - helix_shift_x) + (Nt1 + Nt2)*(4*Rx-TSx) + Rx + 2*pi()*Rc + (Nt1 + Nt2)*(4*Rx - TSx) + Rx - TSx + Rx + 2*AAX + Nn*AAX - Rx - (4*Rx - TSx) + 2*AAX;
        strut_y_rev = -1*(strut_y - helix_shift_y) -SL + Sy0 + 2*cr_yscales*Rx;
        featurex_rev = [featurex_rev; strut_x_rev];
        featurey_rev = [featurey_rev; strut_x_rev]; 

    elseif unit_config(j) == 33
        strut_x = crown_top(:,2) + unit_feature_pos(1,j) + helix_shift_x;
        strut_y = crown_top(:,1) + unit_feature_pos(2,j) + helix_shift_y;

        % Generating the reverse strut
        strut_x_rev = -1*(strut_x - helix_shift_x) + (Nt1 + Nt2)*(4*Rx-TSx) + Rx + 2*pi()*Rc + (Nt1 + Nt2)*(4*Rx - TSx) + Rx - TSx + Rx + 2*AAX + Nn*AAX - Rx - (4*Rx - TSx) + 2*AAX;
        strut_y_rev = -1*(strut_y - helix_shift_y) -SL + Sy0 + 2*cr_yscales*Rx;
        featurex_rev = [featurex_rev; strut_x_rev];
        featurey_rev = [featurey_rev; strut_x_rev]; 

    elseif unit_config(j) == 43
        strut_x = crown_top_w(:,2) + unit_feature_pos(1,j) + helix_shift_x;
        strut_y = crown_top_w(:,1) + unit_feature_pos(2,j) + helix_shift_y;

        % Generating the reverse strut
        strut_x_rev = -1*(strut_x - helix_shift_x) + (Nt1 + Nt2)*(4*Rx-TSx) + Rx + 2*pi()*Rc + (Nt1 + Nt2)*(4*Rx - TSx) + Rx - TSx + Rx + 2*AAX + Nn*AAX - Rx - (4*Rx - TSx) + 2*AAX;
        strut_y_rev = -1*(strut_y - helix_shift_y) -SL + Sy0 + 2*cr_yscales*Rx;
        featurex_rev = [featurex_rev; strut_x_rev];
        featurey_rev = [featurey_rev; strut_y_rev]; 

    elseif unit_config(j) == 32
        strut_x = short_strut(:,1) + unit_feature_pos(1,j) + helix_shift_x;
        strut_y = short_strut(:,2).*(SSy_r2(k)/SSy) + unit_feature_pos(2,j) + helix_shift_y;
        k = k + 1;

        % Generating the reverse strut
        strut_x_rev = -1*(strut_x - helix_shift_x) + (Nt1 + Nt2)*(4*Rx-TSx) + Rx + 2*pi()*Rc + (Nt1 + Nt2)*(4*Rx - TSx) + Rx - TSx + Rx + 2*AAX + Nn*AAX - Rx - (4*Rx - TSx) + 2*AAX;
        strut_y_rev = -1*(strut_y - helix_shift_y) -SL + Sy0 + 2*cr_yscales*Rx;
        featurex_rev = [featurex_rev; strut_x_rev];
        featurey_rev = [featurey_rev; strut_y_rev]; 

    elseif unit_config(j) == 34
        strut_x = ts_full_strut(:,1) + unit_feature_pos(1,j) + helix_shift_x ;
        strut_y = ts_full_strut(:,2).*(TSy_r2(k)/TSy) + unit_feature_pos(2,j) + helix_shift_y;

% Generating the reverse strut
        strut_x_rev = -1*(strut_x - helix_shift_x) + (Nt1 + Nt2)*(4*Rx-TSx) + Rx + 2*pi()*Rc + (Nt1 + Nt2)*(4*Rx - TSx) + Rx - TSx + Rx + 2*AAX + Nn*AAX - Rx - (4*Rx - TSx) + 2*AAX;
        strut_y_rev = -1*(strut_y - helix_shift_y) -SL + Sy0 + 2*cr_yscales*Rx;

        featurex_rev = [featurex_rev; strut_x_rev];
        featurey_rev = [featurey_rev; strut_y_rev]; 
        if j == len(2) && ~isempty(unit_transition_connection_point)
            strut_x = [unit_transition_connection_point(1,1); strut_x];
            strut_y = [unit_transition_connection_point(1,2); strut_y];
        end
    end
end
unit_strut = [featurex_rev featurey_rev];
row2_rev = unit_strut;
row2_rev_sw_name = reverse_strut_sw_name;
end

function [transition_strut, transition_strut_rev, transition_strut_sw_name, transition_strut_rev_sw_name, unit_transition_connection_point] = generate_transition_unit_strut(ts_corrected_strut, short_strut, crown_top, crown_bottom, crown_top_w, crown_bottom_w, cr_yscale, unit_feature_pos, unit_config, helix_shift_x, helix_shift_y, Nn, AAX, AAY, TSy, SSy, TSx, Rx, cr_yscales, alpha, SL, Sy0, beta, Nt1, Nt2, Rc, Sr, distance_thresh, main_helix_connection_point, helix_sw_name)
len = size(unit_config);
featurex = [];
featurey = [];
featurex_rev = [];
featurey_rev = [];

for j = 1:len(2)
    if unit_config(j) == 1
        strut_x = crown_bottom(:,2) + unit_feature_pos(1,j) + helix_shift_x;
        strut_y = crown_bottom(:,1) + unit_feature_pos(2,j) + helix_shift_y;
        featurex = [featurex; strut_x];
        featurey = [featurey; strut_y];
        strut = [strut_x strut_y];
        strut = flipud(strut);
        prev_strut_end = strut_common_point;
        strut = [prev_strut_end;strut];
        strut_close = removing_close_neibours(strut, distance_thresh);
        strut_common_point = strut_close(end - 1,:);
        sw_data_points = spline_data_points_3d(strut_close, Rc);
        spline_name = py.sld_interface.create_spline_3d(sw_data_points, prev_spline_name);
        prev_spline_name = spline_name;
        new_sweep = py.sld_interface.create_circular_sweep_profile(sw_data_points(1:6), spline_name, Sr);
        unit_strut_sw = py.sld_interface.row_combine(unit_strut_sw, new_sweep);% I created this function to combine rows in IR1 and reusing it. That is why the name of the function is row_combine. It just combines solid bodies        

    elseif unit_config(j) == 11
        strut_x = crown_bottom_w(:,2) + unit_feature_pos(1,j) + helix_shift_x;
        strut_y = crown_bottom_w(:,1) + unit_feature_pos(2,j) + helix_shift_y;
        featurex = [featurex; strut_x];
        featurey = [featurey; strut_y];
        strut = [strut_x strut_y];
        strut = flipud(strut);
        if j ~= 1
            prev_strut_end = strut_common_point;
            strut = [prev_strut_end;strut];
        else 
            prev_spline_name = "START";
        end
        strut_close = removing_close_neibours(strut, distance_thresh);
        strut_common_point = strut_close(end - 1,:);
        sw_data_points = spline_data_points_3d(strut_close, Rc);
        spline_name = py.sld_interface.create_spline_3d(sw_data_points, prev_spline_name);
        prev_spline_name = spline_name;
        new_sweep = py.sld_interface.create_circular_sweep_profile(sw_data_points(1:6), spline_name, Sr);
        if j == 1
            unit_strut_sw = new_sweep;
            unit_transition_connection_point = strut_close(1,:);
            continue
        end
        unit_strut_sw = py.sld_interface.row_combine(unit_strut_sw, new_sweep);% I created this function to combine rows in IR1 and reusing it. That is why the name of the function is row_combine. It just combines solid bodies

    elseif unit_config(j) == 3
        strut_x = crown_top(:,2) + unit_feature_pos(1,j) + helix_shift_x;
        strut_y = crown_top(:,1) + unit_feature_pos(2,j) + helix_shift_y;
        featurex = [featurex; strut_x];
        featurey = [featurey; strut_y];
        strut = [strut_x strut_y];
        prev_strut_end = strut_common_point;
        strut = [prev_strut_end;strut];
        strut_close = removing_close_neibours(strut, distance_thresh);
        strut_common_point = strut_close(end - 1,:);
        sw_data_points = spline_data_points_3d(strut_close, Rc);
        spline_name = py.sld_interface.create_spline_3d(sw_data_points, prev_spline_name);
        prev_spline_name = spline_name;
        new_sweep = py.sld_interface.create_circular_sweep_profile(sw_data_points(1:6), spline_name, Sr);
        unit_strut_sw = py.sld_interface.row_combine(unit_strut_sw, new_sweep);% I created this function to combine rows in IR1 and reusing it. That is why the name of the function is row_combine. It just combines solid bodies        

    elseif unit_config(j) == 13
        strut_x = crown_top_w(:,2) + unit_feature_pos(1,j) + helix_shift_x;
        strut_y = crown_top_w(:,1) + unit_feature_pos(2,j) + helix_shift_y;
        featurex = [featurex; strut_x];
        featurey = [featurey; strut_y];
        strut = [strut_x strut_y];
        prev_strut_end = strut_common_point;
        strut = [prev_strut_end;strut];
        strut_close = removing_close_neibours(strut, distance_thresh);
        strut_common_point = strut_close(end - 1,:);
        sw_data_points = spline_data_points_3d(strut_close, Rc);
        spline_name = py.sld_interface.create_spline_3d(sw_data_points, prev_spline_name);
        prev_spline_name = spline_name;
        new_sweep = py.sld_interface.create_circular_sweep_profile(sw_data_points(1:6), spline_name, Sr);
        unit_strut_sw = py.sld_interface.row_combine(unit_strut_sw, new_sweep);% I created this function to combine rows in IR1 and reusing it. That is why the name of the function is row_combine. It just combines solid bodies             

    elseif unit_config(j) == 2
        strut_x = short_strut(:,1) + unit_feature_pos(1,j) + helix_shift_x;
        strut_y = short_strut(:,2) + unit_feature_pos(2,j) + helix_shift_y;
        featurex = [featurex; strut_x];
        featurey = [featurey; strut_y];
        strut = [strut_x strut_y];
        prev_strut_end = strut_common_point;
        strut = [prev_strut_end;strut];
        strut_close = removing_close_neibours(strut, distance_thresh);
        strut_common_point = strut_close(end - 1,:);
        sw_data_points = spline_data_points_3d(strut_close, Rc);
        spline_name = py.sld_interface.create_spline_3d(sw_data_points, prev_spline_name);
        prev_spline_name = spline_name;
        new_sweep = py.sld_interface.create_circular_sweep_profile(sw_data_points(1:6), spline_name, Sr);
        unit_strut_sw = py.sld_interface.row_combine(unit_strut_sw, new_sweep);% I created this function to combine rows in IR1 and reusing it. That is why the name of the function is row_combine. It just combines solid bodies        

    elseif unit_config(j) == 4
        strut_x = ts_corrected_strut(:,1) + unit_feature_pos(1,j) + helix_shift_x ;
        strut_y = ts_corrected_strut(:,2) + unit_feature_pos(2,j) + helix_shift_y;
        featurex = [featurex; strut_x];
        featurey = [featurey; strut_y];
        strut = [strut_x strut_y];
        strut = flipud(strut);
        if j == len(2) && ~isempty(main_helix_connection_point)
            strut = [strut; main_helix_connection_point];
        end
        prev_strut_end = strut_common_point;
        strut = [prev_strut_end;strut];
        strut_close = removing_close_neibours(strut, distance_thresh);
        strut_common_point = strut_close(end - 1,:);
        sw_data_points = spline_data_points_3d(strut_close, Rc);
        spline_name = py.sld_interface.create_spline_3d(sw_data_points, prev_spline_name);
        prev_spline_name = spline_name;
        new_sweep = py.sld_interface.create_circular_sweep_profile(sw_data_points(1:6), spline_name, Sr);
        unit_strut_sw = py.sld_interface.row_combine(unit_strut_sw, new_sweep);% I created this function to combine rows in IR1 and reusing it. That is why the name of the function is row_combine. It just combines solid bodies        
        if j == len(2)
            reverse_ref_pt_1 = sw_data_points(end-2:end);
            ref_X_pt_1 = strut_close(end,1);
        end
    end
end
unit_strut = [featurex featurey];
transition_strut = unit_strut;
transition_strut_sw_name = unit_strut_sw;

% Creating method to generate the reverse struts
z_axis_name = "Axis1";
strut_end_reverse_x = -1*(strut_close(end,1) - helix_shift_x) + 1*helix_shift_x + 2*AAX + (Nn+1)*AAX + 4*(4*Rx - TSx) - 2*Rx + TSx;%- 4*Rx + TSx;% Check if 2 Tsx or 1 Tsx
strut_end_reverse_y = -1*(strut_close(end,2) - helix_shift_y) - SL + Sy0 + 2*Rx*cr_yscales + tand(beta)*(Nt1+Nt2)*(4*Rx-TSx) + 2*Rx*cr_yscales + TSy;
theta_shift = (strut_end_reverse_x - strut_close(end,1))/Rc;
reverse_ref_pt_2 = spline_data_points_3d([strut_end_reverse_x strut_end_reverse_y], Rc);
p1_name = py.sld_interface.generate_ref_pts(reverse_ref_pt_1);
p2_name = py.sld_interface.generate_ref_pts([0 0 reverse_ref_pt_1(end)]);
rot_axis_name = py.sld_interface.generate_axis_between_pts(p1_name, p2_name);
reverse_strut_sw_name =  py.sld_interface.generate_reverse_pattern(transition_strut_sw_name, rot_axis_name, z_axis_name, reverse_ref_pt_2(end) - reverse_ref_pt_1(end), theta_shift);
p3_name = py.sld_interface.generate_ref_pts(reverse_ref_pt_2);

for j = 1:len(2)
    if unit_config(j) == 1
        strut_x = crown_bottom(:,2) + unit_feature_pos(1,j) + helix_shift_x;
        strut_y = crown_bottom(:,1) + unit_feature_pos(2,j) + helix_shift_y;

        % Generating the reverse strut
        strut_x_rev = -1*(strut_x - helix_shift_x) + 1*helix_shift_x + 2*AAX + (Nn+1)*AAX + 4*(4*Rx - TSx)- 2*Rx + TSx;%- 4*Rx + TSx + 2*AAX;
        strut_y_rev = -1*(strut_y - helix_shift_y) - SL + Sy0 + 2*Rx*cr_yscales + tand(beta)*(Nt1+Nt2)*(4*Rx-TSx) + 2*Rx*cr_yscales + TSy;;
        featurex_rev = [featurex_rev; strut_x_rev];
        featurey_rev = [featurey_rev; strut_y_rev];

    elseif unit_config(j) == 11
        strut_x = crown_bottom_w(:,2) + unit_feature_pos(1,j) + helix_shift_x;
        strut_y = crown_bottom_w(:,1) + unit_feature_pos(2,j) + helix_shift_y;

        % Generating the reverse strut
        strut_x_rev = -1*(strut_x - helix_shift_x) + 1*helix_shift_x + 2*AAX + (Nn+1)*AAX + 4*(4*Rx - TSx)- 2*Rx + TSx;%- 4*Rx + TSx;
        strut_y_rev = -1*(strut_y - helix_shift_y) - SL + Sy0 + 2*Rx*cr_yscales + tand(beta)*(Nt1+Nt2)*(4*Rx-TSx) + 2*Rx*cr_yscales + TSy;
        featurex_rev = [featurex_rev; strut_x_rev];
        featurey_rev = [featurey_rev; strut_y_rev];

    elseif unit_config(j) == 3
        strut_x = crown_top(:,2) + unit_feature_pos(1,j) + helix_shift_x;
        strut_y = crown_top(:,1) + unit_feature_pos(2,j) + helix_shift_y;

        % Generating the reverse strut
        strut_x_rev = -1*(strut_x - helix_shift_x) + 1*helix_shift_x + 2*AAX + (Nn+1)*AAX + 4*(4*Rx - TSx) - 2*Rx + TSx;%- 4*Rx + TSx;
        strut_y_rev = -1*(strut_y - helix_shift_y) - SL + Sy0 + 2*Rx*cr_yscales + tand(beta)*(Nt1+Nt2)*(4*Rx-TSx) + 2*Rx*cr_yscales + TSy;; 
        featurex_rev = [featurex_rev; strut_x_rev];
        featurey_rev = [featurey_rev; strut_y_rev];

    elseif unit_config(j) == 13
        strut_x = crown_top_w(:,2) + unit_feature_pos(1,j) + helix_shift_x;
        strut_y = crown_top_w(:,1) + unit_feature_pos(2,j) + helix_shift_y;

        % Generating the reverse strut
        strut_x_rev = -1*(strut_x - helix_shift_x) + 1*helix_shift_x + 2*AAX + (Nn+1)*AAX + 4*(4*Rx - TSx) - 2*Rx + TSx;%- 4*Rx + TSx;
        strut_y_rev = -1*(strut_y - helix_shift_y) - SL + Sy0 + 2*Rx*cr_yscales + tand(beta)*(Nt1+Nt2)*(4*Rx-TSx) + 2*Rx*cr_yscales + TSy;
        featurex_rev = [featurex_rev; strut_x_rev];
        featurey_rev = [featurey_rev; strut_y_rev];
        
    elseif unit_config(j) == 2
        strut_x = short_strut(:,1) + unit_feature_pos(1,j) + helix_shift_x;
        strut_y = short_strut(:,2) + unit_feature_pos(2,j) + helix_shift_y;

        % Generating the reverse strut
        strut_x_rev = -1*(strut_x - helix_shift_x) + 1*helix_shift_x + 2*AAX + (Nn+1)*AAX + 4*(4*Rx - TSx) - 2*Rx + TSx;%- 4*Rx + TSx;
        strut_y_rev = -1*(strut_y - helix_shift_y) - SL + Sy0 + 2*Rx*cr_yscales + tand(beta)*(Nt1+Nt2)*(4*Rx-TSx) + 2*Rx*cr_yscales + TSy;
        featurex_rev = [featurex_rev; strut_x_rev];
        featurey_rev = [featurey_rev; strut_y_rev];

    elseif unit_config(j) == 4
        strut_x = ts_corrected_strut(:,1) + unit_feature_pos(1,j) + helix_shift_x ;
        strut_y = ts_corrected_strut(:,2) + unit_feature_pos(2,j) + helix_shift_y;
        if j == len(2) && ~isempty(main_helix_connection_point)
            strut_x = [main_helix_connection_point(1,1);strut_x];
            strut_y = [main_helix_connection_point(1,2);strut_y];
        end

        % Generating the reverse strut
        strut_x_rev = -1*(strut_x - helix_shift_x) + 1*helix_shift_x + 2*AAX + (Nn+1)*AAX + 4*(4*Rx - TSx) - 2*Rx + TSx;%- 4*Rx + TSx;% Check if 2 Tsx or 1 Tsx
        strut_y_rev = -1*(strut_y - helix_shift_y) - SL + Sy0 + 2*Rx*cr_yscales + tand(beta)*(Nt1+Nt2)*(4*Rx-TSx) + 2*Rx*cr_yscales + TSy;
        featurex_rev = [featurex_rev; strut_x_rev];
        featurey_rev = [featurey_rev; strut_y_rev];
    end
end

unit_strut_rev = [featurex_rev featurey_rev];
transition_strut_rev = unit_strut_rev;
transition_strut_rev_sw_name =reverse_strut_sw_name;
end

function [unit_strut, unit_strut_sw_name, unit_strut_connection_point] = generate_unit_strut(ts_corrected_strut, short_strut, crown_top, crown_bottom, crown_top_w, crown_bottom_w, cr_yscale, unit_feature_pos, unit_config, helix_shift_x, helix_shift_y, last_unit_flag, Rc, Sr, distance_thresh, sw_unit_strut_flag)
len = size(unit_config);
featurex = [];
featurey = [];
unit_strut_connection_point = [];
for j = 1:len(2)
    if unit_config(j) == 1
        strut_x = crown_bottom(:,2) + unit_feature_pos(1,j) + helix_shift_x;
        strut_y = crown_bottom(:,1) + unit_feature_pos(2,j) + helix_shift_y;
        featurex = [featurex; strut_x];
        featurey = [featurey; strut_y];
        if sw_unit_strut_flag == 1
            strut = [strut_x strut_y];
            strut = flipud(strut);
            prev_strut_end = strut_common_point;
            strut = [prev_strut_end;strut];
            strut_close = removing_close_neibours(strut, distance_thresh);
            strut_common_point = strut_close(end - 1,:);
            sw_data_points = spline_data_points_3d(strut_close, Rc);
            %py.sld_interface.create_new_sketch();% Unit feature starts here so create new sketch at Front Plane
            spline_name = py.sld_interface.create_spline_3d(sw_data_points, prev_spline_name);
            %pause(1);
            prev_spline_name = spline_name;
            new_sweep = py.sld_interface.create_circular_sweep_profile(sw_data_points(1:6), spline_name, Sr);
            unit_strut_sw = py.sld_interface.row_combine(unit_strut_sw, new_sweep);% I created this function to combine rows in IR1 and reusing it. That is why the name of the function is row_combine. It just combines solid bodies
        end

    elseif unit_config(j) == 11
        strut_x = crown_bottom_w(:,2) + unit_feature_pos(1,j) + helix_shift_x;
        strut_y = crown_bottom_w(:,1) + unit_feature_pos(2,j) + helix_shift_y;
        featurex = [featurex; strut_x];
        featurey = [featurey; strut_y];
        if sw_unit_strut_flag == 1
            strut = [strut_x strut_y];
            strut = flipud(strut);
            if j ~= 1
                prev_strut_end = strut_common_point;
                strut = [prev_strut_end;strut];
            else
                prev_spline_name = "START";
            end
            strut_close = removing_close_neibours(strut, distance_thresh);
            strut_common_point = strut_close(end - 1,:);
            sw_data_points = spline_data_points_3d(strut_close, Rc);
            %py.sld_interface.create_new_sketch();% Unit feature starts here so create new sketch at Front Plane
            spline_name = py.sld_interface.create_spline_3d(sw_data_points, prev_spline_name);
            new_sweep = py.sld_interface.create_circular_sweep_profile(sw_data_points(1:6), spline_name, Sr);
            prev_spline_name = spline_name;
            if j == 1
                unit_strut_sw = new_sweep;
                unit_strut_connection_point = strut_close(2,:);
            end
        end

    elseif unit_config(j) == 3
        strut_x = crown_top(:,2) + unit_feature_pos(1,j) + helix_shift_x;
        strut_y = crown_top(:,1) + unit_feature_pos(2,j) + helix_shift_y;
        featurex = [featurex; strut_x];
        featurey = [featurey; strut_y];
        if sw_unit_strut_flag == 1
            strut = [strut_x strut_y];
            prev_strut_end = strut_common_point;
            strut = [prev_strut_end;strut];
            strut_close = removing_close_neibours(strut, distance_thresh);     
            strut_common_point = strut_close(end - 1,:);
            sw_data_points = spline_data_points_3d(strut_close, Rc);
            %py.sld_interface.create_new_sketch();% Unit feature starts here so create new sketch at Front Plane
            spline_name = py.sld_interface.create_spline_3d(sw_data_points, prev_spline_name);
            prev_spline_name = spline_name;
            new_sweep = py.sld_interface.create_circular_sweep_profile(sw_data_points(1:6), spline_name, Sr);
            unit_strut_sw = py.sld_interface.row_combine(unit_strut_sw, new_sweep);% I created this function to combine rows in IR1 and reusing it. That is why the name of the function is row_combine. It just combines solid bodies
        end

    elseif unit_config(j) == 13
        strut_x = crown_top_w(:,2) + unit_feature_pos(1,j) + helix_shift_x;
        strut_y = crown_top_w(:,1) + unit_feature_pos(2,j) + helix_shift_y;
        featurex = [featurex; strut_x];
        featurey = [featurey; strut_y];
        if sw_unit_strut_flag == 1
            strut = [strut_x strut_y];
            prev_strut_end = strut_common_point;
            strut = [prev_strut_end;strut];
            strut_close = removing_close_neibours(strut, distance_thresh); 
            strut_common_point = strut_close(end - 1,:);
            sw_data_points = spline_data_points_3d(strut_close, Rc);
            %py.sld_interface.create_new_sketch();% Unit feature starts here so create new sketch at Front Plane
            spline_name = py.sld_interface.create_spline_3d(sw_data_points, prev_spline_name);
            prev_spline_name = spline_name;
            new_sweep = py.sld_interface.create_circular_sweep_profile(sw_data_points(1:6), spline_name, Sr);
            unit_strut_sw = py.sld_interface.row_combine(unit_strut_sw, new_sweep);% I created this function to combine rows in IR1 and reusing it. That is why the name of the function is row_combine. It just combines solid bodies
        end  
        if last_unit_flag == 1
            break;
        end

    elseif unit_config(j) == 2
        strut_x = short_strut(:,1) + unit_feature_pos(1,j) + helix_shift_x;
        strut_y = short_strut(:,2) + unit_feature_pos(2,j) + helix_shift_y;
        featurex = [featurex; strut_x];
        featurey = [featurey; strut_y];
        if sw_unit_strut_flag == 1
            strut = [strut_x strut_y];
            prev_strut_end = strut_common_point;
            strut = [prev_strut_end;strut];
            strut_close = removing_close_neibours(strut, distance_thresh);
            strut_common_point = strut_close(end - 1,:);
            sw_data_points = spline_data_points_3d(strut_close, Rc);
            %py.sld_interface.create_new_sketch();% Unit feature starts here so create new sketch at Front Plane
            spline_name = py.sld_interface.create_spline_3d(sw_data_points, prev_spline_name);
            prev_spline_name = spline_name;
            new_sweep = py.sld_interface.create_circular_sweep_profile(sw_data_points(1:6), spline_name, Sr);
            unit_strut_sw = py.sld_interface.row_combine(unit_strut_sw, new_sweep);% I created this function to combine rows in IR1 and reusing it. That is why the name of the function is row_combine. It just combines solid bodies        
        end

    elseif unit_config(j) == 4
        strut_x = ts_corrected_strut(:,1) + unit_feature_pos(1,j) + helix_shift_x;
        strut_y = ts_corrected_strut(:,2) + unit_feature_pos(2,j) + helix_shift_y;
        featurex = [featurex; strut_x];
        featurey = [featurey; strut_y];
        if sw_unit_strut_flag == 1
            strut = [strut_x strut_y];
            strut = flipud(strut);
            prev_strut_end = strut_common_point;
            strut = [prev_strut_end;strut];
            strut_close = removing_close_neibours(strut, distance_thresh); 
            strut_common_point = strut_close(end - 1,:);
            sw_data_points = spline_data_points_3d(strut_close, Rc);
            %py.sld_interface.create_new_sketch();% Unit feature starts here so create new sketch at Front Plane
            spline_name = py.sld_interface.create_spline_3d(sw_data_points, prev_spline_name);
            prev_spline_name = spline_name;
            new_sweep = py.sld_interface.create_circular_sweep_profile(sw_data_points(1:6), spline_name, Sr);
            unit_strut_sw = py.sld_interface.row_combine(unit_strut_sw, new_sweep);% I created this function to combine rows in IR1 and reusing it. That is why the name of the function is row_combine. It just combines solid bodies
        end
    end
end
if sw_unit_strut_flag == 1
    unit_strut_sw_name = unit_strut_sw;
else
    unit_strut_sw_name = [];
end
unit_strut = [featurex featurey];
end

function [unit_feature_pos, TSy_r1, SSy_r1] = feature_pos_unit_top1(unit_config, Rx, cr_yscale, cr_yscales, TSx, TSy, SSy, alpha, beta, tantheta, Sy0)
ufpx = [];
ufpy = [];
len = size(unit_config);
y_shift = Rx*cr_yscales;
x_shift = -Rx;

% Deriving the lengths of each of the row 2 TSy's and SSys
SSy11 = Sy0;
SSy12 = SSy11 + (4*Rx - TSx)*tand(beta);
TSy11 = SSy12;
SSy13 = SSy12 + (4*Rx - TSx)*tand(beta);
TSy12 = SSy13;
SSy14 = SSy13 + (4*Rx - TSx)*tand(beta);
TSy13 = SSy14;
SSy15 = SSy14 + (4*Rx - TSx)*tand(beta);
TSy14 = SSy15;
SSy16 = SSy15 + (4*Rx - TSx)*tand(beta);
TSy15 = SSy16;
SSy131 = SSy16 + (4*Rx - TSx)*tand(beta);
TSy16 = SSy131;
TSy131 = SSy131 + (4*Rx - TSx)*tand(alpha);
SSy132 = TSy131 - (4*Rx - TSx)*tantheta;
TSy132 = SSy132 + (4*Rx - TSx)*tand(alpha);
SSy133 = TSy132 - (4*Rx - TSx)*tantheta;
SSy_r1 = [SSy11 SSy12 SSy13 SSy14 SSy15 SSy16 SSy131 SSy132 SSy133];
TSy_r1 = [TSy11 TSy12 TSy13 TSy14 TSy15 TSy16 TSy131 TSy132];
k = 1;
for j=1:len(2)
    if unit_config(j) == 41
        y_shift = y_shift - Rx*cr_yscales;
        ufpy(j) = y_shift;
        y_shift = y_shift + Rx*cr_yscales;

    elseif unit_config(j) == 32
        y_shift = y_shift + 0.5*SSy_r1(k);
        ufpy(j) = y_shift;
        y_shift = y_shift + 0.5*SSy_r1(k);

    elseif unit_config(j) == 33
        y_shift = y_shift + Rx*cr_yscale;
        ufpy(j) = y_shift;
        y_shift = y_shift - Rx*cr_yscale;

    elseif unit_config(j) == 34
        y_shift = y_shift - 0.5*TSy_r1(k);
        ufpy(j) = y_shift;
        y_shift = y_shift - 0.5*TSy_r1(k);
        k = k + 1;

    elseif unit_config(j) == 31
        y_shift = y_shift - Rx*cr_yscale;
        ufpy(j) = y_shift;
        y_shift = y_shift + Rx*cr_yscale;

    elseif unit_config(j) == 43
        y_shift = y_shift + Rx*cr_yscales;
        ufpy(j) = y_shift;
        y_shift = y_shift - Rx*cr_yscales;
    end

    if unit_config(j) == 41 || unit_config(j) == 31 || unit_config(j) == 43 || unit_config(j) == 33
        x_shift = x_shift + Rx;
        ufpx(j) = x_shift;
        x_shift = x_shift + Rx;

    elseif unit_config(j) == 32
        ufpx(j) = x_shift;
    
    elseif unit_config(j) == 34
        x_shift = x_shift - 0.5*TSx;
        ufpx(j) = x_shift;
        x_shift = x_shift - 0.5*TSx;
    end
end
unit_feature_pos = [ufpx; ufpy];
end


function [unit_feature_pos, TSy_r2, SSy_r2] = feature_pos_unit_top2(unit_config, Rx, cr_yscale, cr_yscales, TSx, TSy, SSy, alpha, beta)
ufpx = [];
ufpy = [];
len = size(unit_config);
y_shift = -Rx*cr_yscales;
x_shift = -Rx;

% Deriving the lengths of each of the row 2 TSy's and SSys
TSy26 = TSy;
SSy26 = TSy26 - (4*Rx - TSx)*tand(alpha);
TSy25 = SSy26 + (4*Rx - TSx)*tand(beta);
SSy25 = TSy25 - (4*Rx - TSx)*tand(alpha);
TSy24 = SSy25 + (4*Rx - TSx)*tand(beta);
SSy24 = TSy24 - (4*Rx - TSx)*tand(alpha);
TSy23 = SSy24 + (4*Rx - TSx)*tand(beta);
SSy23 = TSy23 - (4*Rx - TSx)*tand(alpha);
TSy22 = SSy23 + (4*Rx - TSx)*tand(beta);
SSy22 = TSy22 - (4*Rx - TSx)*tand(alpha);
TSy21 = SSy22 + (4*Rx - TSx)*tand(beta);
SSy21 = TSy21 - (4*Rx - TSx)*tand(alpha);
TSy132 = SSy21 + (4*Rx - TSx)*tand(beta);
TSy_r2 = [TSy132 TSy21 TSy22 TSy23 TSy24 TSy25 TSy26];
SSy_r2 = [SSy21 SSy22 SSy23 SSy24 SSy25 SSy26];
k = 1;
for j=1:len(2)
    if unit_config(j) == 41
        y_shift = y_shift - Rx*cr_yscales;
        ufpy(j) = y_shift;
        y_shift = y_shift + Rx*cr_yscales;

    elseif unit_config(j) == 32
        y_shift = y_shift + 0.5*SSy_r2(k);
        ufpy(j) = y_shift;
        y_shift = y_shift + 0.5*SSy_r2(k);
        k = k + 1;

    elseif unit_config(j) == 33
        y_shift = y_shift + Rx*cr_yscale;
        ufpy(j) = y_shift;
        y_shift = y_shift - Rx*cr_yscale;

    elseif unit_config(j) == 34
        y_shift = y_shift - 0.5*TSy_r2(k);
        ufpy(j) = y_shift;
        y_shift = y_shift - 0.5*TSy_r2(k);

    elseif unit_config(j) == 31
        y_shift = y_shift - Rx*cr_yscale;
        ufpy(j) = y_shift;
        y_shift = y_shift + Rx*cr_yscale;

    elseif unit_config(j) == 43
        y_shift = y_shift + Rx*cr_yscales;
        ufpy(j) = y_shift;
        y_shift = y_shift - Rx*cr_yscales;
    end

    if unit_config(j) == 41 || unit_config(j) == 31 || unit_config(j) == 43 || unit_config(j) == 33
        x_shift = x_shift + Rx;
        ufpx(j) = x_shift;
        x_shift = x_shift + Rx;

    elseif unit_config(j) == 32
        ufpx(j) = x_shift;
    
    elseif unit_config(j) == 34
        x_shift = x_shift - 0.5*TSx;
        ufpx(j) = x_shift;
        x_shift = x_shift - 0.5*TSx;
    end
end
unit_feature_pos = [ufpx; ufpy];
end

function [unit_feature_pos] = feature_pos_unit(unit_config, Rx, cr_yscale, cr_yscales, TSx, TSy, SSy)
ufpx = [];
ufpy = [];
len = size(unit_config);
y_shift = Rx*cr_yscales;
x_shift = -Rx;
for j=1:len(2)
    if unit_config(j) == 11
        y_shift = y_shift - Rx*cr_yscales;
        ufpy(j) = y_shift;
        y_shift = y_shift + Rx*cr_yscales;

    elseif unit_config(j) == 2
        y_shift = y_shift + 0.5*SSy;
        ufpy(j) = y_shift;
        y_shift = y_shift + 0.5*SSy;

    elseif unit_config(j) == 3
        y_shift = y_shift + Rx*cr_yscale;
        ufpy(j) = y_shift;
        y_shift = y_shift - Rx*cr_yscale;

    elseif unit_config(j) == 4
        y_shift = y_shift - 0.5*TSy;
        ufpy(j) = y_shift;
        y_shift = y_shift - 0.5*TSy;

    elseif unit_config(j) == 1
        y_shift = y_shift - Rx*cr_yscale;
        ufpy(j) = y_shift;
        y_shift = y_shift + Rx*cr_yscale;

    elseif unit_config(j) == 13
        y_shift = y_shift + Rx*cr_yscales;
        ufpy(j) = y_shift;
        y_shift = y_shift - Rx*cr_yscales;
    end

    if unit_config(j) == 11 || unit_config(j) == 1 || unit_config(j) == 13 || unit_config(j) == 3
        x_shift = x_shift + Rx;
        ufpx(j) = x_shift;
        x_shift = x_shift + Rx;

    elseif unit_config(j) == 2
        ufpx(j) = x_shift;
    
    elseif unit_config(j) == 4
        x_shift = x_shift - 0.5*TSx;
        ufpx(j) = x_shift;
        x_shift = x_shift - 0.5*TSx;
    end
end
unit_feature_pos = [ufpx; ufpy];
end

function [crown_points] = crown_aggregator_V3(parsec_points_left,parsec_points_right, circular_points, type, cr_yscale)
    if type == "standard"
        circular_points(:,1) = circular_points(:,1).*cr_yscale;% Scaling circular points in y direction. Parsec has already been shifted accordingly
        crown_points_right = [circular_points;parsec_points_right];
        unflipped_crown_points_left = [circular_points];
        crown_points_left = flip(unflipped_crown_points_left);
        crown_points_left(:,2) = -1*crown_points_left(:,2);
        crown_points = [crown_points_left; crown_points_right];
    end
end

function [crown_points] = crown_aggregator_V4(parsec_points_left,parsec_points_right, circular_points, type, cr_yscale_shifted, cr_yscale, Rx)
    if type == "shifted"
        circular_points(:,1) = circular_points(:,1).*cr_yscale_shifted;% Scaling circular points in y direction. Parsec has already been shifted accordingly
        circular_points(:,1) = circular_points(:,1) - (cr_yscale - cr_yscale_shifted)*Rx; % shifts circular points down to accomodate the required shift
        crown_points_right = [circular_points;parsec_points_right];
        unflipped_crown_points_left = [circular_points];
        crown_points_left = flip(unflipped_crown_points_left);
        crown_points_left(:,2) = -1*crown_points_left(:,2);
        crown_points = [crown_points_left; crown_points_right];
    end
end

function [crown_points] = generate_crown_parsec(A, type, np_crown, x_te, rc, cr_yscale)
if type == "top"
    x = rc:((x_te-rc)/np_crown):x_te;
    x_mat = [x.^(0.5);x.^(1.5);x.^(2.5);x.^(3.5);x.^(4.5);x.^(5.5)];
    A_mat = A';
    z = A_mat*x_mat;
    x = x + rc*(cr_yscale - 1);% Shifing y direction parsec to accomodate for scaling y direction circular points
    crown_points = [-x' z(1,:)'];
end

if type == "bottom"
    x = rc:((x_te-rc)/np_crown):x_te;
    x_mat = [x.^(0.5);x.^(1.5);x.^(2.5);x.^(3.5);x.^(4.5);x.^(5.5)];
    A_mat = A';
    z = A_mat*x_mat;
    x = x + rc*(cr_yscale - 1);% Shifing y direction parsec to accomodate for scaling y direction circular points
    crown_points = [-x' z(1,:)'];
end
end

function [A] = parsec(rc, x_te, z_te, z_te_slope)
M = [1 0 0 0 0 0; 0.5*(rc)^-0.5 1.5*(rc)^0.5 2.5*(rc)^1.5 3.5*(rc)^2.5 4.5*(rc)^3.5 5.5*(rc)^4.5;(0.5)*(-0.5)*rc^(-1.5) (1.5)*(0.5)*rc^(-0.5) (2.5)*(1.5)*rc^(0.5) (3.5)*(2.5)*rc^(1.5) (4.5)*(3.5)*rc^(2.5) (5.5)*(4.5)*rc^(3.5); rc^0.5 rc^1.5 rc^2.5 rc^3.5 rc^4.5 rc^5.5; x_te^0.5 x_te^1.5 x_te^2.5 x_te^3.5 x_te^4.5 x_te^5.5; 0.5*(x_te)^-0.5 1.5*(x_te)^0.5 2.5*(x_te)^1.5 3.5*(x_te)^2.5 4.5*(x_te)^3.5 5.5*(x_te)^4.5];
Q = [sqrt(2*rc) -sqrt(2*rc); 0 0;(-1/rc) (-1/rc); rc rc; z_te z_te; z_te_slope z_te_slope];
A = inv(M)*Q;
end

function [corrected_strut, corrected_slope] = strut_slope_correction(strut, slopes, cutoff)
[~,top_index] = min(abs(strut(:,2) + cutoff));
[~,bottom_index] = min(abs(strut(:,2) - cutoff));
corrected_strut = strut(top_index:bottom_index,:);
corrected_slope = slopes(top_index - 1: bottom_index);
end

function slopes = strut_slopes(strut, np_strut)
    slopes = (strut([2:np_strut+1],2) - strut([1:np_strut],2))./(strut([2:np_strut+1],1) - strut([1:np_strut],1));
end

function [short_strut] = short_strut_generator(SSy, np_strut)
    x_pts = zeros(1,np_strut+1);
    y_pts = 0:SSy/np_strut:SSy;
    y_pts = y_pts - 0.5*SSy;
    short_strut = [x_pts' y_pts'];
end

function [tall_strut] = tall_strut_generator(TSx, TSy, np_strut)
    x_pts = 0:(TSx/(np_strut)):TSx;
    y_pts = x_pts.*(TSy/TSx);
    x_pts = x_pts - 0.5*TSx;
    y_pts = y_pts - 0.5*TSy;
    tall_strut = [x_pts' y_pts'];
end


function [crown_top] = crown_top_circular(cr, np_crown)
    x_mid_pts = 0:(cr/(np_crown)):cr;
    y_mid_pts = sqrt(cr^2 - x_mid_pts.^2);
    y_mid_pts = y_mid_pts - cr;
    crown_top = [y_mid_pts' x_mid_pts'];
end


%% Solidworks integration related functions
function [data_points] = spline_data_points_3d(spline_pts, Rc) % Convert points in [x1 y1;x2 y2] format to [x1 yz1 z1 x2 y2 z2] format
    theta = spline_pts(:,1)/Rc;
    spline_pts_mod(:,1) = Rc*sin(theta);
    spline_pts_mod(:,2) = -Rc*cos(theta);
    spline_pts_mod(:,3) = spline_pts(:,2);
    spline_pts_trans = spline_pts_mod';
    data_points = reshape(spline_pts_trans, 1, []);
    data_points = data_points/1000;
    data_points = round(data_points, 7);
end
function [data_points] = spline_data_points(spline_pts) % Convert points in [x1 y1;x2 y2] format to [x1 yz1 z1 x2 y2 z2] format
    spline_pts(:,3) = 0;
    spline_pts_trans = spline_pts';
    data_points = reshape(spline_pts_trans, 1, []);
    data_points = data_points/1000;
end

function [cleaned_spline_points] = removing_close_neibours(spline_points, distance_thresh)
cleaned_spline_points = [];
cleaned_spline_points(1,:) = spline_points(1,:);
last_point = spline_points(end,:);
for i=2:1:(size(spline_points,1) - 1)
    distance_prev = sqrt((spline_points(i,1) - cleaned_spline_points(end,1))^2 + (spline_points(i,2) - cleaned_spline_points(end,2))^2);
    distance_last = sqrt((spline_points(i,1) - last_point(1,1))^2 + (spline_points(i,2) - last_point(1,2))^2);
    if distance_prev > distance_thresh && distance_last > distance_thresh
        cleaned_spline_points(end+1,:) = spline_points(i,:);
    end
end
cleaned_spline_points(end+1,:) = last_point;
end