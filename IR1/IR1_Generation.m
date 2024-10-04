function IR1_Generation(x_c, result_location)
%% Generate x for candidate generation 
% Base values
x_b(1) = 0.07; % Row 1 NURBS Control Points Width 
x_b(2) = 0.8; % Row 1 NURBS Control Points Height
x_b(3) = 0.61317;%ww(1)
x_b(4) = 0.09022;%ww(2)
x_b(5) = 0.86459;%ww(3)
x_b(6) = 0.36024;%ww(4)
x_b(7) = 0.458; %Stent Crimped radius
x_b(8) = 2; % Half of (Crown multiple of Second row unit structure(C-SS-TS) to generate the number of crowns of first row faced by each crown of second row - HAS TO BE EVEN to generate a viable top row, therefore taking half here, also only 1, 2 are viable options (3) is too conjusted. and () is too big, fixing it to 2
x_b(9) = 0.9; % Row 2 Tall Strut height 
x_b(10) = 0.7; %Row 2 short strut height 
x_b(11) = 21/17; %Row 2 Short Strut top crown radius / Row 2 other crown radius
x_b(12) = 0.095; % Connector lateral crown radius
x_b(13) = 0.16; % Connector lateral length
x_b(14) = 3; % Total Number of repittions of SS-TS-C unit
x_b(15) = 18; % Stent Length
x_b(16) = 13; % Number of rows in the stent
x_b(17) = 1;% Stent width factor (1) will be multiplied with base value before sending to the geometry creater
x_b(18) = 1;% Stent thicknesss factor (1) will be multiplied with base value before sending to the geometry creater

% wt dimensions
wt_base = 0.08;
wt_tight = 0.065;

% Separating into independent Variables to be be optimised and fixed
[~, x_fixed] = Divide_xx(x_b);
x_candidate_ind = x_c;
x_candidate_ind(13) = round(x_candidate_ind(13));
x_candidate_ind(15) = round(x_candidate_ind(15));
x_candidate = Recreate_xx(x_candidate_ind, x_fixed);
if x_candidate(14) == 4
    x_candidate(17:18) =  wt_tight*x_candidate(17:18);
else
    x_candidate(17:18) =  wt_base*x_candidate(17:18);
end
x = x_candidate;

%% Start Solidworks
py.sld_interface.start();

%% Associating input variables to  variable names
stent_width = x(17);
stent_length = x(15); % Total length of the stent
stent_thickness = x(18); 
rnos = x(16); % Total number of strut rows in the stent
crimped_cylinder_radius = x(7);
ncon = x(14);

% Fixed parameters
np_strut = 500; % Default was 500. Changed to 100. Can be varied to see if model works better
np_crown = 30; % Default was 201. Changed to 41. Can be varied to see if model works better

%% Row 1 Creation 
% Row 1 dependent variables 
cnos = x(8)*2*x(14); % Total Number of crowns in row 1
cr = 2*pi*x(7)/(cnos*2); % crown radius for top crowns
N_rep_r1 = cnos/2;% for the first row, the number of repitions of base sketch will be half of the number of crowns in the row
nc = cnos/2; % Numbe of crowns on each side
rc = cr;
cut_off_distance_for_struts = 0.05; %Can be varied to see if model works better. 0.05 was initial

% Main Nurbs parameters
k = 7;
n = 7;
nurbs_half_width = x(1)/2;
nurbs_half_height = x(2)/2;
P = [0 nurbs_half_width nurbs_half_width nurbs_half_width -1*nurbs_half_width -1*nurbs_half_width -1*nurbs_half_width 0; -1*nurbs_half_height -1*nurbs_half_height -1*nurbs_half_height/2 0 0 nurbs_half_height/2 nurbs_half_height nurbs_half_height];
t = [0 0 0 0 0 0 0 1 2 3 3 3 3 3 3 3];
x_nurbs = [x(3) x(4) x(5) x(6)];

% Recreating Row 1 strut
strut_half_height = nurbs_half_height;
strut_length_r1 = 2*strut_half_height;
slope_correction_cutoff_value = strut_half_height - cut_off_distance_for_struts;
stent_r1_config = Stent_r1_builder(nc);
D1_strut = Simple_NURBS(x_nurbs,P,t,k,n,np_strut);
D1_slopes = strut_slopes(D1_strut, np_strut);
[D1_corrected_strut, D1_corrected_slope] = strut_slope_correction(D1_strut, D1_slopes, slope_correction_cutoff_value);
D2_strut = [-1*D1_strut(:,1) D1_strut(:,2)];
D2_corrected_strut = [-1*D1_corrected_strut(:,1) D1_corrected_strut(:,2)];

% generating the parsec points for row 1 top crowns 
z_te = D1_corrected_strut(end,1) + rc;
x_te = -(D1_corrected_strut(end,2) - strut_half_height - rc);
z_te_slope = -1/D1_corrected_slope(end,1);
parsec_coeffs_top_A = parsec(rc, x_te, z_te, z_te_slope);
crown_points_parsec_top = generate_crown_parsec(parsec_coeffs_top_A, "top", np_crown/4,x_te,rc);

% generationg the semicircular crown points for top
crown_points_circular = crown_top_circular(rc, 4*np_crown);

% Aggregating parsec and semicircular crown points
crown_top = crown_aggregator_V2(crown_points_parsec_top, crown_points_parsec_top, crown_points_circular,'standard');

% generating the parsec points for bottom crown
z_te_b = -D1_corrected_strut(1,1) + rc;
x_te_b = (D1_corrected_strut(1,2) + strut_half_height + rc);
z_te_slope_b = -1/D1_corrected_slope(1,1);
parsec_coeffs_bottom = parsec(rc, x_te_b, z_te_b, z_te_slope_b);
crown_points_parsec_bottom = generate_crown_parsec(parsec_coeffs_bottom, "bottom", np_crown/4,x_te_b,rc);
crown_points_circular_bottom = crown_top_circular(rc, 4*np_crown);
crown_bottom_rev = crown_aggregator_V2(crown_points_parsec_bottom, crown_points_parsec_bottom, crown_points_circular_bottom,'standard');
crown_bottom = [-1*crown_bottom_rev(:,1) crown_bottom_rev(:,2)];

% Building Row 1
feature_pos = feature_position(stent_r1_config,cr);% centre point for crowns and exact x point for strut
FR = first_row(stent_r1_config, feature_pos, D1_corrected_strut, D2_corrected_strut, crown_top, crown_bottom,strut_half_height,cr); 

%% Rows 2-end creation
% Row 2 parameters
% Terminology; SS - Short Strut system , SSS - Short Strut System Strut, 
% TS - Tall Strut system, TSS - Tall Strut System Strut
% C - Connector D1 - Strut Direction 1, D2 - Strut Direction 2, TC - Top Crown, 
% BC - Bottom Crow, CS - Connector Strut, CC1 - Connector curve 1,   
% CS1 - COnnector Straight 1, CC - Connector Crow, CS2 - Connector straight
% CC2 - Connector Crown 2, OC - Other crowns (Except SS TC), 
% alpha - Radius(SSTC)/Radius (other crowns) ~ 21/17 for ir1 stent, 
% beta - SSS_length / Row 1 strut length ~ 7/8 for ir1 stent
% gamma - Row 1 strut length / TSS_length ~ 8/9 for ir1 stent
% con_Rad - Connector Crown Radius

% Parameters for Row 2
cnos2 = x(14)*5; % Total number of strut crowns in row 2
N_rep = cnos2/5; % Total number of expected repitions of SS-TS-C pairs
cnos2_sstc = N_rep; % Total number of short strut top crowns
cnos2_oc = cnos2 - cnos2_sstc; % Total Number of other crowns
alpha = x(11);% Defined in parameter comments
beta = x(10)/x(2);% Defined in parameter comments 
gamma = x(2)/x(9); % Defined in parameter comments
con_Rad = x(12);% Defined in parameter comments
con_straight_length = x(13);% Length of straight lateral portion in the connector.
row_1_strut_length = 2*strut_half_height;
OC_Rad = cr*cnos/(cnos2_oc + alpha*cnos2_sstc);
SSTC_Rad = OC_Rad*alpha;

% Create Short Strut
r2_SSD1_slope_correction_cutoff_value = strut_half_height*beta - cut_off_distance_for_struts;
r2_SSD1_strut = new_length_strut_builder(D1_strut, beta);
r2_SSD1_slopes = strut_slopes(r2_SSD1_strut, np_strut);

% Code 214
[r2_SSD1_corrected_strut, r2_SSD1_corrected_slope] = strut_slope_correction(r2_SSD1_strut, r2_SSD1_slopes, r2_SSD1_slope_correction_cutoff_value);

% Code - 212
r2_SSD2_corrected_strut = [-1*r2_SSD1_corrected_strut(:,1) r2_SSD1_corrected_strut(:,2)];

% generating the parsec for row 2 short strut top crowns 
r2_SS_z_te = r2_SSD1_corrected_strut(end,1) + SSTC_Rad;
r2_SS_x_te = -(r2_SSD1_corrected_strut(end,2) - strut_half_height*beta - SSTC_Rad);
r2_SS_z_te_slope = -1/r2_SSD1_corrected_slope(end,1);
r2_SS_parsec_coeffs_top_A = parsec(SSTC_Rad, r2_SS_x_te, r2_SS_z_te, r2_SS_z_te_slope);
r2_SS_crown_points_parsec_top = generate_crown_parsec(r2_SS_parsec_coeffs_top_A, "top", np_crown/4,r2_SS_x_te,SSTC_Rad);

% generationg the semicircular crown points for top
r2_SS_crown_points_circular = crown_top_circular(SSTC_Rad, 4*np_crown);

% Aggregating parsec and semicircular crown points for top crown 
% Code - 213
r2_SS_crown_top = crown_aggregator_V2(r2_SS_crown_points_parsec_top, r2_SS_crown_points_parsec_top, r2_SS_crown_points_circular, 'standard');

% generating the parsec for Short Strut bottom.
r2_SS_z_te_b = -r2_SSD1_corrected_strut(1,1) + OC_Rad;
r2_SS_x_te_b = (r2_SSD1_corrected_strut(1,2) + strut_half_height*beta + OC_Rad);
r2_SS_z_te_slope_b = -1/r2_SSD1_corrected_slope(1,1);
r2_SS_parsec_coeffs_bottom = parsec(OC_Rad, r2_SS_x_te_b, r2_SS_z_te_b, r2_SS_z_te_slope_b);

% generationg the semicircular crown points for standard bottom crown 
r2_SS_crown_points_parsec_bottom = generate_crown_parsec(r2_SS_parsec_coeffs_bottom, "bottom", np_crown/4,r2_SS_x_te_b,OC_Rad);
r2_crown_points_circular_bottom = crown_top_circular(OC_Rad, 4*np_crown);

% Create tall strut 
r2_TSD1_slope_correction_cutoff_value = (strut_half_height/gamma) - cut_off_distance_for_struts;
r2_TSD1_strut = new_length_strut_builder(D1_strut, 1/gamma);
r2_TSD1_slopes = strut_slopes(r2_TSD1_strut, np_strut);

% Code 223
[r2_TSD1_corrected_strut, r2_TSD1_corrected_slope] = strut_slope_correction(r2_TSD1_strut, r2_TSD1_slopes, r2_TSD1_slope_correction_cutoff_value);

% Code - 221
r2_TSD2_corrected_strut = [-1*r2_TSD1_corrected_strut(:,1) r2_TSD1_corrected_strut(:,2)];

% generating the parsec for row 2 tall strut top crowns 
r2_TS_z_te = r2_TSD1_corrected_strut(end,1) + OC_Rad;
r2_TS_x_te = -(r2_TSD1_corrected_strut(end,2) - (strut_half_height/gamma) - OC_Rad);
r2_TS_z_te_slope = -1/r2_TSD1_corrected_slope(end,1);
r2_TS_parsec_coeffs_top_A = parsec(OC_Rad, r2_TS_x_te, r2_TS_z_te, r2_TS_z_te_slope);
r2_TS_crown_points_parsec_top = generate_crown_parsec(r2_TS_parsec_coeffs_top_A, "top", np_crown/4,r2_TS_x_te,OC_Rad);

% generationg the semicircular crown points for top
r2_TS_crown_points_circular = crown_top_circular(OC_Rad, 4*np_crown);

% Aggregating parsec and semicircular crown points for top crown 
% Code - 222
r2_TS_crown_top = crown_aggregator_V2(r2_TS_crown_points_parsec_top, r2_TS_crown_points_parsec_top, r2_TS_crown_points_circular, 'standard');

% generating the parsec for Tall Strut bottom
r2_TS_z_te_b = -r2_TSD1_corrected_strut(1,1) + OC_Rad;
r2_TS_x_te_b = (r2_TSD1_corrected_strut(1,2) + (strut_half_height/gamma) + OC_Rad);
r2_TS_z_te_slope_b = -1/r2_TSD1_corrected_slope(1,1);
r2_TS_parsec_coeffs_bottom = parsec(OC_Rad, r2_TS_x_te_b, r2_TS_z_te_b, r2_TS_z_te_slope_b);

% generationg the semicircular crown points for Tall strut bottom crown 
r2_TS_crown_points_parsec_bottom = generate_crown_parsec(r2_TS_parsec_coeffs_bottom, "bottom", np_crown/4,r2_TS_x_te_b,OC_Rad);
r2_SS_TS_crown_bottom_rev = crown_aggregator_V2(r2_SS_crown_points_parsec_bottom, r2_TS_crown_points_parsec_bottom, r2_crown_points_circular_bottom, 'standard');

% Code - 215 
% Crown bottom between SS and TS
r2_SS_TS_crown_bottom = [-1*r2_SS_TS_crown_bottom_rev(:,1) r2_SS_TS_crown_bottom_rev(:,2)];

% Create Connector
connector_strut_length = ((stent_length - 2*rc - strut_length_r1)/(rnos - 1)) - 4*con_Rad - OC_Rad;
connector_strut = new_length_strut_builder(D1_strut, (connector_strut_length/(2*strut_half_height)));
connector_strut_slopes = strut_slopes(connector_strut, np_strut);
connector_slope_correction_cutoff_value = (connector_strut_length/2) - cut_off_distance_for_struts;

% Code 231
[connector_corrected_strut, connector_corrected_slope] = strut_slope_correction(connector_strut, connector_strut_slopes, connector_slope_correction_cutoff_value);

% generating the parsec for connector strut top half crown 
connector_top_z_te = connector_corrected_strut(end,1) + con_Rad;
connector_top_x_te = -(connector_corrected_strut(end,2) - (connector_strut_length/2) - con_Rad);
connector_top_z_te_slope = -1/connector_corrected_slope(end,1);
connector_top_parsec_coeffs_top_A = parsec(con_Rad, connector_top_x_te, connector_top_z_te, connector_top_z_te_slope);
connector_top_crown_points_parsec_top = generate_crown_parsec(connector_top_parsec_coeffs_top_A, "top", np_crown/4,connector_top_x_te,con_Rad);
connector_top_crown_points_parsec_top_left = connector_top_crown_points_parsec_top;
connector_top_crown_points_parsec_top_left(:,2) = connector_top_crown_points_parsec_top_left(:,2) - 2*con_Rad;

% generationg the semicircular crown points for top
connector_top_crown_points_circular = crown_top_circular(con_Rad, 4*np_crown);

% Aggregating parsec and semicircular crown points for top crown 
% Code - 232
connector_top_cc1_top = crown_aggregator_V2(connector_top_crown_points_parsec_top_left, [], connector_top_crown_points_circular, 'connector_cc1');
connector_top_cc2_rev = crown_aggregator_V2([], [], connector_top_crown_points_circular, 'connector_cc2');

% Code - 236
connector_top_cc2 = [-1*connector_top_cc2_rev(:,1) connector_top_cc2_rev(:,2)];

% Code - 233 and 235
connector_straight_cs1_cs2 = connector_top_straight(con_straight_length,np_strut);
connector_lateral_crown_rotated = crown_aggregator_V2([],[],connector_top_crown_points_circular, 'standard');

% Code 234
connector_lateral_crown = [connector_lateral_crown_rotated(:,2) connector_lateral_crown_rotated(:,1)];

% Create bottom crowns between connectors and struts - Codes 211 and 224
% generating the parsec for Connector bottom crown
connector_z_te_b = -connector_corrected_strut(1,1) + OC_Rad;
connector_x_te_b = (connector_corrected_strut(1,2) + (connector_strut_length/2) + OC_Rad);
connector_z_te_slope_b = -1/connector_corrected_slope(1,1);
connector_parsec_coeffs_bottom = parsec(OC_Rad, connector_x_te_b, connector_z_te_b, connector_z_te_slope_b);

% generationg the semicircular crown points for standard bottom crown 
connector_crown_points_parsec_bottom = generate_crown_parsec(connector_parsec_coeffs_bottom, "bottom", np_crown/4,connector_x_te_b,OC_Rad);
TS_connector_crown_bottom_rev = crown_aggregator_V2(r2_TS_crown_points_parsec_bottom,r2_TS_crown_points_parsec_bottom, r2_crown_points_circular_bottom, "ts_connector");

% Code 224
TS_connector_crown_bottom = [-1*TS_connector_crown_bottom_rev(:,1) TS_connector_crown_bottom_rev(:,2)];

% Code 211
connector_SS_crown_bottom_rev = crown_aggregator_V2(connector_crown_points_parsec_bottom, r2_SS_crown_points_parsec_bottom,r2_crown_points_circular_bottom, 'standard');
connector_SS_crown_bottom = [-1*connector_SS_crown_bottom_rev(:,1) connector_SS_crown_bottom_rev(:,2)];

% Creating classes with helix structures for second row
SS = row2_SS(connector_SS_crown_bottom, r2_SSD2_corrected_strut, r2_SS_crown_top, r2_SSD1_corrected_strut, r2_SS_TS_crown_bottom, SSTC_Rad, OC_Rad, beta);
TS = row2_TS(r2_TSD2_corrected_strut, r2_TS_crown_top, r2_TSD1_corrected_strut, TS_connector_crown_bottom, OC_Rad, gamma);
C = row2_C(connector_corrected_strut, connector_top_cc1_top, connector_straight_cs1_cs2, connector_lateral_crown, connector_straight_cs1_cs2, connector_top_cc2, connector_strut_length,con_straight_length, con_Rad, OC_Rad);

% Building all rows
row_2_con_x_pos = 3*cr; % X position of the first connector of the row 2
row_2_bottom_y_pos = -strut_half_height - cr - 4*C.tc_r - C.length - C.bc_r;% Y Positiion of the bottom of the row 2
stent_r2_config = Stent_r2_builder(cnos2);
feature_pos_r2 = feature_position_r2(stent_r2_config, SSTC_Rad, OC_Rad);
MR = main_rows(stent_r2_config,feature_pos_r2, SS, TS, C);
RowFig = figure();
RowAx = axes('Parent', RowFig);
Stent = rows_builder(MR,FR, rnos, row_2_con_x_pos, row_2_bottom_y_pos, stent_width, stent_thickness,stent_length, crimped_cylinder_radius, N_rep_r1, N_rep, result_location, ncon, cnos, RowAx);

%% Saving figure, stent solidworks file, parasolid file and Mesh STL file 

% Saving 2D image
RowAx.XLim = [-1 13];
RowAx.XTick = -1:1:13;
RowAx.YLim = [-20,2];
RowAx.YTick = -20:2:2;
set(RowFig, 'Position', [100, 100, 1400*1.5, 2200*1.5]);
%figname = "2d_stent.png";
%fullpath = fullfile(result_location,figname);
%saveas(RowFig,fullpath);
close(RowFig);

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
py.sld_interface.ps_to_stl(fullpath_ps,fullpath_stl);
py.sld_interface.exit_sw();

% Stent built. 
end

%% Helper functions

function reloadPy()
clear classes;
m = py.importlib.import_module('sld_interface');
py.importlib.reload(m);
end

function [stent_data] = rows_builder(MR, FR, rnos, row_2_con_x_pos, row_2_bottom_y_pos, stent_width, stent_thickness, stent_length, crimped_cylinder_radius, N_rep_r1, N_rep,result_location, ncon, cnos, RowAx)
x_shift = row_2_con_x_pos;
y_shift = row_2_bottom_y_pos;
stent_data = [];

% Building row 2 through wrapping 
cylinder_name = py.sld_interface.create_base_cylinder(crimped_cylinder_radius/1000, stent_length/1000);% Create the base crimped cylinder of crimped radius
py.sld_interface.create_base_cylinder_axis();% Create cylinder axis
[row, wrap_name] = main_rows_builder(MR.stent_r2_config,MR.feature_pos_r2,MR.SS,MR.TS,MR.C,FR.strut_half_height, x_shift, y_shift, stent_width, stent_thickness, crimped_cylinder_radius, ncon, 1, RowAx);
x_shift = x_shift + 3*MR.SS.bc_r + 2*MR.SS.tc_r;
y_shift = y_shift - 4*MR.C.tc_r - MR.C.length - MR.C.bc_r;
rel_y_shift = - 4*MR.C.tc_r - MR.C.length - MR.C.bc_r;% y shift in the subsequent row as compared to previous row.
rel_theta_shift = 2*pi*(3*MR.SS.bc_r + 2*MR.SS.tc_r)/(2*pi*crimped_cylinder_radius);% Theta shift in subsequent rows as compared to previous rows in Radians
stent_data = [stent_data;row];
py.sld_interface.zoom_fit();
cylinder_name_2 = py.sld_interface.create_base_cylinder(crimped_cylinder_radius/1000, stent_length/1000);
stent_part = py.sld_interface.combine_subtract(cylinder_name_2,wrap_name);
current_row = py.sld_interface.circpattern(stent_part,"Axis1",N_rep);
full_stent = current_row;
previous_row = current_row;
py.sld_interface.trimetric_zoom_fit();
for i=3:rnos% goes to rnos usually
    [row, wrap_name] = main_rows_builder(MR.stent_r2_config,MR.feature_pos_r2,MR.SS,MR.TS,MR.C,FR.strut_half_height, x_shift, y_shift, stent_width, stent_thickness,crimped_cylinder_radius, ncon, 0, RowAx);
    x_shift = x_shift + 3*MR.SS.bc_r + 2*MR.SS.tc_r;
    y_shift = y_shift - 4*MR.C.tc_r - MR.C.length - MR.C.bc_r;
    stent_data = [stent_data;row];
    py.sld_interface.zoom_fit();
    current_row = py.sld_interface.copy_move(previous_row, "Axis1", rel_y_shift/1000,rel_theta_shift);
    if i ~= 3
        full_stent = py.sld_interface.row_combine(previous_row, full_stent);
    end
    previous_row = current_row;
    py.sld_interface.trimetric_zoom_fit();
end
full_stent = py.sld_interface.row_combine(previous_row, full_stent);

% Building row 1 now
cylinder_name = py.sld_interface.create_base_cylinder(crimped_cylinder_radius/1000, stent_length/1000);% Create the base crimped cylinder of crimped radius
[row, wrap_name] = row_1_builder(FR.stent_r1_config, FR.feature_pos, FR.D1_corrected_strut, FR.D2_corrected_strut, FR.crown_top, FR.crown_bottom,FR.strut_half_height,FR.cr, stent_width, stent_thickness, cnos, crimped_cylinder_radius, RowAx);
stent_data = [row;stent_data];
py.sld_interface.zoom_fit()
py.sld_interface.trimetric_zoom_fit()
cylinder_name_2 = py.sld_interface.create_base_cylinder(crimped_cylinder_radius/1000, stent_length/1000);
stent_one = py.sld_interface.combine_subtract(cylinder_name_2,wrap_name);
current_row = py.sld_interface.circpattern(stent_one,"Axis1",N_rep_r1);
current_row = py.sld_interface.row_combine(full_stent, current_row);
end

function [row, wrap_name] = main_rows_builder(stent_r2_config, feature_pos, SS, TS, C, strut_half_height, x_shift, y_shift, stent_width, stent_thickness,crimped_cylinder_radius, ncon, sw_cons, RowAx)
len = size(stent_r2_config);
row = [];
first_point = [];
last_point = [];
sketch_name_index = 1;
sketch_name = strings(0);
for j=1:len(2)
    % Starting with the Short Strut
    if stent_r2_config(j) == 211
        feature_x = SS.con_ss_bc(:,2) + feature_pos(j) + x_shift;
        feature_y = SS.con_ss_bc(:,1) + y_shift;        
        if sw_cons == 1
            sw_data_points = spline_data_points([feature_x feature_y]);
            py.sld_interface.create_spline(sw_data_points);
            first_point = [feature_x(1)/1000 feature_y(1)/1000];
            py.sld_interface.create_line(first_point,last_point_23);
            last_point = [feature_x(end)/1000 feature_y(end)/1000];
        end

    elseif stent_r2_config(j) == 212
        feature_x = SS.ss_d2(:,1) + feature_pos(j) + x_shift;
        feature_y = SS.ss_d2(:,2) + strut_half_height*SS.beta + SS.bc_r + y_shift;
        if sw_cons == 1
            sw_data_points = spline_data_points([feature_x feature_y]);
            py.sld_interface.create_spline(sw_data_points);
            first_point = [feature_x(1)/1000 feature_y(1)/1000];
            py.sld_interface.create_line(first_point,last_point);% Create joining line
            last_point = [feature_x(end)/1000 feature_y(end)/1000];
        end

    elseif stent_r2_config(j) == 213
        feature_x = SS.ss_tc(:,2) + feature_pos(j) + x_shift;
        feature_y = SS.ss_tc(:,1) + SS.bc_r + SS.tc_r + 2*strut_half_height*SS.beta + y_shift;
        if sw_cons == 1
            sw_data_points = spline_data_points([feature_x feature_y]);
            py.sld_interface.create_spline(sw_data_points);
            first_point = [feature_x(1)/1000 feature_y(1)/1000];
            py.sld_interface.create_line(first_point,last_point);% Create joining line
            last_point = [feature_x(end)/1000 feature_y(end)/1000];
        end

    elseif stent_r2_config(j) == 214
        feature_x = SS.ss_d1(:,1) + feature_pos(j) + x_shift;
        feature_y = SS.ss_d1(:,2) + strut_half_height*SS.beta + SS.bc_r + y_shift;
        if sw_cons == 1
            sw_data_points = spline_data_points([feature_x feature_y]);
            py.sld_interface.create_spline(sw_data_points);
            first_point = [feature_x(end)/1000 feature_y(end)/1000];
            py.sld_interface.create_line(first_point,last_point);% Create joining line
            last_point = [feature_x(1)/1000 feature_y(1)/1000];
        end

    elseif stent_r2_config(j) == 215
        feature_x = SS.ss_ts_bc(:,2) + feature_pos(j) + x_shift;
        feature_y = SS.ss_ts_bc(:,1) + y_shift;
        if sw_cons == 1
            sw_data_points = spline_data_points([feature_x feature_y]);
            py.sld_interface.create_spline(sw_data_points);
            first_point = [feature_x(1)/1000 feature_y(1)/1000];
            py.sld_interface.create_line(first_point,last_point);% Create joining line
            last_point = [feature_x(end)/1000 feature_y(end)/1000];
        end

    % Moving on to the tall strut
    elseif stent_r2_config(j) == 221
        feature_x = TS.ts_d2(:,1) + feature_pos(j) + x_shift;
        feature_y = TS.ts_d2(:,2) + (strut_half_height/TS.gamma) + TS.r + y_shift;
        if sw_cons == 1
            sw_data_points = spline_data_points([feature_x feature_y]);
            py.sld_interface.create_spline(sw_data_points);
            first_point = [feature_x(1)/1000 feature_y(1)/1000];
            py.sld_interface.create_line(first_point,last_point);% Create joining line
            last_point = [feature_x(end)/1000 feature_y(end)/1000];
        end

    elseif stent_r2_config(j) == 222
        feature_x = TS.ts_tc(:,2) + feature_pos(j) + x_shift;
        feature_y = TS.ts_tc(:,1) + 2*TS.r + 2*(strut_half_height/TS.gamma) + y_shift;
        if sw_cons == 1
            sw_data_points = spline_data_points([feature_x feature_y]);
            py.sld_interface.create_spline(sw_data_points);
            first_point = [feature_x(1)/1000 feature_y(1)/1000];
            py.sld_interface.create_line(first_point,last_point);% Create joining line
            last_point = [feature_x(end)/1000 feature_y(end)/1000];
        end

    elseif stent_r2_config(j) == 223
        feature_x = TS.ts_d1(:,1) + feature_pos(j) + x_shift;
        feature_y = TS.ts_d1(:,2) + (strut_half_height/TS.gamma) + TS.r + y_shift;
        if sw_cons == 1
            sw_data_points = spline_data_points([feature_x feature_y]);
            py.sld_interface.create_spline(sw_data_points);
            first_point = [feature_x(end)/1000 feature_y(end)/1000];
            py.sld_interface.create_line(first_point,last_point);% Create joining line
            last_point = [feature_x(1)/1000 feature_y(1)/1000];
        end

    elseif stent_r2_config(j) == 224
        feature_x = TS.ts_con_bc(:,2) + feature_pos(j) + x_shift;
        feature_y = TS.ts_con_bc(:,1) + y_shift;
        if sw_cons == 1
            sw_data_points = spline_data_points([feature_x feature_y]);
            py.sld_interface.create_spline(sw_data_points);
            first_point = [feature_x(1)/1000 feature_y(1)/1000];
            py.sld_interface.create_line(first_point,last_point);% Create joining line
            last_point = [feature_x(end)/1000 feature_y(end)/1000];
            last_point_23t(1) = last_point_23(1) + (2*pi()*crimped_cylinder_radius/(1000*ncon)); 
            last_point_23t(2) = last_point_23(2);
            py.sld_interface.create_line(last_point,last_point_23t);% Create the other joining line
            py.sld_interface.create_stent_segment_width(spline_id,stent_width/1000);%Create sketch width
            activesketch = py.sld_interface.get_active_sketch_name();
            sketch_name(sketch_name_index) = activesketch;
            sketch_name_index = sketch_name_index + 1;
            py.sld_interface.end_current_sketch();% Unit feature ends here so exit previous sketch
        end

    elseif stent_r2_config(j) == 231
        feature_x = C.cs(:,1) + feature_pos(j) + x_shift;
        feature_y = C.cs(:,2) + (C.length/2) + C.bc_r + y_shift;
        if sw_cons == 1
            sw_data_points = spline_data_points([feature_x feature_y]);
            py.sld_interface.create_new_sketch()% Unit feature starts here so create new sketch at Front Plane        
            spline_id = py.sld_interface.create_spline(sw_data_points);
            % straight line connection code here%
            last_point = [feature_x(end)/1000 feature_y(end)/1000];% last point for the connector side chain
            last_point_23 = [feature_x(1)/1000 feature_y(1)/1000];% last point for strut side chain
        end

    elseif stent_r2_config(j) == 232
        feature_x = C.cc1(:,2) + feature_pos(j) + C.tc_r + x_shift;
        feature_y = C.cc1(:,1) + (C.length) + C.bc_r + C.tc_r + y_shift;
        if sw_cons == 1
            sw_data_points = spline_data_points([feature_x feature_y]);
            py.sld_interface.create_spline(sw_data_points);
            first_point = [feature_x(end)/1000 feature_y(end)/1000];
            py.sld_interface.create_line(first_point,last_point);
            % straight line connection code here%
            %last_point = [feature_x(1)/1000 feature_y(1)/1000]; 
        end

    elseif stent_r2_config(j) == 233
        feature_x = C.cs1(:,1) + feature_pos(j) + C.tc_r + x_shift;
        feature_y = C.cs1(:,2) + C.length + C.bc_r + C.tc_r + y_shift;
        if sw_cons == 1
            sw_data_points = spline_data_points([feature_x feature_y]);
            py.sld_interface.create_spline(sw_data_points);
        end

    elseif stent_r2_config(j) == 234
        feature_x = C.cc(:,2) + feature_pos(j) + 2*C.tc_r + C.straight_length + x_shift;
        feature_y = C.cc(:,1) + (C.length) + C.bc_r + 2*C.tc_r + y_shift;
        if sw_cons == 1
            sw_data_points = spline_data_points([feature_x feature_y]);
            py.sld_interface.create_spline(sw_data_points);
        end

     elseif stent_r2_config(j) == 235
        feature_x = C.cs2(:,1) + feature_pos(j) + C.tc_r + x_shift;
        feature_y = C.cs2(:,2) + C.length + C.bc_r + 3*C.tc_r + y_shift;
        if sw_cons == 1
            sw_data_points = spline_data_points([feature_x feature_y]);
            py.sld_interface.create_spline(sw_data_points);
        end

     elseif stent_r2_config(j) == 236
        feature_x = C.cc2(:,2) + feature_pos(j) + C.tc_r + x_shift;
        feature_y = C.cc2(:,1) + (C.length) + C.bc_r + 3*C.tc_r + y_shift;
        if sw_cons == 1
            sw_data_points = spline_data_points([feature_x feature_y]);
            py.sld_interface.create_spline(sw_data_points); 
        end
    end

   feature = [feature_x feature_y];
   plot(RowAx, feature_x, feature_y,'LineWidth',1);
   hold(RowAx, "on")
   row = [row feature'];
end
    if sw_cons == 1
        no_sketches = size(sketch_name,2);
        wrapped = 0;
        for k=1:1:no_sketches
            if wrapped ~= 0
                dummy_var = "Dummy_var";
            else
                printing_sketc_name = sketch_name(k);
                [wrap_name] = py.sld_interface.wrap_sketch_around_base_cylinder_v2(sketch_name(k), stent_thickness/1000);
                if wrap_name ~= "NOWRAP"
                    wrapped = 1;
                end
            end
        end
    end
    if sw_cons == 0
        wrap_name = "NOSW";
    end    
row = row';
end

function [feature_pos_r2] = feature_position_r2(config,sstc_rad, oc_rad)
SSTC = [213];
OC = [211 215 222 224];
Con_constants = [232 233 234 235 236];
len = size(config);
feature_pos_r2(1) = 0;
for j=2:len(2)
    if ismember(config(j),SSTC) || ismember(config(j-1),SSTC)
        feature_pos_r2(j) = feature_pos_r2(j-1) + sstc_rad;
    elseif ismember(config(j),OC) || ismember(config(j-1),OC)
        feature_pos_r2(j) = feature_pos_r2(j-1) + oc_rad;
    elseif ismember(config(j),Con_constants)
        feature_pos_r2(j) = feature_pos_r2(j-1);
    else
        error('Mistake in setting feature position for second row')
    end
end
end

function [stent_r2_config] = Stent_r2_builder(cnos2)
%Nomeclature
% Short Strut SS --> BC2(211) - SSD2(212) - SSTC (213) - BC1(215)
% Tall Strut --> TSD2 (221) - TSTC(222) - TSD1(223) - BC3(224)
% Connector --> CS (231) [To Be Added later] - CC1 (232) - CS1(233) - CC
% (234) - CS2(235) CC2(236)

Nreps = cnos2/5; % 5 crowns in each base unit
SS = [211 212 213 214 215];
TS = [221 222 223 224];
C = [231 232 233 234 235 236];
base_unit = [C SS TS];
row_2 = [];
for i=1:Nreps
    row_2 = [row_2 base_unit]; 
end
stent_r2_config = row_2;
end
function [r2_SSD1_strut] = new_length_strut_builder(D1_strut, beta)
r2_SSD1_strut = [D1_strut(:,1) D1_strut(:,2)*beta];
end

function [crown_points] = crown_aggregator_V2(parsec_points_left,parsec_points_right, circular_points, type)
    sz = size(parsec_points_right);
    if sz(1) ~= 0
        parsec_points_right(1,:) = [];
    end
    sz2 = size(parsec_points_left);
    if sz(1) ~= 0
        parsec_points_left(1,:) = [];
    end
if type == "standard"
    crown_points_right = [circular_points;parsec_points_right];
    unflipped_crown_points_left = [circular_points;parsec_points_left];
    crown_points_left = flip(unflipped_crown_points_left);
    crown_points_left(:,2) = -1*crown_points_left(:,2);
    crown_points = [crown_points_left;crown_points_right];
end
if type == "ts_connector"
    parsec_points_left(1,:) = [];
    crown_points_right = [circular_points];
    crown_points_left = flip([circular_points;parsec_points_left]);
    crown_points_left(:,2) = -1*crown_points_left(:,2);
    crown_points = [crown_points_left;crown_points_right];
end
if type == "connector_cc1"
    parsec_points_left(1,:) = [];
    crown_points_left = circular_points;
    crown_points_left(:,2) = -1*crown_points_left(:,2);
    crown_points = [crown_points_left;parsec_points_left];
end
if type == "connector_cc2"
    crown_points_left = flip(circular_points);
    crown_points_left(:,2) = -1*crown_points_left(:,2);
    crown_points = crown_points_left;
end
end

function [crown_points] = generate_crown_parsec(A, type, np_crown, x_te, rc)
if type == "top"
    x = rc:((x_te-rc)/np_crown):x_te;
    x_mat = [x.^(0.5);x.^(1.5);x.^(2.5);x.^(3.5);x.^(4.5);x.^(5.5)];
    A_mat = A';
    z = A_mat*x_mat;
    %figure()
    %plot(z(1,:),-x);
    crown_points = [-x' z(1,:)'];
end

if type == "bottom"
    x = rc:((x_te-rc)/np_crown):x_te;
    x_mat = [x.^(0.5);x.^(1.5);x.^(2.5);x.^(3.5);x.^(4.5);x.^(5.5)];
    A_mat = A';
    z = A_mat*x_mat;
    %figure()
    %plot(z(1,:),-x);
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
slopes = (strut([2:np_strut],2) - strut([1:np_strut-1],2))./(strut([2:np_strut],1) - strut([1:np_strut-1],1));
end

function [row_1, wrap_name] = row_1_builder(stent_r1_config, feature_pos, D1_strut, D2_strut, crown_top, crown_bottom,strut_half_height,cr, stent_width, stent_thickness, cnos, crimped_cylinder_radius, RowAx)
len = size(stent_r1_config);
row_1 = [];
sketch_name_index = 1;
sketch_name = strings(0);
for j=1:len(2)
   if stent_r1_config(j) == 1
        feature_x = crown_top(:,2) + feature_pos(j);
        feature_y = crown_top(:,1) + strut_half_height + cr;
        sw_data_points = spline_data_points([feature_x feature_y]);
        py.sld_interface.create_new_sketch()% Unit feature starts here so create new sketch at Front Plane
        spline_id = py.sld_interface.create_spline(sw_data_points);
        first_point = [feature_x(1)/1000 feature_y(1)/1000];
        if j ~= 1
            py.sld_interface.create_line(first_point,last_point);% Create joining line
        else
            feature_x_prev = D2_strut(:,1) + feature_pos(j+3);
            feature_y_prev = D2_strut(:,2);
            last_point_prev = [feature_x_prev(end)/1000 feature_y_prev(end)/1000];
            last_point_prev_transformed(1) = last_point_prev(1) - (2*2*pi()*crimped_cylinder_radius/(cnos*1000));
            last_point_prev_transformed(2) = last_point_prev(2);
            py.sld_interface.create_line(first_point,last_point_prev_transformed);% Create joining line
        end

        last_point = [feature_x(end)/1000 feature_y(end)/1000];
        extra_last_point = [(feature_x(1) + 4*cr)/1000 feature_y(1)/1000];% Extra point to add at the end of the strut to ensure there is no 0 thickness combination in solidworks
        
   elseif stent_r1_config(j) == 3
       feature_x = crown_bottom(:,2) + feature_pos(j);
       feature_y = crown_bottom(:,1) - strut_half_height - cr;
       sw_data_points = spline_data_points([feature_x feature_y]);
       py.sld_interface.create_spline(sw_data_points);
       first_point = [feature_x(1)/1000 feature_y(1)/1000];
       py.sld_interface.create_line(first_point,last_point);% Create joining line
       last_point = [feature_x(end)/1000 feature_y(end)/1000];

   elseif stent_r1_config(j) == 2
        feature_x = D1_strut(:,1) + feature_pos(j);
        feature_y = D1_strut(:,2);
        sw_data_points = spline_data_points([feature_x feature_y]);
        py.sld_interface.create_spline(sw_data_points);   
        first_point = [feature_x(end)/1000 feature_y(end)/1000];
        py.sld_interface.create_line(first_point,last_point);% Create joining line
        last_point = [feature_x(1)/1000 feature_y(1)/1000];

   elseif stent_r1_config(j) == 4
        feature_x = D2_strut(:,1) + feature_pos(j);
        feature_y = D2_strut(:,2);       
        sw_data_points = spline_data_points([feature_x feature_y]);
        py.sld_interface.create_spline(sw_data_points);
        first_point = [feature_x(1)/1000 feature_y(1)/1000];
        py.sld_interface.create_line(first_point,last_point);% Create joining line
        last_point = [feature_x(end)/1000 feature_y(end)/1000];
        py.sld_interface.create_line(last_point, extra_last_point);% Create an extra joining line to avoid a 0 thickness geometry when combining.
        py.sld_interface.create_stent_segment_width(spline_id,stent_width/1000);%Create sketch width
        activesketch = py.sld_interface.get_active_sketch_name();
        sketch_name(sketch_name_index) = activesketch;
        sketch_name_index = sketch_name_index + 1;
        py.sld_interface.end_current_sketch();
   end
   feature = [feature_x feature_y];
   plot(RowAx, feature_x, feature_y,'LineWidth',1);
   hold(RowAx, "on");
   row_1 = [row_1 feature'];
end
    no_sketches = size(sketch_name,2);
    wrapped = 0;
    for k=1:1:no_sketches
        if wrapped ~= 0
            dummy_var = "dummy_var_2";
        else
            printing_sketc_name = sketch_name(k);
            [wrap_name] = py.sld_interface.wrap_sketch_around_base_cylinder_v2(sketch_name(k), stent_thickness/1000);
            if wrap_name ~= "NOWRAP"
                wrapped = 1;
            end            
        end
    end
row_1 = row_1';
end

function [straight_line] = connector_top_straight(length, np)
    x_pts = 0:(length/np):length;
    y_pts = x_pts*0;
    straight_line = [x_pts' y_pts'];
end

function [crown_top] = crown_top_circular(cr, np_crown)
x_mid_pts = 0:(cr/(np_crown)):cr;
y_mid_pts = sqrt(cr^2 - x_mid_pts.^2);
y_mid_pts = y_mid_pts - cr;
crown_top = [y_mid_pts' x_mid_pts'];
end

function [feature_pos] = feature_position(config,cr)
len = size(config);
for j=1:len(2)
    if j == 1
        feature_pos(j) = cr;
    else
        feature_pos(j) = feature_pos(j-1) + cr;
    end
end
end

function [stent_r1_config] = Stent_r1_builder(nc)
%Nomeclature
% crown top - 1
% Strut direction 1 - 2
% crown bottom - 3
% strut direction 2 - 4
base_unit = [1 2 3 4];
row_1 = [];
for i=1:nc
    row_1 = [row_1 base_unit]; 
end
stent_r1_config = row_1;
end

function [B] = Simple_NURBS(ww,P,t,k,n,np)
% This sends ww which is the vector of variables sought
w = [ww(1) ww(2) ww(3) ww(4) ww(4) ww(3) ww(2) ww(1) ];
u=linspace(min(t),max(t),np);
N = cell(n + k,k);
for i=1:k
    if i==1
        for j=1:n+k
            for l=1:length(u)
                if u(l)>t(j) && u(l)<t(j+1)
                    N{j,i}(l)=1;
                else
                    N{j,i}(l)=0;
                end
            end
        end
    else
        for j=1:(length(N(:,i-1))-(i-1))

            for l=1:length(u)
                if (t(j+i-1)-t(j))==0
                    temp1=0;
                else
                    temp1=((u(l)-t(j))*N{j,i-1}(l))/(t(j+i-1)-t(j));
                end
                if (t(j+i)-t(j+1))==0
                    temp2=0;
                else
                    temp2=((t(j+i)-u(l))*N{j+1,i-1}(l))/(t(j+i)-t(j+1));
                end
                N{j,i}(l)=temp1+temp2;
            end
        end
    end
end
R = cell(n + 1,1);
denom = zeros(1,length(u));
for i = 1:n+1
    denom = w(i) * N{i,k} + denom;
end
for i = 1:n+1
    num = w(i) * N{i,k};
    R{i,1} = num./denom;
end
B = zeros(2,length(u)); % Defing a matrix to store all points
for i = 1:n+1
    B(1,:) = P(1,i) * R{i,1} + B(1,:);
end
for j = 1:n+1
    B(2,:) = P(2,j) * R{j,1} + B(2,:);
end
% To add first and last control points, to get continous curve
B(:,1) = P(:,1);
B(:,end) = P(:,end);
B=B';
end

% Solidworks Helper Functions 
function [data_points] = spline_data_points(spline_pts) % Convert points in [x1 y1;x2 y2] format to [x1 yz1 z1 x2 y2 z2] format
    spline_pts(:,3) = 0;
    spline_pts_trans = spline_pts';
    data_points = reshape(spline_pts_trans, 1, []);
    data_points = data_points/1000;
end