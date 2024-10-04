function HS1_Generation(x, resolution_threshold, result_location)
%% Start Solidworks
py.sld_interface.start();

% Associating inputs with variable names
NS = x(1); % Number of unit struts in a helix x(1)
crimped_radius = x(2);% x(2)
alpha = x(3); %Helix angle x(3)
k = x(4);% Approximate number of units in a  % x(4)
SL = x(5);% Stent length x(5)
L1 = x(6); %End ROw height (6)
wc = x(7); %x(7) Width of the S Connector
strut_width = x(19);
strut_thickness = x(20);

% Other independent variables - Weights for different NURBS
ww_ab(1) = x(8);
ww_ab(2) = x(9);
ww_ab(3) = x(10);
ww_ab(4) = x(11);
ww_da(1) = x(12);
ww_da(2) = x(13);
ww_da(3) = x(14);
ww_con(1) = x(15);
ww_con(2) = x(16); 
ww_con(3) = x(17);
ww_con(4) = x(18);

% Associating repair variables with names
cp_width_ab = x(21);% 0.3 was default x(17)
cp_width_da = x(22); % x(18)
cp_height_da = x(23);% x(19)
% Inputs end here

%Fixed parameters 
no_helix = 2; % Number of helixes in the model
cut_off_distance_for_struts = 0.05; % Cut off location for parsecs
np_strut = 500;% Number of points in struts
np_crown = 200;% Number of points in crown
distance_thresh = resolution_threshold(1); % Distance threshold for strut neighbourhood points 
distance_thresh_conn = resolution_threshold(2);% Distance threshold for connectors
distance_thresh_er = resolution_threshold(3);% Distance threshold for end row struts
distance_thresh_tr = distance_thresh_er*10;% Distance threshold for transition row 
extra_point_thresh = 1E-3;% Distance at which an extra point is added at the end of strut wraps to allow combine operation in solidworks
extrema_width_factor = 1.6/1.1;% Extrema width factor applied to get the extrema width from the strut width for the connector

% Calculating dependent variables
NS2 = (2*NS + 1)/2;
HH = SL - 2*L1;
R = crimped_radius;
WDAN = cp_width_da;
WAB = cp_width_ab;
HDAN = cp_height_da;
crx = (1/(4*k))*(2*pi()*R - (k - 1)*WDAN -2*wc - 2*k*WAB);
cr_da = crx;
AAX = 2*WAB + 4*crx + WDAN;
HAB = -1*NS2*AAX*tand(alpha) + (SL - 2*L1);
cp_height_ab = HAB;
cry = 0.25*(2*HAB -2*HDAN -AAX*tand(alpha));
cr_yscale = cry/crx;
AAY = tand(alpha)*AAX;
hc = (k*AAY/2); % Connector height - constrained due to C7 - Helical ring alignment vertical constraint. 
stent_length = SL; % From now on SL is stent_length because I mistakenly reused SL symbol for an another variable
end_row_height = L1;

%% Generating Main Unit Strut
% Generating AB strut
end_bottom_ab = [-1.*cp_width_ab/2 -1.*cp_height_ab/2];
end_top_ab = [cp_width_ab/2 cp_height_ab/2];
CP_ab=[end_bottom_ab ;end_top_ab(1)/4 end_bottom_ab(2); end_top_ab(1)/4 end_bottom_ab(2)/2; end_top_ab(1)/4 0; end_bottom_ab(1)/4 0; end_bottom_ab(1)/4 end_top_ab(2)/2; end_bottom_ab(1)/4 end_top_ab(2); end_top_ab]; % Control points
P_ab=[CP_ab(:,1)';CP_ab(:,2)'];n_ab = size(P_ab,2) - 1;k_ab = 7;t_ab = [0 0 0 0 0 0 0 1 2 3 3 3 3 3 3 3];
[B_ab] = NURBS_AB(ww_ab,P_ab,t_ab,k_ab,n_ab, np_strut);
strut_ab_untransformed = flip(B_ab,1); % Still needs to be transformed so that it starts from 0,0
strut_ab(:,1) = strut_ab_untransformed(:,1) - cp_width_ab/2;
strut_ab(:,2) = strut_ab_untransformed(:,2) - cp_height_ab/2;
% Code 1 - Strut AB
strut_ab = removing_close_neibours(strut_ab, distance_thresh); 

% Code 3 - Strut CD. 
% C-D strut is same as AB but with x and y shifts
strut_cd = strut_ab;

% Generating DA strut
WDA = WDAN + 2*crx;
HDA = HDAN +2*cry;
end_bottom_da = [cp_width_da/2 -1*cp_height_da/2];
end_top_da = [-1*cp_width_da/2 cp_height_da/2];
CP_da=[end_bottom_da ;end_bottom_da(1) end_bottom_da(2)/2; end_bottom_da(1) 0; end_top_da(1) 0; end_top_da(1) end_top_da(2)/2; end_top_da]; % Control points
P_da=[CP_da(:,1)';CP_da(:,2)'];n_da = size(P_da,2) - 1;k_da = 5;t_da = [0 0 0 0 0 1 2 3 3 3 3 3];
[B_da] = NURBS_DA(ww_da,P_da,t_da,k_da,n_da, np_strut);

% Generating crown points for DA Strut
crown_points_circular_da = crown_top_circular(cr_da, 4*np_crown);

% Generating DA strut 
strut_da_untransformed = strut_da_generator(B_da, crown_points_circular_da, cr_da, cp_width_da, cp_height_da, cr_yscale);% Still needs to be transformed so that it starts from 0,0
strut_da(:,1) = strut_da_untransformed(:,1) - cr_da - cp_width_da/2;
strut_da(:,2) = strut_da_untransformed(:,2) + cr_da*cr_yscale + cp_height_da/2;

% Code 4 - Strut DA
strut_da = removing_close_neibours(strut_da, distance_thresh);

% Generating BC Strut 
% Modifying DA  to BC Strut - BC strut is basically same as DA strut but
% without the width in the NURBS section. The linear correction applied to
% DA strut thus requires an additional parsec to maintain continuity and
% differentiability
HBC = HDA;
WBC = 2*crx;
strut_half_height = cp_height_da/2;
slope_correction_cutoff_value = strut_half_height - cut_off_distance_for_struts;
B_bc = DA_to_BC(B_da, cp_width_da, cp_height_da, np_strut);
D1_strut = B_bc;
D1_slopes = strut_slopes(D1_strut, np_strut);
[D1_corrected_strut, D1_corrected_slope] = strut_slope_correction(D1_strut, D1_slopes, slope_correction_cutoff_value);

% generating the parsec points for top crowns 
z_te = D1_corrected_strut(end,1) + cr_da;
x_te = -(D1_corrected_strut(end,2) - strut_half_height - cr_da);
z_te_slope = -1/D1_corrected_slope(end,1);
parsec_coeffs_top_A = parsec(cr_da, x_te, z_te, z_te_slope);
crown_points_parsec_top = generate_crown_parsec(parsec_coeffs_top_A, "top", np_crown/4,x_te,cr_da, cr_yscale);

% generationg the semicircular crown points for top
crown_points_circular = crown_top_circular(cr_da, 4*np_crown);
crown_top = crown_aggregator_V3(crown_points_parsec_top, crown_points_parsec_top, crown_points_circular,'standard', cr_yscale);

% generating the parsec points for bottom crown
z_te_b = -D1_corrected_strut(1,1) + cr_da;
x_te_b = (D1_corrected_strut(1,2) + strut_half_height + cr_da);
z_te_slope_b = -1/D1_corrected_slope(1,1);
parsec_coeffs_bottom = parsec(cr_da, x_te_b, z_te_b, z_te_slope_b);
crown_points_parsec_bottom = generate_crown_parsec(parsec_coeffs_bottom, "bottom", np_crown/4,x_te_b,cr_da, cr_yscale);
crown_points_circular_bottom = crown_top_circular(cr_da, 4*np_crown);
crown_bottom_rev = crown_aggregator_V3(crown_points_parsec_bottom, crown_points_parsec_bottom, crown_points_circular_bottom,'standard',cr_yscale);
crown_bottom = [-1.*crown_bottom_rev(:,1) -1.*crown_bottom_rev(:,2)];

% Generating BC strut 
strut_bc_untransformed = strut_bc_generator(D1_corrected_strut, crown_top, crown_bottom, cr_da, cp_height_da, cr_yscale); % Still needs to be transformed so that it starts from 0,0
strut_bc(:,1) = strut_bc_untransformed(:,1) - cr_da;
strut_bc(:,2) = strut_bc_untransformed(:,2) + cp_height_da/2 + cr_da*cr_yscale;

% Code 3 - Strut BC
strut_bc = removing_close_neibours(strut_bc, distance_thresh);

% Generating Main connector 
extrema_width = extrema_width_factor*wc;
extrema_half_width = extrema_width*0.5;
end_bottom_con = [wc/2 -hc/2];
end_top_con = [-wc/2 hc/2];
CP_con=[end_bottom_con ;(-1*extrema_half_width) end_bottom_con(2); -1*extrema_half_width end_bottom_con(2)/2; -1*extrema_half_width 0; extrema_half_width 0; extrema_half_width end_top_con(2)/2; extrema_half_width end_top_con(2); end_top_con]; % Control points
P_con=[CP_con(:,1)';CP_con(:,2)'];n_con = size(P_con,2) - 1;k_con = 7;t_con = [0 0 0 0 0 0 0 1 2 3 3 3 3 3 3 3];
[B_con] = NURBS_AB(ww_con,P_con,t_con,k_con,n_con, np_strut);% The simple NUrbs formulation for AB and connector is the same, just the parameters entred are different
connector_untransformed = flip(B_con,1); % Flip so that we start from top. Still needs to be transformed so that it starts from 0,0
connector(:,1) = connector_untransformed(:,1) + wc/2;
connector(:,2) = connector_untransformed(:,2) - hc/2;

% Code 5 - connector
connector = removing_close_neibours(connector, distance_thresh_conn);

% Generating unit strut
unit_strut = generate_unit_strut(strut_ab, strut_bc, strut_cd, strut_da, connector, cp_width_ab, cp_height_ab, cp_width_da, cp_height_da, cr_da, cr_yscale);

%% Generating Transition row Strut
% AR strut segment Psrameters
cr_ar = (1/5)*(((pi()*R)/no_helix) - WAB);% Radius of Crowns in the AR section of transition strut 

% AP Strut segment
psi = 0.84/0.95;
HAP = psi*end_row_height;% Height of straight portion of AP strut
HAPS = HAP - 2*cr_ar;
S = straight_strut_generator(HAPS, np_strut);
crown_points_circular_ap = crown_top_circular(cr_ar, 4*np_crown);
strut_ap_untransformed = strut_ap_generator(S, crown_points_circular_ap, cr_ar, HAPS);
WAP = 2*cr_ar;
strut_ap(:,1) = strut_ap_untransformed(:,1) + cr_ar;
strut_ap(:,2) = strut_ap_untransformed(:,2) - cr_ar - (HAPS/2);

% Code 21 - Strut AP
strut_ap = removing_close_neibours(strut_ap, distance_thresh_tr);

% PQ strut segment
alpha_ratio = 0.91/0.95;
HPQ = alpha_ratio*end_row_height;
HPQS = HPQ - 2*cr_ar;% Height of the straight portion of PQ
WPQS = cr_ar; % Width of the straight portion of PQ
WPQ = 2*cr_ar-WPQS;
SL = slanted_strut_generator(HPQS, WPQS, np_strut);
D1PQ_strut(:,1) = -1*SL(:,1);
D1PQ_strut(:,2) = SL(:,2);
D1PQ_slopes = strut_slopes(D1PQ_strut, np_strut);
slope_correction_cutoff_valuePQ = 0.5*HPQS - cut_off_distance_for_struts;
[D1PQ_corrected_strut, D1PQ_corrected_slope] = strut_slope_correction(D1PQ_strut, D1PQ_slopes, slope_correction_cutoff_valuePQ);

% generating the parsec points for top crowns 
z_tePQ = D1PQ_corrected_strut(end,1) + cr_ar - (WPQS/2);
x_tePQ = -(D1PQ_corrected_strut(end,2) - 0.5*HPQS - cr_ar);
z_te_slopePQ = -1/D1PQ_corrected_slope(end,1);
parsec_coeffs_top_APQ = parsec(cr_ar, x_tePQ, z_tePQ, z_te_slopePQ);
crown_points_parsec_topPQ = generate_crown_parsec(parsec_coeffs_top_APQ, "top", np_crown/4,x_tePQ,cr_ar, 1);

% generationg the semicircular crown points for top
crown_points_circularPQ = crown_top_circular(cr_ar, 4*np_crown);
crown_topPQ = crown_aggregator_V3(crown_points_parsec_topPQ, crown_points_parsec_topPQ, crown_points_circularPQ,'standard', 1);

% generating the parsec points for bottom crown
z_te_bPQ = -D1PQ_corrected_strut(1,1) + cr_ar - (WPQS/2);
x_te_bPQ = (D1PQ_corrected_strut(1,2) + 0.5*HPQS + cr_ar);
z_te_slope_bPQ = -1/D1PQ_corrected_slope(1,1);
parsec_coeffs_bottomPQ = parsec(cr_ar, x_te_bPQ, z_te_bPQ, z_te_slope_bPQ);
crown_points_parsec_bottomPQ = generate_crown_parsec(parsec_coeffs_bottomPQ, "bottom", np_crown/4,x_te_bPQ,cr_ar, 1);
crown_points_circular_bottomPQ = crown_top_circular(cr_ar, 4*np_crown);
crown_bottom_revPQ = crown_aggregator_V3(crown_points_parsec_bottomPQ, crown_points_parsec_bottomPQ, crown_points_circular_bottomPQ,'standard',1);
crown_bottomPQ = [-1.*crown_bottom_revPQ(:,1) -1.*crown_bottom_revPQ(:,2)];
strut_pq_untransformed = strut_pq_generator(D1PQ_corrected_strut, crown_topPQ, crown_bottomPQ, cr_ar, HPQS, WPQS);
strut_pq(:,1) = strut_pq_untransformed(:,1) + cr_ar - (WPQS/2);
strut_pq(:,2) = strut_pq_untransformed(:,2) + cr_ar + (HPQS/2);

% Code 22 - Strut PQ
strut_pq = removing_close_neibours(strut_pq, distance_thresh_tr);

% QR Strut segment
beta = 0.62/0.95;
HQRS = beta*end_row_height - 2*cr_ar;% Height of straight portion of QR strut
HQR = HQRS + 2*cr_ar;
WQR = 2*cr_ar;
SQR = straight_strut_generator(HQRS, np_strut);
crown_points_circular_qr = crown_top_circular(cr_ar, 4*np_crown);
strut_qr_untransformed = strut_ap_generator(SQR, crown_points_circular_qr, cr_ar, HQRS); % QR generator is same as ap generator with different input parameters
strut_qr(:,1) = strut_qr_untransformed(:,1) + cr_ar;
strut_qr(:,2) = strut_qr_untransformed(:,2) - cr_ar - (HQRS/2);

% code 23
strut_qr = removing_close_neibours(strut_qr, distance_thresh_tr);

% RC Strut segment
sigma = 9/7; % Proportion of width w.r.t cr_rc covered by TU middle portion
delta = 5/7;% Poprtion of the width w.r.t. cr_rc covered by RS Nurbs
cr_rc = (1/(7 + sigma + delta))*((pi()*R/no_helix) - 2*crx);

% RS Strut
HRSM = beta*end_row_height - 2*cr_rc;
omega1 = 0.09/0.48;% Portion of height of RSM covered by the bottom straight part
omega2 = 0.12/0.48;% Portion of height of RSM covered by the middle slanted part
omega3 = 0.27/0.48;% Portion of height of RSM covered by the upper straight part
HRSM1 = omega1*HRSM;% Height of the bottom straight part of HRSM
HRSM2 = omega2*HRSM;% Height of the middle slanted part of HRSM
HRSM3 = omega3*HRSM;% Height of the top straight part of HRSM
WRSM2 = delta*cr_rc;
SRSM1 = straight_strut_generator(HRSM1, round(np_strut/3));
SRSM2 = slanted_strut_generator(HRSM2, WRSM2, round(np_strut/3));
SRSM2(:,1) = -1*SRSM2(:,1);
SRSM3 = straight_strut_generator(HRSM3, round(np_strut/3));
crown_points_circular_rs = crown_top_circular(cr_rc, 4*np_crown);
D1SRSM2_strut = SRSM2;
D1SRSM2_slopes = strut_slopes(D1SRSM2_strut, round(np_strut/3));
slope_correction_cutoff_valueSRSM2 = 0.5*HRSM2 - (cut_off_distance_for_struts/2);
[D1SRSM2_corrected_strut, D1SRSM2_corrected_slope] = strut_slope_correction(D1SRSM2_strut, D1SRSM2_slopes, slope_correction_cutoff_valueSRSM2);

% generating the parsec points for top straight line 
z_teSRSM2 = D1SRSM2_corrected_strut(end,1);
x_teSRSM2 = -(D1SRSM2_corrected_strut(end,2) - 0.5*HRSM2 - 0.5*WRSM2);
z_te_slopeSRSM2 = -1/D1SRSM2_corrected_slope(end,1);
parsec_coeffs_top_ASRSM2 = parsec(0.5*WRSM2, x_teSRSM2, z_teSRSM2, z_te_slopeSRSM2);
crown_points_parsec_topSRSM2 = generate_crown_parsec(parsec_coeffs_top_ASRSM2, "top", np_crown/4,x_teSRSM2,0.5*WRSM2, 1);

% generationg the semicircular crown points for top
parsec_topSRSM2 = crown_points_parsec_topSRSM2;

% generating the parsec points for bottom intersection
z_te_bSRSM2 = -D1SRSM2_corrected_strut(1,1);
x_te_bSRSM2 = (D1SRSM2_corrected_strut(1,2) + 0.5*HRSM2 + 0.5*WRSM2);
z_te_slope_bSRSM2 = -1/D1SRSM2_corrected_slope(1,1);
parsec_coeffs_bottomSRSM2 = parsec(0.5*WRSM2, x_te_bSRSM2, z_te_bSRSM2, z_te_slope_bSRSM2);
crown_points_parsec_bottomSSRM2 = generate_crown_parsec(parsec_coeffs_bottomSRSM2, "bottom", np_crown/4,x_te_bSRSM2,0.5*WRSM2, 1);
parsec_bottom_revSRSM2 = crown_points_parsec_bottomSSRM2;
parsec_bottomSRSM2 = [-1.*parsec_bottom_revSRSM2(:,1) -1.*parsec_bottom_revSRSM2(:,2)];
RSM2_untransformed = strut_pq_generator(D1SRSM2_corrected_strut, parsec_topSRSM2, parsec_bottomSRSM2, 0.5*WRSM2, HRSM2, WRSM2);
RSM2_untransformed(:,1) = -1.*RSM2_untransformed(:,1);% Because the strut_pq_generator generates the mirror image along x, we have to take the mirror image again. 
HRS = HRSM+2*cr_rc;
WRS = WRSM2 + 2*cr_rc;
strut_rs = rs_generator(SRSM1, RSM2_untransformed, SRSM3, crown_points_circular_rs, HRSM1, HRSM2, HRSM3, WRSM2, cr_rc);

% code 24
strut_rs = removing_close_neibours(strut_rs, distance_thresh_tr);

% ST Strut segment
gamma = 0.58/0.95;
HSTS = gamma*end_row_height - 2*cr_rc;
crown_points_circular_st = crown_top_circular(cr_rc, 4*np_crown);
SST = straight_strut_generator(HSTS, np_strut);
strut_st_untransformed = strut_ap_generator(SST, crown_points_circular_st, cr_rc, HSTS); % Strut AP generator can be used to create ST strut as they are both similar with different dimentions
HST = 2*cr_rc + HSTS;
WST = 2*cr_rc;
strut_st(:,1) = strut_st_untransformed(:,1) + cr_rc;
strut_st(:,2) = strut_st_untransformed(:,2) - cr_rc - (HSTS/2);

% code 25
strut_st = removing_close_neibours(strut_st, distance_thresh_tr);

% Generating TU Strut
% TU Strut is shaped like a mirror image of DA NURBS. THerefore the method to
% generate TU Strut would be to produce DA NURBS with modified dimensions
% and then take a mirror image
sigma = 9/7;
HTUN = HSTS;
WTUN = sigma*cr_rc;
end_bottom_tu = [WTUN/2 -1*HTUN/2];
end_top_tu = [-1*WTUN/2 HTUN/2];
CP_tu=[end_bottom_tu ;end_bottom_tu(1) end_bottom_tu(2)/2; end_bottom_tu(1) 0; end_top_tu(1) 0; end_top_tu(1) end_top_tu(2)/2; end_top_tu]; % Control points
P_tu=[CP_tu(:,1)';CP_tu(:,2)'];n_tu = size(P_tu,2) - 1;k_tu = 5;t_tu = [0 0 0 0 0 1 2 3 3 3 3 3];
[B_tu] = NURBS_DA(ww_da,P_tu,t_tu,k_tu,n_tu, np_strut);% Using NURBS_DA to generate TU strut. also using same weights as DA

% Generating crown points for TU Strut
crown_points_circular_tu = crown_top_circular(cr_rc, 4*np_crown);

% Generating TU strut 
strut_tu_untransformed = strut_da_generator(B_tu, crown_points_circular_tu, cr_rc, WTUN, HTUN, 1);% Still needs to be transformed so that it starts from 0,0
HTU = 2*cr_rc + HTUN;
WTU = 2*cr_rc + WTUN;
strut_tu(:,1) = strut_tu_untransformed(:,1) - cr_rc - 0.5*WTUN;
strut_tu(:,2) = strut_tu_untransformed(:,2) + cr_rc + 0.5*HTUN;
strut_tu(:,1) = -1.*strut_tu(:,1);

% Code 26
strut_tu = removing_close_neibours(strut_tu, distance_thresh_tr);

% Generating strut UC
HUCS = HAB + (alpha_ratio - psi)*end_row_height - 2*cry - HDAN - cr_rc;
crown_points_circular_uc = crown_top_circular(cr_rc, 4*np_crown);
SUC = straight_strut_generator(HUCS, np_strut);
HUC = HUCS + cr_rc;
WUC = cr_rc;
strut_uc = strut_uc_generator(SUC, crown_points_circular_st, cr_rc, HUCS); % Strut AP generator can be used to create ST strut as they are both similar with different dimentions

% code 27
strut_uc = removing_close_neibours(strut_uc, distance_thresh_tr);

% Generating Transition strut
transition_unit_strut = generate_transition_unit_strut(strut_ap, strut_pq, strut_qr, strut_rs, strut_st, strut_tu, strut_uc, HAP,WAP, HPQ, WPQ, HQR, WQR, HRS, WRS, HST, WST, HTU, WTU);

%% Generating End Row Strut requirements
phi = 0.71/0.95; % Proportion of end row height taken by End row struts
N_er_crowns = 16;
cr_er = pi()*crimped_radius/N_er_crowns;
HERS = phi*end_row_height - 2*cr_er;% Height of straight part of end row strut 
SER = straight_strut_generator(HERS, np_strut);
crown_points_circular_er = crown_top_circular_v2(cr_er, 4*np_crown);
[crown_er_top, crown_er_bottom] = crown_er_generator(crown_points_circular_er, HERS, cr_er);
crown_er_bottom = removing_close_neibours(crown_er_bottom, distance_thresh_er);
crown_er_top = removing_close_neibours(crown_er_top, distance_thresh_er);
strut_er = SER;
strut_er(:,2) = strut_er(:,2) + 0.5*HERS + cr_er;
stent_er_config = er_config_builder(N_er_crowns); % End row config builder
feature_pos_er = feature_position_er(stent_er_config,cr_er);

%% Generating End Row Connectors
TCH = (1-phi)*end_row_height + HAB - cry -0.5*HDAN; %Height of tall connector
cr_tc = crx;
TCHS = TCH - cr_tc; % Height of straight portion of tall connector
STC = straight_strut_generator(TCHS,np_strut);
crown_points_circular_tc = crown_top_circular(cr_tc, 4*np_crown);
tall_connector = end_connector_generator(STC, crown_points_circular_tc, TCHS, cr_tc);
tall_connector = removing_close_neibours(tall_connector, distance_thresh_conn);
SCH = (1-phi)*end_row_height + HAP - HPQ + 0.5*HQR; % Height of short connector
cr_sc = -WAB + (2*pi*crimped_radius/(2*no_helix)) - (WAP + WPQ + 0.5*WQR);% Reference point is A which was taken to be 0,0
SCHS = SCH - cr_sc; % height of straight portion of short connector 
SSC = straight_strut_generator(SCHS,np_strut);
crown_points_circular_sc = crown_top_circular(cr_sc, 4*np_crown);
short_connector = end_connector_generator(SSC, crown_points_circular_sc, SCHS, cr_sc);
short_connector = removing_close_neibours(short_connector, distance_thresh_conn);

%% Building helical stent 
stent_helix_config = helix_config_builder(NS);
feature_pos_helix = feature_position(stent_helix_config,AAX, AAY);
RowFig = figure();
RowAx = axes('Parent', RowFig);
% Building helix including connectors between helix
[helix, sw_helix_name] = helix_builder(strut_ab, strut_bc, strut_cd, strut_da, connector, cp_width_ab, cp_height_ab, cp_width_da, cp_height_da, cr_da, stent_helix_config, feature_pos_helix, crimped_radius, no_helix, cr_yscale, NS, k, stent_length, strut_thickness, strut_width, extra_point_thresh, AAX, AAY, WDA, HDA, RowAx);

% Building tansition row struts
[transition_row, transition_row_rev, sw_helix_name] = transition_row_builder(strut_ap, strut_pq, strut_qr, strut_rs, strut_st, strut_tu, strut_uc, HAP, WAP, HPQ, WPQ, HQR, WQR, HRS, WRS, HST, WST, HTU, WTU, crimped_radius, no_helix, AAX, AAY, HAB, WAB, crx, cry, HDAN, NS, sw_helix_name, WBC, HBC, extra_point_thresh, strut_width, strut_thickness, stent_length, RowAx);

% Building connectors between end rows and main stent
[end_connector, end_connector_rev, sw_helix_name] = end_connector_builder(tall_connector, short_connector, WAB, HAB, cry, HDAN, no_helix, crimped_radius, cr_tc, WAP, WPQ, WQR, HAP, HPQ, HQR, AAX, AAY, NS, crx, sw_helix_name, stent_length, strut_width, strut_thickness, RowAx);

% Building end rows
[er, er_rev, sw_helix_name] = er_row_builder(crown_er_bottom, crown_er_top, strut_er, stent_er_config, feature_pos_er, end_row_height, phi, crimped_radius, no_helix, cr_er, AAX, NS, WAB, crx, AAY, HAB, cry, HDAN, sw_helix_name, extra_point_thresh, N_er_crowns, strut_width, strut_thickness, stent_length, RowAx);

%% Saving figure, stent solidworks file, parasolid file and Mesh STL file 

% Saving 2D image
RowAx.XLim = [-25 5];
RowAx.XTick = -25:1:5;
RowAx.YLim = [-20,2];
RowAx.YTick = -20:2:2;
set(RowFig, 'Position', [100, 100, 1400*1.5, 2200*1.5]);
close(RowFig);
close all;

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

%% Required helper functions

function [end_connector, end_connector_rev, sw_helix_name] = end_connector_builder(tall_connector, short_connector, WAB, HAB, cry, HDAN, no_helix, crimped_radius, cr_tc, WAP, WPQ, WQR, HAP, HPQ, HQR, AAX, AAY, NS, crx, sw_helix_name, stent_length, stent_width, stent_thickness, RowAx)
featurex = [];
featurey = [];
helix_shift_x = 0;
end_conn_iter = [];
end_connector = [];
featurex_rev = [];
featurey_rev = [];
helix_shift_x = 0;
end_conn_iter_rev = [];
end_connector_rev = [];
unit_ec_strut_list = string();
sw_ec_strut_flag = 0;
sw_ec_rev_strut_flag = 0;
sw_ec_rev_strut_flag = 0;

for j = 1:1:no_helix % Creating the top connector
    % Adding tall connector 
    strut_x = tall_connector(:,1) - WAB - cr_tc + helix_shift_x;
    strut_y = tall_connector(:,2) - HAB + cry + 0.5*HDAN; 
    plot(strut_x, strut_y, 'color', 'black');hold on;
    featurex = [featurex; strut_x];
    featurey = [featurey; strut_y];
    end_conn_iter = [end_conn_iter; j*ones(size(strut_x))];
    if sw_ec_strut_flag == 1
        ec_unit_new = py.sld_interface.copy_move(unit_ec_strut_list(1), "Axis1", 0,pi());
        unit_ec_strut_list(end+1) = ec_unit_new;
    elseif sw_ec_strut_flag == 0
        cylinder_name = py.sld_interface.create_base_cylinder(crimped_radius/1000, stent_length/1000);
        ec_strut = [strut_x strut_y];
        sw_data_points = spline_data_points(ec_strut);
        py.sld_interface.create_new_sketch();% create a new sketch at front plane
        spline_id = py.sld_interface.create_spline(sw_data_points);
        py.sld_interface.create_stent_segment_width(spline_id,stent_width/1000);%Create sketch width
        activesketch = py.sld_interface.get_active_sketch_name();
        sketch_name = activesketch;
        py.sld_interface.end_current_sketch();% Unit feature ends here so exit previous sketch        
        wrap_name = py.sld_interface.wrap_sketch_around_base_cylinder_v3(sketch_name, stent_thickness/1000);
        if wrap_name == "NOWRAP"
            error("Could not generate tall end connector wrap")
        end
        cylinder_name_2 = py.sld_interface.create_base_cylinder(crimped_radius/1000, stent_length/1000);
        ec_unit = py.sld_interface.combine_subtract(cylinder_name_2,wrap_name);
        unit_ec_strut_list(1) = ec_unit;
    end
    % Adding Short Connector
    strut_x = short_connector(:,1) + WAP + WPQ + 0.5*WQR + helix_shift_x;
    strut_y = short_connector(:,2) - HAP+ HPQ - 0.5*HQR;
    plot(strut_x, strut_y, 'color', 'red');hold on;
    featurex = [featurex; strut_x];
    featurey = [featurey; strut_y];
    end_conn_iter = [end_conn_iter; j*ones(size(strut_x))];
    if sw_ec_strut_flag == 1
        ec_unit_new = py.sld_interface.copy_move(unit_ec_strut_list(2), "Axis1", 0,pi());
        unit_ec_strut_list(end+1) = ec_unit_new;
    elseif sw_ec_strut_flag == 0
        cylinder_name = py.sld_interface.create_base_cylinder(crimped_radius/1000, stent_length/1000);
        ec_strut = [strut_x strut_y];
        sw_data_points = spline_data_points(ec_strut);
        py.sld_interface.create_new_sketch();% create a new sketch at front plane
        spline_id = py.sld_interface.create_spline(sw_data_points);
        py.sld_interface.create_stent_segment_width(spline_id,stent_width/1000);%Create sketch width
        activesketch = py.sld_interface.get_active_sketch_name();
        sketch_name = activesketch;
        py.sld_interface.end_current_sketch();% Unit feature ends here so exit previous sketch        
        wrap_name = py.sld_interface.wrap_sketch_around_base_cylinder_v3(sketch_name, stent_thickness/1000);
        if wrap_name == "NOWRAP"
            error("Could not generate short end connector wrap")
        end
        cylinder_name_2 = py.sld_interface.create_base_cylinder(crimped_radius/1000, stent_length/1000);
        ec_unit = py.sld_interface.combine_subtract(cylinder_name_2,wrap_name);
        unit_ec_strut_list(end+1) = ec_unit;
        sw_ec_strut_flag = 1;
    end
angle_shift = 2*pi()/no_helix;
helix_shift_x = helix_shift_x - angle_shift*crimped_radius;
end

helix = sw_helix_name;

for k = 1:1:size(unit_ec_strut_list,2)
    helix = py.sld_interface.row_combine(helix, unit_ec_strut_list(k));% I created this function to combine rows in IR1 and reusing it. That is why the name of the function is row_combine. It just combines solid bodies
end

unit_ec_rev_strut_list = string();
helix_shift_x = 0;

for j = 1:1:no_helix % Creating the bottom connector
    % Adding tall connector 
    strut_x = tall_connector(:,1) - WAB - cr_tc + helix_shift_x;
    strut_y = tall_connector(:,2) - HAB + cry + 0.5*HDAN; 
    % Generating the reverse strut
    strut_x_rev = -1*strut_x - AAX*NS - 2*WAB - 2*crx + 2*helix_shift_x;
    strut_y_rev = -1*strut_y - AAY*NS - 2*HAB + 2*cry + HDAN;
    plot(strut_x_rev, strut_y_rev,'color','black');hold on;
    end_conn_iter_rev = [end_conn_iter_rev;j*ones(size(strut_x_rev))];
    featurex_rev = [featurex_rev; strut_x_rev];
    featurey_rev = [featurey_rev; strut_y_rev];
    if sw_ec_rev_strut_flag == 1
        ec_unit_new = py.sld_interface.copy_move(unit_ec_rev_strut_list(1), "Axis1", 0,pi());
        unit_ec_rev_strut_list(end+1) = ec_unit_new;
    elseif sw_ec_rev_strut_flag == 0
        cylinder_name = py.sld_interface.create_base_cylinder(crimped_radius/1000, stent_length/1000);
        strut_x_rev_initial = strut_x_rev - (-AAX*NS - 2*WAB - 2*crx + 2*helix_shift_x);
        strut_y_rev_initial = strut_y_rev - (- AAY*NS - 2*HAB + 2*cry + HDAN);
        strut_sp_initial = [strut_x_rev_initial strut_y_rev_initial];
        sw_data_points = spline_data_points(strut_sp_initial);
        py.sld_interface.create_new_sketch();% create a new sketch at front plane
        spline_id = py.sld_interface.create_spline(sw_data_points);
        py.sld_interface.create_stent_segment_width(spline_id,stent_width/1000);%Create sketch width
        activesketch = py.sld_interface.get_active_sketch_name();
        sketch_name = activesketch;
        py.sld_interface.end_current_sketch();% Unit feature ends here so exit previous sketch        
        wrap_name = py.sld_interface.wrap_sketch_around_base_cylinder_v3(sketch_name, stent_thickness/1000);
        if wrap_name == "NOWRAP"
            error("Could not generate tall bottom end connector wrap")
        end
        cylinder_name_2 = py.sld_interface.create_base_cylinder(crimped_radius/1000, stent_length/1000);
        ec_unit = py.sld_interface.combine_subtract(cylinder_name_2,wrap_name);
        y_shift_3d_strut = - AAY*NS - 2*HAB + 2*cry + HDAN;
        theta_shift_3d_strut = (-AAX*NS - 2*WAB - 2*crx + 2*helix_shift_x)/crimped_radius;
        ec_unit = py.sld_interface.move(ec_unit, "Axis1", y_shift_3d_strut/1000,theta_shift_3d_strut);
        unit_ec_rev_strut_list(1) = ec_unit;
    end
    % Adding Short Connector
    strut_x = short_connector(:,1) + WAP + WPQ + 0.5*WQR + helix_shift_x;
    strut_y = short_connector(:,2) - HAP+ HPQ - 0.5*HQR;
    % Generating the reverse strut
    strut_x_rev = -1*strut_x - AAX*NS - 2*WAB - 2*crx + 2*helix_shift_x;
    strut_y_rev = -1*strut_y - AAY*NS - 2*HAB + 2*cry + HDAN;
    plot(strut_x_rev, strut_y_rev,'color','black');hold on;
    end_conn_iter_rev = [end_conn_iter_rev;j*ones(size(strut_x_rev))];
    featurex_rev = [featurex_rev; strut_x_rev];
    featurey_rev = [featurey_rev; strut_y_rev];
    if sw_ec_rev_strut_flag == 1
        ec_unit_new = py.sld_interface.copy_move(unit_ec_rev_strut_list(2), "Axis1", 0,pi());
        unit_ec_rev_strut_list(end+1) = ec_unit_new;
    elseif sw_ec_rev_strut_flag == 0
        cylinder_name = py.sld_interface.create_base_cylinder(crimped_radius/1000, stent_length/1000);
        strut_x_rev_initial = strut_x_rev - (-AAX*NS - 2*WAB - 2*crx + 2*helix_shift_x);
        strut_y_rev_initial = strut_y_rev - (- AAY*NS - 2*HAB + 2*cry + HDAN);
        strut_sp_initial = [strut_x_rev_initial strut_y_rev_initial];        
        sw_data_points = spline_data_points(strut_sp_initial);
        py.sld_interface.create_new_sketch();% create a new sketch at front plane
        spline_id = py.sld_interface.create_spline(sw_data_points);
        py.sld_interface.create_stent_segment_width(spline_id,stent_width/1000);%Create sketch width
        activesketch = py.sld_interface.get_active_sketch_name();
        sketch_name = activesketch;
        py.sld_interface.end_current_sketch();% Unit feature ends here so exit previous sketch        
        wrap_name = py.sld_interface.wrap_sketch_around_base_cylinder_v3(sketch_name, stent_thickness/1000);
        if wrap_name == "NOWRAP"
            error("Could not generate short bottom end connector wrap")
        end
        cylinder_name_2 = py.sld_interface.create_base_cylinder(crimped_radius/1000, stent_length/1000);
        ec_unit = py.sld_interface.combine_subtract(cylinder_name_2,wrap_name);
        y_shift_3d_strut = - AAY*NS - 2*HAB + 2*cry + HDAN;
        theta_shift_3d_strut = (-AAX*NS - 2*WAB - 2*crx + 2*helix_shift_x)/crimped_radius;
        ec_unit = py.sld_interface.move(ec_unit, "Axis1", y_shift_3d_strut/1000,theta_shift_3d_strut);
        unit_ec_rev_strut_list(end+1) = ec_unit;
        sw_ec_rev_strut_flag = 1;
    end    
angle_shift = 2*pi()/no_helix;
helix_shift_x = helix_shift_x - angle_shift*crimped_radius;
end

for k = 1:1:size(unit_ec_rev_strut_list,2)
    helix = py.sld_interface.row_combine(helix, unit_ec_rev_strut_list(k));% I created this function to combine rows in IR1 and reusing it. That is why the name of the function is row_combine. It just combines solid bodies
end

end_connector = [end_conn_iter featurex featurey];
sw_helix_name = helix;
end

function [er, er_rev, sw_helix_name] = er_row_builder(crown_er_bottom, crown_er_top, strut_er, stent_er_config, feature_pos_er, end_row_height, phi, crimped_radius, no_helix, cr_er, AAX, NS, WAB, crx, AAY, HAB, cry, HDAN, sw_helix_name, extra_point_thresh, N_er_crowns, stent_width, stent_thickness, stent_length,RowAx)
er = [];
er_rev = [];
featurex = [];
featurey = [];
featurex_rev = [];
featurey_rev = [];
sw_er_strut_flag = 0;% 0 means that the strut has ti be built. 1 means the base strut is built and only the copy function is left, 2 means all required struts are built.
sw_er_rev_strut_flag = 0;% Same meaning as above

for j = 1:1:size(stent_er_config,2)% Creating the top end row
    if stent_er_config(j) == 1
        strut_x = crown_er_bottom(:,1) - feature_pos_er(j) + 2*pi*crimped_radius/(no_helix) + cr_er - WAB;
        strut_y = crown_er_bottom(:,2) + (1-phi)*end_row_height;
        plot(RowAx, strut_x, strut_y,'LineWidth',1);hold on;
    if sw_er_strut_flag == 0
            cylinder_name = py.sld_interface.create_base_cylinder(crimped_radius/1000, stent_length/1000);        
            strut_sp_initial = [strut_x strut_y];
            sw_data_points = spline_data_points(strut_sp_initial);% Adding an extra point at the start to maintain 0 slope at start
            py.sld_interface.create_new_sketch();% Unit feature starts here so create new sketch at Front Plane
            spline_id = py.sld_interface.create_spline(sw_data_points);
            last_point = [strut_x(1)/1000 strut_y(1)/1000];% last point for connection to strut 2
    end        
    elseif stent_er_config(j) == 2
        strut_x = strut_er(:,1) - feature_pos_er(j) + 2*pi*crimped_radius/(no_helix) + cr_er - WAB;
        strut_y = strut_er(:,2)+ (1-phi)*end_row_height;
        plot(RowAx, strut_x, strut_y,'LineWidth',1);hold on;
    if sw_er_strut_flag == 0
        sw_data_points = spline_data_points([strut_x strut_y]);
        py.sld_interface.create_spline(sw_data_points);
        first_point = [strut_x(1)/1000 strut_y(1)/1000];
        py.sld_interface.create_line(first_point,last_point);
        % straight line connection code here%
        last_point = [strut_x(end)/1000 strut_y(end)/1000];
    end
    elseif stent_er_config(j) == 3
        strut_x = crown_er_top(:,1) - feature_pos_er(j) + 2*pi*crimped_radius/(no_helix) + cr_er - WAB;
        strut_y = crown_er_top(:,2) + (1-phi)*end_row_height;
        plot(RowAx, strut_x, strut_y,'LineWidth',1);hold on;
    if sw_er_strut_flag == 0
        sw_data_points = spline_data_points([strut_x strut_y]);
        py.sld_interface.create_spline(sw_data_points);
        first_point = [strut_x(end)/1000 strut_y(end)/1000];
        py.sld_interface.create_line(first_point,last_point);
            % straight line connection code here%
        last_point = [strut_x(1)/1000 strut_y(1)/1000];
    end            
    elseif stent_er_config(j) == 4
        strut_x = strut_er(:,1) - feature_pos_er(j) + 2*pi*crimped_radius/(no_helix) + cr_er - WAB;
        strut_y = strut_er(:,2)+ (1-phi)*end_row_height;
        plot(RowAx, strut_x, strut_y,'LineWidth',1);hold on;
        if sw_er_strut_flag == 0
        extra_end_point = [strut_x(1) (strut_y(1)-0.5*extra_point_thresh)];
        strut_sp_initial = [strut_x strut_y];
        strut = [extra_end_point; strut_sp_initial];
        sw_data_points = spline_data_points(strut_sp_initial);% Adding an extra point at the end to maintain 0 slope at start
        py.sld_interface.create_spline(sw_data_points);
        first_point = [strut_x(end)/1000 strut_y(end)/1000];
        py.sld_interface.create_line(first_point,last_point);
        py.sld_interface.create_stent_segment_width(spline_id,stent_width/1000);%Create sketch width
        activesketch = py.sld_interface.get_active_sketch_name();
        sketch_name = activesketch;
        py.sld_interface.end_current_sketch();% Unit feature ends here so exit previous sketch        
        wrap_name = py.sld_interface.wrap_sketch_around_base_cylinder_v3(sketch_name, stent_thickness/1000);
        if wrap_name == "NOWRAP"
            error("Could not generate end row wrap")
        end
        cylinder_name_2 = py.sld_interface.create_base_cylinder(crimped_radius/1000, stent_length/1000);
        er_row_unit = py.sld_interface.combine_subtract(cylinder_name_2,wrap_name);
        er_row_unit = py.sld_interface.circpattern(er_row_unit,"Axis1",0.5*N_er_crowns);
        helix = sw_helix_name;
        helix = py.sld_interface.row_combine(helix, er_row_unit);% I created this function to combine rows in IR1 and reusing it. That is why the name of the function is row_combine. It just combines solid bodies
        sw_er_strut_flag = 1;
        end
    end        
    featurex = [featurex; strut_x];
    featurey = [featurey; strut_y];
    end

for j = 1:1:size(stent_er_config,2)% Creating the bottom End Row
    if stent_er_config(j) == 1
        strut_x = crown_er_bottom(:,1) - feature_pos_er(j) + 2*pi*crimped_radius/(no_helix) + cr_er - WAB;
        strut_y = crown_er_bottom(:,2) + (1-phi)*end_row_height;
        strut_x_rev = -1*strut_x - AAX*NS - 2*WAB - 2*crx - 2*pi*crimped_radius/(no_helix);
        strut_y_rev = -1*strut_y - AAY*NS - 2*HAB + 2*cry + HDAN;
        plot(RowAx, strut_x_rev, strut_y_rev,'LineWidth',1);hold on;
    if sw_er_rev_strut_flag == 0
            cylinder_name = py.sld_interface.create_base_cylinder(crimped_radius/1000, stent_length/1000);        
            strut_x_rev_initial = strut_x_rev - (-AAX*NS - 2*WAB - 2*crx - 2*pi*crimped_radius/(no_helix));
            strut_y_rev_initial = strut_y_rev - (- AAY*NS - 2*HAB + 2*cry + HDAN);
            strut_sp_initial = [strut_x_rev_initial strut_y_rev_initial];             
            sw_data_points = spline_data_points(strut_sp_initial);% Adding an extra point at the start to maintain 0 slope at start
            py.sld_interface.create_new_sketch();% Unit feature starts here so create new sketch at Front Plane
            spline_id = py.sld_interface.create_spline(sw_data_points);
            last_point = [strut_x_rev_initial(1)/1000 strut_y_rev_initial(1)/1000];% last point for connection to strut 2
    end
    elseif stent_er_config(j) == 2
        strut_x = strut_er(:,1) - feature_pos_er(j) + 2*pi*crimped_radius/(no_helix) + cr_er - WAB;
        strut_y = strut_er(:,2)+ (1-phi)*end_row_height;
        strut_x_rev = -1*strut_x - AAX*NS - 2*WAB - 2*crx - 2*pi*crimped_radius/(no_helix);
        strut_y_rev = -1*strut_y - AAY*NS - 2*HAB + 2*cry + HDAN;
        plot(RowAx, strut_x_rev, strut_y_rev,'LineWidth',1);hold on;
    if sw_er_rev_strut_flag == 0
        strut_x_rev_initial = strut_x_rev - (-AAX*NS - 2*WAB - 2*crx - 2*pi*crimped_radius/(no_helix));
        strut_y_rev_initial = strut_y_rev - (- AAY*NS - 2*HAB + 2*cry + HDAN);
        strut_sp_initial = [strut_x_rev_initial strut_y_rev_initial]; 
        sw_data_points = spline_data_points(strut_sp_initial);
        py.sld_interface.create_spline(sw_data_points);
        first_point = [strut_x_rev_initial(1)/1000 strut_y_rev_initial(1)/1000];
        py.sld_interface.create_line(first_point,last_point);
            % straight line connection code here%
        last_point = [strut_x_rev_initial(end)/1000 strut_y_rev_initial(end)/1000];
    end
    elseif stent_er_config(j) == 3
        strut_x = crown_er_top(:,1) - feature_pos_er(j) + 2*pi*crimped_radius/(no_helix) + cr_er - WAB;
        strut_y = crown_er_top(:,2) + (1-phi)*end_row_height;        
        strut_x_rev = -1*strut_x - AAX*NS - 2*WAB - 2*crx - 2*pi*crimped_radius/(no_helix);
        strut_y_rev = -1*strut_y - AAY*NS - 2*HAB + 2*cry + HDAN;        
        plot(RowAx, strut_x_rev, strut_y_rev,'LineWidth',1);hold on;
    if sw_er_rev_strut_flag == 0
        strut_x_rev_initial = strut_x_rev - (-AAX*NS - 2*WAB - 2*crx - 2*pi*crimped_radius/(no_helix));
        strut_y_rev_initial = strut_y_rev - (- AAY*NS - 2*HAB + 2*cry + HDAN);
        strut_sp_initial = [strut_x_rev_initial strut_y_rev_initial]; 
        sw_data_points = spline_data_points(strut_sp_initial);
        py.sld_interface.create_spline(sw_data_points);
        first_point = [strut_x_rev_initial(end)/1000 strut_y_rev_initial(end)/1000];
        py.sld_interface.create_line(first_point,last_point);
            % straight line connection code here%
        last_point = [strut_x_rev_initial(1)/1000 strut_y_rev_initial(1)/1000];
    end
    elseif stent_er_config(j) == 4
        strut_x = strut_er(:,1) - feature_pos_er(j) + 2*pi*crimped_radius/(no_helix) + cr_er - WAB;
        strut_y = strut_er(:,2)+ (1-phi)*end_row_height;
        strut_x_rev = -1*strut_x - AAX*NS - 2*WAB - 2*crx - 2*pi*crimped_radius/(no_helix);
        strut_y_rev = -1*strut_y - AAY*NS - 2*HAB + 2*cry + HDAN;
        plot(RowAx, strut_x_rev, strut_y_rev,'LineWidth',1);hold on; 
        if sw_er_rev_strut_flag == 0
        strut_x_rev_initial = strut_x_rev - (-AAX*NS - 2*WAB - 2*crx - 2*pi*crimped_radius/(no_helix));
        strut_y_rev_initial = strut_y_rev - (- AAY*NS - 2*HAB + 2*cry + HDAN);
        extra_end_point = [strut_x_rev_initial(1) (strut_y_rev_initial(1)+extra_point_thresh)];
        strut_sp_initial = [strut_x_rev_initial strut_y_rev_initial];
        strut = [extra_end_point; strut_sp_initial];
        sw_data_points = spline_data_points(strut_sp_initial);% Adding an extra point at the end to maintain 0 slope at start
        py.sld_interface.create_spline(sw_data_points);
        first_point = [strut_x_rev_initial(end)/1000 strut_y_rev_initial(end)/1000];
        py.sld_interface.create_line(first_point,last_point);
        py.sld_interface.create_stent_segment_width(spline_id,stent_width/1000);%Create sketch width
        activesketch = py.sld_interface.get_active_sketch_name();
        sketch_name = activesketch;
        py.sld_interface.end_current_sketch();% Unit feature ends here so exit previous sketch        
        wrap_name = py.sld_interface.wrap_sketch_around_base_cylinder_v3(sketch_name, stent_thickness/1000);
        if wrap_name == "NOWRAP"
            error("Could not generate end row bottom wrap")
        end
        cylinder_name_2 = py.sld_interface.create_base_cylinder(crimped_radius/1000, stent_length/1000);
        er_row_unit = py.sld_interface.combine_subtract(cylinder_name_2,wrap_name);
        er_row_unit = py.sld_interface.circpattern(er_row_unit,"Axis1",0.5*N_er_crowns);
        y_shift_3d_strut = - AAY*NS - 2*HAB + 2*cry + HDAN;
        theta_shift_3d_strut = (-AAX*NS - 2*WAB - 2*crx - 2*pi*crimped_radius/(no_helix))/crimped_radius;
        er_row_unit = py.sld_interface.move(er_row_unit, "Axis1", y_shift_3d_strut/1000,theta_shift_3d_strut);
        helix = py.sld_interface.row_combine(helix, er_row_unit);% I created this function to combine rows in IR1 and reusing it. That is why the name of the function is row_combine. It just combines solid bodies
        sw_er_rev_strut_flag = 1;
        end        
    end
    featurex_rev = [featurex_rev; strut_x_rev];
    featurey_rev = [featurey_rev; strut_y_rev];
end

    sw_helix_name = helix;
    er = [featurex featurey];
    er_rev = [featurex_rev featurey_rev];
end 

function [transition_row, transition_row_rev, sw_helix_name] = transition_row_builder(strut_ap, strut_pq, strut_qr, strut_rs, strut_st, strut_tu, strut_uc, HAP, WAP, HPQ, WPQ, HQR, WQR, HRS, WRS, HST, WST, HTU, WTU, crimped_radius, no_helix, AAX, AAY, HAB, WAB, crx, cry, HDAN, NS, sw_helix_name, WBC, HBC, extra_point_thresh, stent_width, stent_thickness, stent_length, RowAx)
transition_row = [];% will call transition row as tr from now on in this function
transition_row_rev = [];
helix_shift_x = 0;% The x shift in TR to accommodate more than one helix in the stent
featurex = [];
featurey = [];
tr_no = [];
featurex_rev = [];
featurey_rev = [];
tr_no_rev = [];
unit_tr_strut_list = string();
sw_tr_strut_flag = 0;
sw_tr_rev_strut_flag = 0;

for i = 1:1:no_helix% To generate top transition rows
% starting with AP Code 21
    strut_x = strut_ap(:,1) + helix_shift_x;
    strut_y = strut_ap(:,2);
    plot(RowAx, strut_x, strut_y,'LineWidth',1);hold on;
    tr_no = [tr_no;i*ones(size(strut_x))];
    featurex = [featurex; strut_x];
    featurey = [featurey; strut_y];
    if sw_tr_strut_flag == 0
        cylinder_name = py.sld_interface.create_base_cylinder(crimped_radius/1000, stent_length/1000);        
        extra_point_start = [-1*extra_point_thresh 0];
        strut_sp_initial = [strut_x strut_y];
        strut = [strut_sp_initial; extra_point_start];
        sw_data_points = spline_data_points(strut);% Adding an extra point at the start to maintain 0 slope at start
            % sw_data_points = spline_data_points([strut_x strut_y]);
        py.sld_interface.create_new_sketch();% Unit feature starts here so create new sketch at Front Plane
        spline_id = py.sld_interface.create_spline(sw_data_points);
            % straight line connection code here%
        last_point = [strut_x(1)/1000 strut_y(1)/1000];% last point for connection to strut 2
    end    
% Code 22 Strut PQ
    strut_x = strut_pq(:,1) + WAP +helix_shift_x;
    strut_y = strut_pq(:,2) - HAP;   
    plot(RowAx, strut_x, strut_y,'LineWidth',1);hold on;
    tr_no = [tr_no;i*ones(size(strut_x))];
    featurex = [featurex; strut_x];
    featurey = [featurey; strut_y];
    if sw_tr_strut_flag == 0
        sw_data_points = spline_data_points([strut_x strut_y]);
        py.sld_interface.create_spline(sw_data_points);
        first_point = [strut_x(1)/1000 strut_y(1)/1000];
        py.sld_interface.create_line(first_point,last_point);
            % straight line connection code here%
        last_point = [strut_x(end)/1000 strut_y(end)/1000];
    end
    % Code 23 Strut QR
    strut_x = strut_qr(:,1) + WAP + WPQ + helix_shift_x;
    strut_y = strut_qr(:,2) - HAP + HPQ;   
    plot(RowAx, strut_x, strut_y,'LineWidth',1);hold on;
    tr_no = [tr_no;i*ones(size(strut_x))];
    featurex = [featurex; strut_x];
    featurey = [featurey; strut_y];
    if sw_tr_strut_flag == 0
        sw_data_points = spline_data_points([strut_x strut_y]);
        py.sld_interface.create_spline(sw_data_points);
        first_point = [strut_x(end)/1000 strut_y(end)/1000];
        py.sld_interface.create_line(first_point,last_point);
            % straight line connection code here%
        last_point = [strut_x(1)/1000 strut_y(1)/1000];
    end
    % Code 24 Strut RS
    strut_x = strut_rs(:,1) + WAP + WPQ + WQR + helix_shift_x;
    strut_y = strut_rs(:,2) - HAP + HPQ - HQR;   
    plot(RowAx, strut_x, strut_y,'LineWidth',1);hold on;
    tr_no = [tr_no;i*ones(size(strut_x))];
    featurex = [featurex; strut_x];
    featurey = [featurey; strut_y];
    if sw_tr_strut_flag == 0
        sw_data_points = spline_data_points([strut_x strut_y]);
        py.sld_interface.create_spline(sw_data_points);
        first_point = [strut_x(1)/1000 strut_y(1)/1000];
        py.sld_interface.create_line(first_point,last_point);
            % straight line connection code here%
        last_point = [strut_x(end)/1000 strut_y(end)/1000];
    end
    % code 25 Strut ST
    strut_x = strut_st(:,1) + WAP + WPQ + WQR + WRS + helix_shift_x;
    strut_y = strut_st(:,2) - HAP + HPQ - HQR + HRS;   
    plot(RowAx, strut_x, strut_y,'LineWidth',1);hold on;
    tr_no = [tr_no;i*ones(size(strut_x))];
    featurex = [featurex; strut_x];
    featurey = [featurey; strut_y];
    if sw_tr_strut_flag == 0
        sw_data_points = spline_data_points([strut_x strut_y]);
        py.sld_interface.create_spline(sw_data_points);
        first_point = [strut_x(end)/1000 strut_y(end)/1000];
        py.sld_interface.create_line(first_point,last_point);
            % straight line connection code here%
        last_point = [strut_x(1)/1000 strut_y(1)/1000];
    end    
    % code 26 Strut TU
    strut_x = strut_tu(:,1) + WAP + WPQ + WQR + WRS + WST + helix_shift_x;
    strut_y = strut_tu(:,2) - HAP + HPQ - HQR + HRS - HST;   
    plot(RowAx, strut_x, strut_y,'LineWidth',1);hold on;
    tr_no = [tr_no;i*ones(size(strut_x))];
    featurex = [featurex; strut_x];
    featurey = [featurey; strut_y];
    if sw_tr_strut_flag == 0
        sw_data_points = spline_data_points([strut_x strut_y]);
        py.sld_interface.create_spline(sw_data_points);
        first_point = [strut_x(1)/1000 strut_y(1)/1000];
        py.sld_interface.create_line(first_point,last_point);
            % straight line connection code here%
        last_point = [strut_x(end)/1000 strut_y(end)/1000];
    end    
    % code 27 Strut UC
    strut_x = strut_uc(:,1) + WAP + WPQ + WQR + WRS + WST + WTU + helix_shift_x;
    strut_y = strut_uc(:,2) - HAP + HPQ - HQR + HRS - HST + HTU;   
    plot(RowAx, strut_x, strut_y,'LineWidth',1);hold on;
    tr_no = [tr_no;i*ones(size(strut_x))];
    featurex = [featurex; strut_x];
    featurey = [featurey; strut_y];
    if sw_tr_strut_flag == 0% To add the final three strut components but making the 3D version at the same place at the starting position so as to make sure that the strut is indeed made
        strut_x_initial = strut_x + helix_shift_x; % starting strut from 0,0 to avoid the strut crossing the 2pi boundary
        strut_y_initial = strut_y;
        strut = [strut_x_initial strut_y_initial];
        sw_data_points = spline_data_points(strut);% Adding an extra point at the start to maintain 0 slope at start
        py.sld_interface.create_spline(sw_data_points);
        first_point = [strut_x_initial(end)/1000 strut_y_initial(end)/1000];
        py.sld_interface.create_line(first_point,last_point); 
        py.sld_interface.create_stent_segment_width(spline_id,stent_width/1000);%Create sketch width
        activesketch = py.sld_interface.get_active_sketch_name();
        sketch_name = activesketch;
        py.sld_interface.end_current_sketch();% Unit feature ends here so exit previous sketch        
        wrap_name = py.sld_interface.wrap_sketch_around_base_cylinder_v3(sketch_name, stent_thickness/1000);
        if wrap_name == "NOWRAP"
            error("Could not generate transition row")
        end
        cylinder_name_2 = py.sld_interface.create_base_cylinder(crimped_radius/1000, stent_length/1000);
        tr_unit = py.sld_interface.combine_subtract(cylinder_name_2,wrap_name);
        unit_tr_strut_list(1) = tr_unit;
        sw_tr_strut_flag = 1;
    end
    if sw_tr_strut_flag == 1
        tr_unit_new = py.sld_interface.copy_move(tr_unit, "Axis1", 0,pi());
        unit_tr_strut_list(end+1) = tr_unit_new;
    end
transition_row = [transition_row; tr_no featurex featurey];
angle_shift = 2*pi()/no_helix;
helix_shift_x = helix_shift_x - angle_shift*crimped_radius;
sw_tr_strut_flag = 2;
end

helix_shift_x = 0;

for i = 1:1:no_helix% To generate bottom transition rows
% starting with AP Code 21
    strut_x = strut_ap(:,1) + helix_shift_x;
    strut_y = strut_ap(:,2);
% Generating the reverse strut
    strut_x_rev = -1*strut_x - AAX*NS - 2*WAB - 2*crx + 2*helix_shift_x;
    strut_y_rev = -1*strut_y - AAY*NS - 2*HAB + 2*cry + HDAN;
    plot(RowAx, strut_x_rev, strut_y_rev,'LineWidth',1);hold on;
    tr_no_rev = [tr_no_rev;i*ones(size(strut_x_rev))];
    featurex_rev = [featurex_rev; strut_x_rev];
    featurey_rev = [featurey_rev; strut_y_rev];
    if sw_tr_rev_strut_flag == 0
        cylinder_name = py.sld_interface.create_base_cylinder(crimped_radius/1000, stent_length/1000);        
        extra_point_start = [-1*extra_point_thresh 0];
        strut_x_rev_initial = strut_x_rev - (-AAX*NS - 2*WAB - 2*crx + 2*helix_shift_x);
        strut_y_rev_initial = strut_y_rev - (- AAY*NS - 2*HAB + 2*cry + HDAN);
        strut_sp_initial = [strut_x_rev_initial strut_y_rev_initial];
        strut = [strut_sp_initial; -1.*extra_point_start];
        sw_data_points = spline_data_points(strut);% Adding an extra point at the start to maintain 0 slope at start
            % sw_data_points = spline_data_points([strut_x strut_y]);
        py.sld_interface.create_new_sketch();% Unit feature starts here so create new sketch at Front Plane
        spline_id = py.sld_interface.create_spline(sw_data_points);
            % straight line connection code here%
        last_point = [strut_x_rev_initial(1)/1000 strut_y_rev_initial(1)/1000];% last point for connection to strut 2
    end    
% Code 22 Strut PQ
    strut_x = strut_pq(:,1) + WAP +helix_shift_x;
    strut_y = strut_pq(:,2) - HAP;
% Generating the reverse strut
    strut_x_rev = -1*strut_x - AAX*NS - 2*WAB - 2*crx + 2*helix_shift_x;
    strut_y_rev = -1*strut_y - AAY*NS - 2*HAB + 2*cry + HDAN;   
    plot(RowAx, strut_x_rev, strut_y_rev,'LineWidth',1);hold on;
    tr_no_rev = [tr_no_rev;i*ones(size(strut_x_rev))];
    featurex_rev = [featurex_rev; strut_x_rev];
    featurey_rev = [featurey_rev; strut_y_rev];    
    if sw_tr_rev_strut_flag == 0
        strut_x_rev_initial = strut_x_rev - (-AAX*NS - 2*WAB - 2*crx + 2*helix_shift_x);
        strut_y_rev_initial = strut_y_rev - (- AAY*NS - 2*HAB + 2*cry + HDAN);
        strut_sp_initial = [strut_x_rev_initial strut_y_rev_initial];
        sw_data_points = spline_data_points(strut_sp_initial);
        py.sld_interface.create_spline(sw_data_points);
        first_point = [strut_x_rev_initial(1)/1000 strut_y_rev_initial(1)/1000];
        py.sld_interface.create_line(first_point,last_point);
            % straight line connection code here%
        last_point = [strut_x_rev_initial(end)/1000 strut_y_rev_initial(end)/1000];
    end 
    % Code 23 Strut QR
    strut_x = strut_qr(:,1) + WAP + WPQ + helix_shift_x;
    strut_y = strut_qr(:,2) - HAP + HPQ;   
% Generating the reverse strut
    strut_x_rev = -1*strut_x - AAX*NS - 2*WAB - 2*crx + 2*helix_shift_x;
    strut_y_rev = -1*strut_y - AAY*NS - 2*HAB + 2*cry + HDAN;
    plot(RowAx, strut_x_rev, strut_y_rev,'LineWidth',1);hold on;
    tr_no_rev = [tr_no_rev;i*ones(size(strut_x_rev))];
    featurex_rev = [featurex_rev; strut_x_rev];
    featurey_rev = [featurey_rev; strut_y_rev];
    if sw_tr_rev_strut_flag == 0
        strut_x_rev_initial = strut_x_rev - (-AAX*NS - 2*WAB - 2*crx + 2*helix_shift_x);
        strut_y_rev_initial = strut_y_rev - (- AAY*NS - 2*HAB + 2*cry + HDAN);
        strut_sp_initial = [strut_x_rev_initial strut_y_rev_initial];
        sw_data_points = spline_data_points(strut_sp_initial);
        py.sld_interface.create_spline(sw_data_points);
        first_point = [strut_x_rev_initial(end)/1000 strut_y_rev_initial(end)/1000];
        py.sld_interface.create_line(first_point,last_point);
            % straight line connection code here%
        last_point = [strut_x_rev_initial(1)/1000 strut_y_rev_initial(1)/1000];
    end    
    % Code 24 Strut RS
    strut_x = strut_rs(:,1) + WAP + WPQ + WQR + helix_shift_x;
    strut_y = strut_rs(:,2) - HAP + HPQ - HQR;
% Generating the reverse strut
    strut_x_rev = -1*strut_x - AAX*NS - 2*WAB - 2*crx + 2*helix_shift_x;
    strut_y_rev = -1*strut_y - AAY*NS - 2*HAB + 2*cry + HDAN;
    plot(RowAx, strut_x_rev, strut_y_rev,'LineWidth',1);hold on;
    tr_no_rev = [tr_no_rev;i*ones(size(strut_x_rev))];
    featurex_rev = [featurex_rev; strut_x_rev];
    featurey_rev = [featurey_rev; strut_y_rev];
    if sw_tr_rev_strut_flag == 0
        strut_x_rev_initial = strut_x_rev - (-1*AAX*NS - 2*WAB - 2*crx + 2*helix_shift_x);
        strut_y_rev_initial = strut_y_rev - (- AAY*NS - 2*HAB + 2*cry + HDAN);
        strut_sp_initial = [strut_x_rev_initial strut_y_rev_initial];
        sw_data_points = spline_data_points(strut_sp_initial);
        py.sld_interface.create_spline(sw_data_points);
        first_point = [strut_x_rev_initial(1)/1000 strut_y_rev_initial(1)/1000];
        py.sld_interface.create_line(first_point,last_point);
            % straight line connection code here%
        last_point = [strut_x_rev_initial(end)/1000 strut_y_rev_initial(end)/1000];
    end    
    % code 25 Strut ST
    strut_x = strut_st(:,1) + WAP + WPQ + WQR + WRS + helix_shift_x;
    strut_y = strut_st(:,2) - HAP + HPQ - HQR + HRS;
% Generating the reverse strut
    strut_x_rev = -1*strut_x - AAX*NS - 2*WAB - 2*crx + 2*helix_shift_x;
    strut_y_rev = -1*strut_y - AAY*NS - 2*HAB + 2*cry + HDAN;
    plot(RowAx, strut_x_rev, strut_y_rev,'LineWidth',1);hold on;
    tr_no_rev = [tr_no_rev;i*ones(size(strut_x_rev))];
    featurex_rev = [featurex_rev; strut_x_rev];
    featurey_rev = [featurey_rev; strut_y_rev];
    if sw_tr_rev_strut_flag == 0
        strut_x_rev_initial = strut_x_rev - (-AAX*NS - 2*WAB - 2*crx + 2*helix_shift_x);
        strut_y_rev_initial = strut_y_rev - (- AAY*NS - 2*HAB + 2*cry + HDAN);
        strut_sp_initial = [strut_x_rev_initial strut_y_rev_initial];
        sw_data_points = spline_data_points(strut_sp_initial);
        py.sld_interface.create_spline(sw_data_points);
        first_point = [strut_x_rev_initial(end)/1000 strut_y_rev_initial(end)/1000];
        py.sld_interface.create_line(first_point,last_point);
            % straight line connection code here%
        last_point = [strut_x_rev_initial(1)/1000 strut_y_rev_initial(1)/1000];
    end    
    % code 26 Strut TU
    strut_x = strut_tu(:,1) + WAP + WPQ + WQR + WRS + WST + helix_shift_x;
    strut_y = strut_tu(:,2) - HAP + HPQ - HQR + HRS - HST;
% Generating the reverse strut
    strut_x_rev = -1*strut_x - AAX*NS - 2*WAB - 2*crx + 2*helix_shift_x;
    strut_y_rev = -1*strut_y - AAY*NS - 2*HAB + 2*cry + HDAN;
    plot(RowAx, strut_x_rev, strut_y_rev,'LineWidth',1);hold on;
    tr_no_rev = [tr_no_rev;i*ones(size(strut_x_rev))];
    featurex_rev = [featurex_rev; strut_x_rev];
    featurey_rev = [featurey_rev; strut_y_rev];    
    if sw_tr_rev_strut_flag == 0
        strut_x_rev_initial = strut_x_rev - (-AAX*NS - 2*WAB - 2*crx + 2*helix_shift_x);
        strut_y_rev_initial = strut_y_rev - (- AAY*NS - 2*HAB + 2*cry + HDAN);
        strut_sp_initial = [strut_x_rev_initial strut_y_rev_initial];
        sw_data_points = spline_data_points(strut_sp_initial);
        py.sld_interface.create_spline(sw_data_points);
        first_point = [strut_x_rev_initial(1)/1000 strut_y_rev_initial(1)/1000];
        py.sld_interface.create_line(first_point,last_point);
            % straight line connection code here%
        last_point = [strut_x_rev_initial(end)/1000 strut_y_rev_initial(end)/1000];
    end        
    % code 27 Strut UC
    strut_x = strut_uc(:,1) + WAP + WPQ + WQR + WRS + WST + WTU + helix_shift_x;
    strut_y = strut_uc(:,2) - HAP + HPQ - HQR + HRS - HST + HTU;
% Generating the reverse strut
    strut_x_rev = -1*strut_x - AAX*NS - 2*WAB - 2*crx + 2*helix_shift_x;
    strut_y_rev = -1*strut_y - AAY*NS - 2*HAB + 2*cry + HDAN;
    plot(RowAx, strut_x_rev, strut_y_rev,'LineWidth',1);hold on;
    tr_no_rev = [tr_no_rev;i*ones(size(strut_x_rev))];
    featurex_rev = [featurex_rev; strut_x_rev];
    featurey_rev = [featurey_rev; strut_y_rev];
    if sw_tr_rev_strut_flag == 0% To add the final three strut components but making the 3D version at the same place at the starting position so as to make sure that the strut is indeed made
        strut_x_rev_initial = strut_x_rev - (-AAX*NS - 2*WAB - 2*crx + 2*helix_shift_x);
        strut_y_rev_initial = strut_y_rev - (- AAY*NS - 2*HAB + 2*cry + HDAN);
        strut_sp_initial = [strut_x_rev_initial strut_y_rev_initial];
        sw_data_points = spline_data_points(strut_sp_initial);% Adding an extra point at the start to maintain 0 slope at start
        py.sld_interface.create_spline(sw_data_points);
        first_point = [strut_x_rev_initial(end)/1000 strut_y_rev_initial(end)/1000];
        py.sld_interface.create_line(first_point,last_point); 
        py.sld_interface.create_stent_segment_width(spline_id,stent_width/1000);%Create sketch width
        activesketch = py.sld_interface.get_active_sketch_name();
        sketch_name = activesketch;
        py.sld_interface.end_current_sketch();% Unit feature ends here so exit previous sketch        
        wrap_name = py.sld_interface.wrap_sketch_around_base_cylinder_v3(sketch_name, stent_thickness/1000);
        if wrap_name == "NOWRAP"
            error("Could not generate transition row")
        end
        cylinder_name_2 = py.sld_interface.create_base_cylinder(crimped_radius/1000, stent_length/1000);
        tr_unit = py.sld_interface.combine_subtract(cylinder_name_2,wrap_name);
        y_shift_3d_strut = - AAY*NS - 2*HAB + 2*cry + HDAN;
        theta_shift_3d_strut = (-AAX*NS - 2*WAB - 2*crx + 2*helix_shift_x)/crimped_radius;
        tr_unit = py.sld_interface.move(tr_unit, "Axis1", y_shift_3d_strut/1000,theta_shift_3d_strut);
        unit_tr_strut_list(end+1) = tr_unit;
        sw_tr_rev_strut_flag = 1;
    end
    if sw_tr_rev_strut_flag == 1
        tr_unit_new = py.sld_interface.copy_move(tr_unit, "Axis1", 0,pi());
        unit_tr_strut_list(end+1) = tr_unit_new;
    end
transition_row_rev = [transition_row_rev; tr_no_rev featurex_rev featurey_rev];
angle_shift = 2*pi()/no_helix;
helix_shift_x = helix_shift_x - angle_shift*crimped_radius;
sw_tr_rev_strut_flag = 2;
end

helix = sw_helix_name;

for k=1:1:size(unit_tr_strut_list,2)
    helix = py.sld_interface.row_combine(helix, unit_tr_strut_list(k));% I created this function to combine rows in IR1 and reusing it. That is why the name of the function is row_combine. It just combines solid bodies
end

sw_helix_name = helix;
end

function [helix, sw_helix_name] = helix_builder(strut_ab, strut_bc, strut_cd, strut_da, connector, cp_width_ab, cp_height_ab, cp_width_da, cp_height_da, cr_da, stent_helix_config, feature_pos_heliix, crimped_radius, no_helix, cr_yscale, NS, k, stent_length, stent_thickness, stent_width, extra_point_thresh, AAX, AAY, WDA, HDA, RowAx)
helix = [];
helix_shift_x = 0;% The x shift in helix to accommodate more than one helix in the stent
featurex = [];
featurey = [];
helix_no = [];
sw_unit_strut_flag = 0; % A flag that describes the status of a unit strut - 0 means that the 3D unit strut is not created fully yet. 1 means that 3D unit strut is created and its copies are being created to make the helix. Code 2 means the helix is already constructed and no need to create new unit struts. Code 9 to create the last three strut components to maintain symmetry from both sides
sw_connector_flag = 0;% -1 means that connector has not been generate and should not be generated, 0 means that connector should be generated through sketch, 1 menas that a baseline 3D connector exists and it should be copied. 2 means all connectors are created 
cylinder_name = py.sld_interface.create_base_cylinder(crimped_radius/1000, stent_length/1000);% Create the base crimped cylinder of crimped radius
py.sld_interface.create_base_cylinder_axis();% Create cylinder axis
rel_y_shift = -1*AAY;
rel_theta_shift = -1*(AAX/crimped_radius);

for i = 1:1:no_helix
stent_unit_name_list = string();
connector_unit_name_list = string();
% starting with AB Code 1
N_curr = 0;% Current unit strut being constructed. This is to allow for k = 5 to work .
len = size(stent_helix_config);
for j = 1:1:len(2)
if size(stent_unit_name_list,2) == NS
    sw_unit_strut_flag = 9;
end
if stent_helix_config(j) == 1
    N_curr = N_curr + 1;
    strut_x = strut_ab(:,1) + feature_pos_heliix(j,1) + helix_shift_x;
    strut_y = strut_ab(:,2) + feature_pos_heliix(j,2);
    plot(RowAx, strut_x, strut_y,'LineWidth',1);hold on;
    if sw_unit_strut_flag == 0
            extra_point_start = [extra_point_thresh 0];
            extra_point_end = [(strut_x(end)-extra_point_thresh) strut_y(end)];
            strut_sp_initial = [strut_x strut_y];
            strut = [extra_point_start;strut_sp_initial; extra_point_end];
            [strut1, strut2] = divide_into_two(strut);
            sw_data_points1 = spline_data_points(strut1);% Adding an extra point at the start to maintain 0 slope at start
            py.sld_interface.create_new_sketch();% Unit feature starts here so create new sketch at Front Plane
            spline_id = py.sld_interface.create_spline(sw_data_points1);
            py.sld_interface.create_stent_segment_width(spline_id,stent_width/1000);%Create sketch width
            activesketch = py.sld_interface.get_active_sketch_name();
            sketch_name = activesketch;
            py.sld_interface.end_current_sketch();% Unit feature ends here so exit previous sketch
            wrap_name = py.sld_interface.wrap_sketch_around_base_cylinder_v3(sketch_name, stent_thickness/1000);
            if wrap_name == "NOWRAP"
                error("Could not generate AB top part of base hellical wrap")
            end
            cylinder_name_2 = py.sld_interface.create_base_cylinder(crimped_radius/1000, stent_length/1000);
            stent_unit_AB_top = py.sld_interface.combine_subtract(cylinder_name_2,wrap_name);
            cylinder_name = py.sld_interface.create_base_cylinder(crimped_radius/1000, stent_length/1000);
            sw_data_points2 = spline_data_points(strut2);% Adding an extra point at the start to maintain 0 slope at start
            py.sld_interface.create_new_sketch();% Unit feature starts here so create new sketch at Front Plane
            spline_id = py.sld_interface.create_spline(sw_data_points2);
            py.sld_interface.create_stent_segment_width(spline_id,stent_width/1000);%Create sketch width
            activesketch = py.sld_interface.get_active_sketch_name();
            sketch_name = activesketch;
            py.sld_interface.end_current_sketch();% Unit feature ends here so exit previous sketch
            wrap_name = py.sld_interface.wrap_sketch_around_base_cylinder_v3(sketch_name, stent_thickness/1000);
            if wrap_name == "NOWRAP"
                error("Could not generate AB bottom part of base hellical wrap")
            end
            cylinder_name_2 = py.sld_interface.create_base_cylinder(crimped_radius/1000, stent_length/1000);
            stent_unit_AB_bottom = py.sld_interface.combine_subtract(cylinder_name_2,wrap_name);
            stent_unit_AB = py.sld_interface.row_combine(stent_unit_AB_top, stent_unit_AB_bottom);
    elseif sw_unit_strut_flag == 9% To add the final three strut components but making the 3D version at the same place at the starting position so as to make sure that the strut is indeed made
            cylinder_name = py.sld_interface.create_base_cylinder(crimped_radius/1000, stent_length/1000);        
            strut_x_initial = strut_x - feature_pos_heliix(j,1) - helix_shift_x; % starting strut from 0,0 to avoid the strut crossing the 2pi boundary
            strut_y_initial = strut_y - feature_pos_heliix(j,2); 
            extra_point_start = [extra_point_thresh 0];
            extra_point_end = [(strut_x_initial(end)-extra_point_thresh) strut_y_initial(end)];
            strut_sp_initial = [strut_x_initial strut_y_initial];
            strut = [extra_point_start;strut_sp_initial; extra_point_end];
            [strut1, strut2] = divide_into_two(strut);
            sw_data_points1 = spline_data_points(strut1);% Adding an extra point at the start to maintain 0 slope at start
            py.sld_interface.create_new_sketch();% Unit feature starts here so create new sketch at Front Plane
            spline_id = py.sld_interface.create_spline(sw_data_points1);
            py.sld_interface.create_stent_segment_width(spline_id,stent_width/1000);%Create sketch width
            activesketch = py.sld_interface.get_active_sketch_name();
            sketch_name = activesketch;
            py.sld_interface.end_current_sketch();% Unit feature ends here so exit previous sketch
            wrap_name = py.sld_interface.wrap_sketch_around_base_cylinder_v3(sketch_name, stent_thickness/1000);
            if wrap_name == "NOWRAP"
                error("Could not generate AB top part of bottom hellical wrap")
            end
            cylinder_name_2 = py.sld_interface.create_base_cylinder(crimped_radius/1000, stent_length/1000);
            stent_unit_AB_top = py.sld_interface.combine_subtract(cylinder_name_2,wrap_name);
            cylinder_name = py.sld_interface.create_base_cylinder(crimped_radius/1000, stent_length/1000);
            sw_data_points2 = spline_data_points(strut2);% Adding an extra point at the start to maintain 0 slope at start
            py.sld_interface.create_new_sketch();% Unit feature starts here so create new sketch at Front Plane
            spline_id = py.sld_interface.create_spline(sw_data_points2);
            py.sld_interface.create_stent_segment_width(spline_id,stent_width/1000);%Create sketch width
            activesketch = py.sld_interface.get_active_sketch_name();
            sketch_name = activesketch;
            py.sld_interface.end_current_sketch();% Unit feature ends here so exit previous sketch
            wrap_name = py.sld_interface.wrap_sketch_around_base_cylinder_v3(sketch_name, stent_thickness/1000);
            if wrap_name == "NOWRAP"
                error("Could not generate AB bottom part of end hellical wrap")
            end
            cylinder_name_2 = py.sld_interface.create_base_cylinder(crimped_radius/1000, stent_length/1000);
            stent_unit_AB_bottom = py.sld_interface.combine_subtract(cylinder_name_2,wrap_name);
            stent_unit_AB = py.sld_interface.row_combine(stent_unit_AB_top, stent_unit_AB_bottom);
    end
elseif stent_helix_config(j) == 2
% Code 2 Strut BC
    strut_x = strut_bc(:,1) - cp_width_ab + feature_pos_heliix(j,1) + helix_shift_x;
    strut_y = strut_bc(:,2) - cp_height_ab + feature_pos_heliix(j,2);
    plot(RowAx, strut_x, strut_y,'LineWidth',1);hold on;
    if sw_unit_strut_flag == 0
        extra_point_start = [(strut_x(1)+ extra_point_thresh) strut_y(1)];
        extra_point_end = [(strut_x(end)- extra_point_thresh) strut_y(end)];
        strut_sp_initial = [strut_x strut_y];
        strut = [extra_point_start;strut_sp_initial;extra_point_end];
        [strut1, strut2, strut3] = divide_into_three(strut);        
        % Creating BC bottom wrap 
        cylinder_name = py.sld_interface.create_base_cylinder(crimped_radius/1000, stent_length/1000);
        sw_data_points1 = spline_data_points(strut1);% Adding an extra point at the start to maintain 0 slope at start
        py.sld_interface.create_new_sketch();% Unit feature starts here so create new sketch at Front Plane
        spline_id = py.sld_interface.create_spline(sw_data_points1);
        py.sld_interface.create_stent_segment_width(spline_id,stent_width/1000);%Create sketch width
        activesketch = py.sld_interface.get_active_sketch_name();
        sketch_name = activesketch;
        py.sld_interface.end_current_sketch();% Unit feature ends here so exit previous sketch
        wrap_name = py.sld_interface.wrap_sketch_around_base_cylinder_v3(sketch_name, stent_thickness/1000);
        if wrap_name == "NOWRAP"
            error("Could not generate BC bottom part of base hellical wrap")
        end
        cylinder_name_2 = py.sld_interface.create_base_cylinder(crimped_radius/1000, stent_length/1000);
        stent_unit_BC_bottom = py.sld_interface.combine_subtract(cylinder_name_2,wrap_name);        
        % Creating BC middle wrap 
        cylinder_name = py.sld_interface.create_base_cylinder(crimped_radius/1000, stent_length/1000);
        sw_data_points2 = spline_data_points(strut2);% Adding an extra point at the start to maintain 0 slope at start
        py.sld_interface.create_new_sketch();% Unit feature starts here so create new sketch at Front Plane
        spline_id = py.sld_interface.create_spline(sw_data_points2);
        py.sld_interface.create_stent_segment_width(spline_id,stent_width/1000);%Create sketch width
        activesketch = py.sld_interface.get_active_sketch_name();
        sketch_name = activesketch;
        py.sld_interface.end_current_sketch();% Unit feature ends here so exit previous sketch
        wrap_name = py.sld_interface.wrap_sketch_around_base_cylinder_v3(sketch_name, stent_thickness/1000);
        if wrap_name == "NOWRAP"
            error("Could not generate BC middle part of base hellical wrap")
        end
        cylinder_name_2 = py.sld_interface.create_base_cylinder(crimped_radius/1000, stent_length/1000);
        stent_unit_BC_mid = py.sld_interface.combine_subtract(cylinder_name_2,wrap_name);
        stent_unit_BC_mid_bottom = py.sld_interface.row_combine(stent_unit_BC_mid, stent_unit_BC_bottom);
        % Creating BC top wrap 
        cylinder_name = py.sld_interface.create_base_cylinder(crimped_radius/1000, stent_length/1000);
        sw_data_points3 = spline_data_points(strut3);% Adding an extra point at the start to maintain 0 slope at start
        py.sld_interface.create_new_sketch();% Unit feature starts here so create new sketch at Front Plane
        spline_id = py.sld_interface.create_spline(sw_data_points3);
        py.sld_interface.create_stent_segment_width(spline_id,stent_width/1000);%Create sketch width
        activesketch = py.sld_interface.get_active_sketch_name();
        sketch_name = activesketch;
        py.sld_interface.end_current_sketch();% Unit feature ends here so exit previous sketch
        wrap_name = py.sld_interface.wrap_sketch_around_base_cylinder_v3(sketch_name, stent_thickness/1000);
        if wrap_name == "NOWRAP"
            error("Could not generate BC top part of base hellical wrap")
        end
        cylinder_name_2 = py.sld_interface.create_base_cylinder(crimped_radius/1000, stent_length/1000);
        stent_unit_BC_top = py.sld_interface.combine_subtract(cylinder_name_2,wrap_name);
        stent_unit_BC = py.sld_interface.row_combine(stent_unit_BC_mid_bottom, stent_unit_BC_top);
        stent_unit_AC = py.sld_interface.row_combine(stent_unit_AB, stent_unit_BC);        
    elseif sw_unit_strut_flag == 9% To add the final three strut components but making the 3D version at the same place at the starting position so as to make sure that the strut is indeed made
        strut_x_initial = strut_x - feature_pos_heliix(j,1) - helix_shift_x; % starting strut from 0,0 to avoid the strut crossing the 2pi boundary
        strut_y_initial = strut_y - feature_pos_heliix(j,2);
        extra_point_start = [(strut_x_initial(1)+ extra_point_thresh) strut_y_initial(1)];
        extra_point_end = [(strut_x_initial(end)- extra_point_thresh) strut_y_initial(end)];
        strut_sp_initial = [strut_x_initial strut_y_initial];
        strut = [extra_point_start;strut_sp_initial;extra_point_end];
        [strut1, strut2, strut3] = divide_into_three(strut);        
        % Creating BC bottom wrap 
        cylinder_name = py.sld_interface.create_base_cylinder(crimped_radius/1000, stent_length/1000);
        sw_data_points1 = spline_data_points(strut1);% Adding an extra point at the start to maintain 0 slope at start
        py.sld_interface.create_new_sketch();% Unit feature starts here so create new sketch at Front Plane
        spline_id = py.sld_interface.create_spline(sw_data_points1);
        py.sld_interface.create_stent_segment_width(spline_id,stent_width/1000);%Create sketch width
        activesketch = py.sld_interface.get_active_sketch_name();
        sketch_name = activesketch;
        py.sld_interface.end_current_sketch();% Unit feature ends here so exit previous sketch
        wrap_name = py.sld_interface.wrap_sketch_around_base_cylinder_v3(sketch_name, stent_thickness/1000);
        if wrap_name == "NOWRAP"
            error("Could not generate BC bottom part of end hellical wrap")
        end
        cylinder_name_2 = py.sld_interface.create_base_cylinder(crimped_radius/1000, stent_length/1000);
        stent_unit_BC_bottom = py.sld_interface.combine_subtract(cylinder_name_2,wrap_name);        
        % Creating BC middle wrap 
        cylinder_name = py.sld_interface.create_base_cylinder(crimped_radius/1000, stent_length/1000);
        sw_data_points2 = spline_data_points(strut2);% Adding an extra point at the start to maintain 0 slope at start
        py.sld_interface.create_new_sketch();% Unit feature starts here so create new sketch at Front Plane
        spline_id = py.sld_interface.create_spline(sw_data_points2);
        py.sld_interface.create_stent_segment_width(spline_id,stent_width/1000);%Create sketch width
        activesketch = py.sld_interface.get_active_sketch_name();
        sketch_name = activesketch;
        py.sld_interface.end_current_sketch();% Unit feature ends here so exit previous sketch
        wrap_name = py.sld_interface.wrap_sketch_around_base_cylinder_v3(sketch_name, stent_thickness/1000);
        if wrap_name == "NOWRAP"
            error("Could not generate BC middle part of end hellical wrap")
        end
        cylinder_name_2 = py.sld_interface.create_base_cylinder(crimped_radius/1000, stent_length/1000);
        stent_unit_BC_mid = py.sld_interface.combine_subtract(cylinder_name_2,wrap_name);
        stent_unit_BC_mid_bottom = py.sld_interface.row_combine(stent_unit_BC_mid, stent_unit_BC_bottom);
        % Creating BC top wrap 
        cylinder_name = py.sld_interface.create_base_cylinder(crimped_radius/1000, stent_length/1000);
        sw_data_points3 = spline_data_points(strut3);% Adding an extra point at the start to maintain 0 slope at start
        py.sld_interface.create_new_sketch();% Unit feature starts here so create new sketch at Front Plane
        spline_id = py.sld_interface.create_spline(sw_data_points3);
        py.sld_interface.create_stent_segment_width(spline_id,stent_width/1000);%Create sketch width
        activesketch = py.sld_interface.get_active_sketch_name();
        sketch_name = activesketch;
        py.sld_interface.end_current_sketch();% Unit feature ends here so exit previous sketch
        wrap_name = py.sld_interface.wrap_sketch_around_base_cylinder_v3(sketch_name, stent_thickness/1000);
        if wrap_name == "NOWRAP"
            error("Could not generate BC top part of end hellical wrap")
        end
        cylinder_name_2 = py.sld_interface.create_base_cylinder(crimped_radius/1000, stent_length/1000);
        stent_unit_BC_top = py.sld_interface.combine_subtract(cylinder_name_2,wrap_name);
        stent_unit_BC = py.sld_interface.row_combine(stent_unit_BC_mid_bottom, stent_unit_BC_top);
        stent_unit_AC = py.sld_interface.row_combine(stent_unit_AB, stent_unit_BC);                
    end
elseif stent_helix_config(j) == 3
    % Code 3 Strut CD
    strut_x = strut_cd(:,1) - cp_width_ab - 2*cr_da + feature_pos_heliix(j,1) + helix_shift_x;
    strut_y = strut_cd(:,2) - cp_height_ab + 2*cr_da*cr_yscale + cp_height_da + feature_pos_heliix(j,2);
    plot(RowAx, strut_x, strut_y,'LineWidth',1);hold on;
    if sw_unit_strut_flag == 0
        cylinder_name = py.sld_interface.create_base_cylinder(crimped_radius/1000, stent_length/1000);        
        extra_point_start = [(strut_x(1)+extra_point_thresh) strut_y(1)];
        extra_point_end = [(strut_x(end)-extra_point_thresh) strut_y(end)];
        strut_sp_initial = [strut_x strut_y];
        strut = [extra_point_start;strut_sp_initial; extra_point_end];
        [strut1, strut2] = divide_into_two(strut);
        sw_data_points1 = spline_data_points(strut1);% 
        py.sld_interface.create_new_sketch();% Unit feature starts here so create new sketch at Front Plane
        spline_id = py.sld_interface.create_spline(sw_data_points1);
        py.sld_interface.create_stent_segment_width(spline_id,stent_width/1000);%Create sketch width
        activesketch = py.sld_interface.get_active_sketch_name();
        sketch_name = activesketch;
        py.sld_interface.end_current_sketch();% Unit feature ends here so exit previous sketch
        wrap_name = py.sld_interface.wrap_sketch_around_base_cylinder_v3(sketch_name, stent_thickness/1000);
        if wrap_name == "NOWRAP"
            error("Could not generate CD top part of base hellical wrap")
        end
        cylinder_name_2 = py.sld_interface.create_base_cylinder(crimped_radius/1000, stent_length/1000);
        stent_unit_CD_top = py.sld_interface.combine_subtract(cylinder_name_2,wrap_name);
        cylinder_name = py.sld_interface.create_base_cylinder(crimped_radius/1000, stent_length/1000);
        sw_data_points2 = spline_data_points(strut2);% Adding an extra point at the start to maintain 0 slope at start
        py.sld_interface.create_new_sketch();% Unit feature starts here so create new sketch at Front Plane
        spline_id = py.sld_interface.create_spline(sw_data_points2);
        py.sld_interface.create_stent_segment_width(spline_id,stent_width/1000);%Create sketch width
        activesketch = py.sld_interface.get_active_sketch_name();
        sketch_name = activesketch;
        py.sld_interface.end_current_sketch();% Unit feature ends here so exit previous sketch
        wrap_name = py.sld_interface.wrap_sketch_around_base_cylinder_v3(sketch_name, stent_thickness/1000);
        if wrap_name == "NOWRAP"
            error("Could not generate CD bottom part of base hellical wrap")
        end
        cylinder_name_2 = py.sld_interface.create_base_cylinder(crimped_radius/1000, stent_length/1000);
        stent_unit_CD_bottom = py.sld_interface.combine_subtract(cylinder_name_2,wrap_name);
        stent_unit_CD = py.sld_interface.row_combine(stent_unit_CD_top, stent_unit_CD_bottom);        
    elseif sw_unit_strut_flag == 9% To add the final three strut components but making the 3D version at the same place at the starting position so as to make sure that the strut is indeed made
        strut_x_initial = strut_x - feature_pos_heliix(j,1) - helix_shift_x; % starting strut from 0,0 to avoid the strut crossing the 2pi boundary
        strut_y_initial = strut_y - feature_pos_heliix(j,2);
        cylinder_name = py.sld_interface.create_base_cylinder(crimped_radius/1000, stent_length/1000);        
        extra_point_start = [(strut_x_initial(1)+extra_point_thresh) strut_y_initial(1)];
        extra_point_end = [(strut_x_initial(end)-extra_point_thresh) strut_y_initial(end)];
        strut_sp_initial = [strut_x_initial strut_y_initial];
        strut = [extra_point_start;strut_sp_initial; extra_point_end];
        [strut1, strut2] = divide_into_two(strut);
        sw_data_points1 = spline_data_points(strut1);% 
        py.sld_interface.create_new_sketch();% Unit feature starts here so create new sketch at Front Plane
        spline_id = py.sld_interface.create_spline(sw_data_points1);
        py.sld_interface.create_stent_segment_width(spline_id,stent_width/1000);%Create sketch width
        activesketch = py.sld_interface.get_active_sketch_name();
        sketch_name = activesketch;
        py.sld_interface.end_current_sketch();% Unit feature ends here so exit previous sketch
        wrap_name = py.sld_interface.wrap_sketch_around_base_cylinder_v3(sketch_name, stent_thickness/1000);
        if wrap_name == "NOWRAP"
            error("Could not generate CD top part of end hellical wrap")
        end
        cylinder_name_2 = py.sld_interface.create_base_cylinder(crimped_radius/1000, stent_length/1000);
        stent_unit_CD_top = py.sld_interface.combine_subtract(cylinder_name_2,wrap_name);
        cylinder_name = py.sld_interface.create_base_cylinder(crimped_radius/1000, stent_length/1000);
        sw_data_points2 = spline_data_points(strut2);% Adding an extra point at the start to maintain 0 slope at start
        py.sld_interface.create_new_sketch();% Unit feature starts here so create new sketch at Front Plane
        spline_id = py.sld_interface.create_spline(sw_data_points2);
        py.sld_interface.create_stent_segment_width(spline_id,stent_width/1000);%Create sketch width
        activesketch = py.sld_interface.get_active_sketch_name();
        sketch_name = activesketch;
        py.sld_interface.end_current_sketch();% Unit feature ends here so exit previous sketch
        wrap_name = py.sld_interface.wrap_sketch_around_base_cylinder_v3(sketch_name, stent_thickness/1000);
        if wrap_name == "NOWRAP"
            error("Could not generate CD bottom part of end hellical wrap")
        end
        cylinder_name_2 = py.sld_interface.create_base_cylinder(crimped_radius/1000, stent_length/1000);
        stent_unit_CD_bottom = py.sld_interface.combine_subtract(cylinder_name_2,wrap_name);
        stent_unit_CD = py.sld_interface.row_combine(stent_unit_CD_top, stent_unit_CD_bottom);         
        stent_unit = py.sld_interface.row_combine(stent_unit_AC, stent_unit_CD);
        end_comp_y_shift = feature_pos_heliix(j,2);
        end_comp_theta_shift = feature_pos_heliix(j,1)/(crimped_radius);
        stent_unit_new = py.sld_interface.move(stent_unit, "Axis1", end_comp_y_shift/1000,end_comp_theta_shift);
        stent_unit = stent_unit_new;
        stent_unit_name_list(end+1) = stent_unit;
        sw_unit_strut_flag = 1;
    end
elseif stent_helix_config(j) == 4
% Code 4 Strut DA
    strut_x = strut_da(:,1) - cp_width_ab - 2*cr_da - cp_width_ab + feature_pos_heliix(j,1) + helix_shift_x;
    strut_y = strut_da(:,2) - cp_height_ab + 2*cr_da*cr_yscale + cp_height_da - cp_height_ab + feature_pos_heliix(j,2);
    plot(RowAx, strut_x, strut_y,'LineWidth',1);hold on;
    if sw_unit_strut_flag == 1
            stent_unit_new = py.sld_interface.copy_move(stent_unit, "Axis1", rel_y_shift/1000,rel_theta_shift);
            stent_unit = stent_unit_new;
            stent_unit_name_list(end+1) = stent_unit;    
    elseif sw_unit_strut_flag == 0
        extra_point_start = [(strut_x(1)+ extra_point_thresh) strut_y(1)];
        extra_point_end = [(strut_x(end)- extra_point_thresh) strut_y(end)];
        strut_sp_initial = [strut_x strut_y];
        strut = [extra_point_start;strut_sp_initial;extra_point_end];
        [strut1, strut2, strut3] = divide_into_three(strut);        
        % Creating DA bottom wrap 
        cylinder_name = py.sld_interface.create_base_cylinder(crimped_radius/1000, stent_length/1000);
        sw_data_points1 = spline_data_points(strut1);% Adding an extra point at the start to maintain 0 slope at start
        py.sld_interface.create_new_sketch();% Unit feature starts here so create new sketch at Front Plane
        spline_id = py.sld_interface.create_spline(sw_data_points1);
        py.sld_interface.create_stent_segment_width(spline_id,stent_width/1000);%Create sketch width
        activesketch = py.sld_interface.get_active_sketch_name();
        sketch_name = activesketch;
        py.sld_interface.end_current_sketch();% Unit feature ends here so exit previous sketch
        wrap_name = py.sld_interface.wrap_sketch_around_base_cylinder_v3(sketch_name, stent_thickness/1000);
        if wrap_name == "NOWRAP"
            error("Could not generate DA bottom part of base hellical wrap")
        end
        cylinder_name_2 = py.sld_interface.create_base_cylinder(crimped_radius/1000, stent_length/1000);
        stent_unit_DA_bottom = py.sld_interface.combine_subtract(cylinder_name_2,wrap_name);        
        % Creating DA middle wrap 
        cylinder_name = py.sld_interface.create_base_cylinder(crimped_radius/1000, stent_length/1000);
        sw_data_points2 = spline_data_points(strut2);% Adding an extra point at the start to maintain 0 slope at start
        py.sld_interface.create_new_sketch();% Unit feature starts here so create new sketch at Front Plane
        spline_id = py.sld_interface.create_spline(sw_data_points2);
        py.sld_interface.create_stent_segment_width(spline_id,stent_width/1000);%Create sketch width
        activesketch = py.sld_interface.get_active_sketch_name();
        sketch_name = activesketch;
        py.sld_interface.end_current_sketch();% Unit feature ends here so exit previous sketch
        wrap_name = py.sld_interface.wrap_sketch_around_base_cylinder_v3(sketch_name, stent_thickness/1000);
        if wrap_name == "NOWRAP"
            error("Could not generate DA middle part of base hellical wrap")
        end
        cylinder_name_2 = py.sld_interface.create_base_cylinder(crimped_radius/1000, stent_length/1000);
        stent_unit_DA_mid = py.sld_interface.combine_subtract(cylinder_name_2,wrap_name);
        stent_unit_DA_mid_bottom = py.sld_interface.row_combine(stent_unit_DA_mid, stent_unit_DA_bottom);
        % Creating DA top wrap 
        cylinder_name = py.sld_interface.create_base_cylinder(crimped_radius/1000, stent_length/1000);
        sw_data_points3 = spline_data_points(strut3);% Adding an extra point at the start to maintain 0 slope at start
        py.sld_interface.create_new_sketch();% Unit feature starts here so create new sketch at Front Plane
        spline_id = py.sld_interface.create_spline(sw_data_points3);
        py.sld_interface.create_stent_segment_width(spline_id,stent_width/1000);%Create sketch width
        activesketch = py.sld_interface.get_active_sketch_name();
        sketch_name = activesketch;
        py.sld_interface.end_current_sketch();% Unit feature ends here so exit previous sketch
        wrap_name = py.sld_interface.wrap_sketch_around_base_cylinder_v3(sketch_name, stent_thickness/1000);
        if wrap_name == "NOWRAP"
            error("Could not generate DA top part of base hellical wrap")
        end
        cylinder_name_2 = py.sld_interface.create_base_cylinder(crimped_radius/1000, stent_length/1000);
        stent_unit_DA_top = py.sld_interface.combine_subtract(cylinder_name_2,wrap_name);
        stent_unit_DA = py.sld_interface.row_combine(stent_unit_DA_mid_bottom, stent_unit_DA_top);
        stent_unit_CA = py.sld_interface.row_combine(stent_unit_CD, stent_unit_DA);
        stent_unit = py.sld_interface.row_combine(stent_unit_AC, stent_unit_CA);
        stent_unit_name_list(1) = stent_unit;        
        sw_unit_strut_flag = 1;
    end
elseif stent_helix_config(j) == 5
% Code 5 Connector
    if N_curr > NS - (floor(k/2) - 1)
        continue
    else
    conn_x = connector(:,1) - (cp_width_ab/2) + feature_pos_heliix(j,1) + helix_shift_x;
    conn_y = connector(:,2) - (cp_height_ab/2) + feature_pos_heliix(j,2);
    plot(RowAx, conn_x, conn_y,'LineWidth',1);hold on;
    strut_x = conn_x;
    strut_y = conn_y;
    if sw_connector_flag == 1
        connector_unit_new = py.sld_interface.copy_move(connector_unit, "Axis1", rel_y_shift/1000,rel_theta_shift);
        connector_unit = connector_unit_new;
        connector_unit_name_list(end+1) = connector_unit;        
    elseif sw_connector_flag == 0
        cylinder_name = py.sld_interface.create_base_cylinder(crimped_radius/1000, stent_length/1000);
        conn = [conn_x conn_y];
        sw_data_points = spline_data_points(conn);
        py.sld_interface.create_new_sketch();% create a new sketch at front plane
        spline_id = py.sld_interface.create_spline(sw_data_points);
        py.sld_interface.create_stent_segment_width(spline_id,stent_width/1000);%Create sketch width
        activesketch = py.sld_interface.get_active_sketch_name();
        sketch_name = activesketch;
        py.sld_interface.end_current_sketch();% Unit feature ends here so exit previous sketch        
        wrap_name = py.sld_interface.wrap_sketch_around_base_cylinder_v3(sketch_name, stent_thickness/1000);
        if wrap_name == "NOWRAP"
            error("Could not generate connector wrap")
        end
        cylinder_name_2 = py.sld_interface.create_base_cylinder(crimped_radius/1000, stent_length/1000);
        connector_unit = py.sld_interface.combine_subtract(cylinder_name_2,wrap_name);
        connector_unit_name_list(1) = connector_unit;
        sw_connector_flag = 1;
    end
    end
end
helix_no = [helix_no;i*ones(size(strut_x))];
featurex = [featurex; strut_x];
featurey = [featurey; strut_y];
end
helix = [helix_no featurex featurey];
angle_shift = 2*pi()/no_helix;
helix_shift_x = helix_shift_x - angle_shift*crimped_radius;
if sw_unit_strut_flag == 1
    prev_unit = stent_unit_name_list(1);
    for k=2:1:size(stent_unit_name_list,2)
        prev_unit = py.sld_interface.row_combine(prev_unit, stent_unit_name_list(k));% I created this function to combine rows in IR1 and reusing it. That is why the name of the function is row_combine. It just combines solid bodies
        helix1 = prev_unit;
    end
end
if sw_connector_flag == 1
    for k=1:1:size(connector_unit_name_list,2)
        helix1 = py.sld_interface.row_combine(helix1, connector_unit_name_list(k));% I created this function to combine rows in IR1 and reusing it. That is why the name of the function is row_combine. It just combines solid bodies
    end
end
if i==1
    helix_sw_y_shift = 0;
    helix_sw_theta_shift = pi();
    helix2 = py.sld_interface.copy_move(helix1, "Axis1", helix_sw_y_shift/1000,helix_sw_theta_shift);
    sw_helix_name = py.sld_interface.row_combine(helix1, helix2);
    sw_unit_strut_flag = 2;
    sw_connector_flag = 2;
end 
end
end

function [feature_pos_er] = feature_position_er(config,cr_er)
len = size(config);
current_position = cr_er; 

for j=1:len(2)
    feature_pos_er(j) = current_position;
    current_position = current_position + cr_er;
end
end

function [stent_er_config] = er_config_builder(N_crowns) % End row config builder
er = [];
base_unit = [1 2 3 4];

for i=1:1:0.5*N_crowns
    er = [er base_unit]; 
end

stent_er_config = er;
end

function [feature_pos] = feature_position(config,AAX, AAY)
len = size(config);
shift_x = -1*AAX;
shift_y = -1*AAY;
current_pos_x = AAX;
current_pos_y = AAY;

for j=1:len(2)
    if config(j) == 1
        current_pos_x = current_pos_x + shift_x;
        current_pos_y = current_pos_y + shift_y;
        feature_pos(j,1) = current_pos_x;
        feature_pos(j,2) = current_pos_y;
    else
        feature_pos(j,1) = current_pos_x;
        feature_pos(j,2) = current_pos_y;
    end

end
end

function [stent_helix_config] = helix_config_builder(N_struts)
%Nomeclature
% Strut AB
% Strut BC
% Strut CD
% Strut DA
% Connector
% base_unit = [1 2 3 4 5]; 
helix = [];
base_unit = [1 2 3 4 5];

for i=1:N_struts
    % if N_struts - i < no_connector_units_at_bottom_end
    %     base_unit = [1 2 3 4];
    % else
    %     base_unit = [1 2 3 4 5];
    % end
    helix = [helix base_unit]; 
end

helix = [helix 1 2 3];% Adding the extra AB strut at the end
stent_helix_config = helix;
end

function [transition_unit_strut] = generate_transition_unit_strut(strut_ap, strut_pq, strut_qr, strut_rs, strut_st, strut_tu, strut_uc, HAP,WAP, HPQ, WPQ, HQR, WQR, HRS, WRS, HST, WST, HTU, WTU)
featurex = [];
featurey = [];

% starting with AB Code21
strut_ap_x = strut_ap(:,1);
strut_ap_y = strut_ap(:,2);
featurex = [featurex; strut_ap_x];
featurey = [featurey; strut_ap_y];

% Code 22 Strut PQ
strut_pq_x = strut_pq(:,1) + WAP;
strut_pq_y = strut_pq(:,2) - HAP;
featurex = [featurex; strut_pq_x];
featurey = [featurey; strut_pq_y];

% Code 23 Strut QR
strut_qr_x = strut_qr(:,1) + WAP + WPQ;
strut_qr_y = strut_qr(:,2) - HAP + HPQ;
featurex = [featurex; strut_qr_x];
featurey = [featurey; strut_qr_y];

% Code 24 Strut RS
strut_rs_x = strut_rs(:,1) + WAP + WPQ + WQR;
strut_rs_y = strut_rs(:,2) - HAP + HPQ - HQR;
featurex = [featurex; strut_rs_x];
featurey = [featurey; strut_rs_y];

% Code 25 Strut ST
strut_st_x = strut_st(:,1) + WAP + WPQ + WQR + WRS;
strut_st_y = strut_st(:,2) - HAP + HPQ - HQR + HRS;
featurex = [featurex; strut_st_x];
featurey = [featurey; strut_st_y];

% Code 26 Strut TU
strut_tu_x = strut_tu(:,1) + WAP + WPQ + WQR + WRS + WST;
strut_tu_y = strut_tu(:,2) - HAP + HPQ - HQR + HRS - HST;
featurex = [featurex; strut_tu_x];
featurey = [featurey; strut_tu_y];

% Code 27 Strut UC
strut_uc_x = strut_uc(:,1) + WAP + WPQ + WQR + WRS + WST + WTU;
strut_uc_y = strut_uc(:,2) - HAP + HPQ - HQR + HRS - HST + HTU;
featurex = [featurex; strut_uc_x];
featurey = [featurey; strut_uc_y];
transition_unit_strut = [featurex featurey];
end

function [unit_strut] = generate_unit_strut(strut_ab, strut_bc, strut_cd, strut_da, connector,  cp_width_ab, cp_height_ab, cp_width_da, cp_height_da, cr_da, cr_yscale)
featurex = [];
featurey = [];
feature_connx = [];
feature_conny = [];

% starting with AB Code 1
strut_ab_x = strut_ab(:,1);
strut_ab_y = strut_ab(:,2);
featurex = [featurex; strut_ab_x];
featurey = [featurey; strut_ab_y];

% Code 2 Strut BC
strut_bc_x = strut_bc(:,1) - cp_width_ab;
strut_bc_y = strut_bc(:,2) - cp_height_ab;
featurex = [featurex; strut_bc_x];
featurey = [featurey; strut_bc_y];

% Code 3 Strut CD
strut_cd_x = strut_cd(:,1) - cp_width_ab - 2*cr_da;
strut_cd_y = strut_cd(:,2) - cp_height_ab + 2*cr_da*cr_yscale + cp_height_da;
featurex = [featurex; strut_cd_x];
featurey = [featurey; strut_cd_y];

% Code 4 Strut DA
strut_da_x = strut_da(:,1) - cp_width_ab - 2*cr_da - cp_width_ab;
strut_da_y = strut_da(:,2) - cp_height_ab + 2*cr_da*cr_yscale + cp_height_da - cp_height_ab;
featurex = [featurex; strut_da_x];
featurey = [featurey; strut_da_y];

% Code 5 connector
conn_x = connector(:,1) - (cp_width_ab/2);
conn_y = connector(:,2) - (cp_height_ab/2);
feature_connx = [feature_connx; conn_x];
feature_conny = [feature_conny; conn_y];

unit_strut = [featurex featurey];
unit_conn = [feature_connx feature_conny];
end

function [end_connector] = end_connector_generator(STC, crown_points_circular, CHS, cr_c)
crown_points_conn = [crown_points_circular(:,2) -1.*crown_points_circular(:,1)];
conn_points_straight = STC;
conn_points_straight(:,1) = conn_points_straight(:,1) + cr_c;
conn_points_straight(:,2) = conn_points_straight(:,2) + cr_c + 0.5*CHS;
end_connector = [crown_points_conn;conn_points_straight];
end

function [strut_bc] = strut_bc_generator(D1_corrected_strut, crown_top, crown_bottom, cr, cp_height, cr_yscale)
crown_points_bottom_shifted(:,1) = crown_bottom(:,2) + cr;
crown_points_bottom_shifted(:,2) = crown_bottom(:,1) - cp_height/2 - cr*cr_yscale;% cr_yscale accommodates for y scaling of circular points
crown_points_top_shifted(:,1) = crown_top(:,2) - cr;
crown_points_top_shifted(:,2) = crown_top(:,1) + cp_height/2 + cr*cr_yscale;
crown_points_top_shifted_fliiped = flip(crown_points_top_shifted,1);
strut_bc = [crown_points_bottom_shifted; D1_corrected_strut; crown_points_top_shifted_fliiped];
end

function [crown_points] = crown_aggregator_V3(parsec_points_left,parsec_points_right, circular_points, type, cr_yscale)
if type == "standard"
    circular_points(:,1) = circular_points(:,1).*cr_yscale;% Scaling circular points in y direction. Parsec has already been shifted accordingly
    crown_points_right = [circular_points;parsec_points_right];
    crown_points = [crown_points_right];
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

function slopes = strut_slopes(strut, np_strut)
slopes = (strut([2:np_strut],2) - strut([1:np_strut-1],2))./(strut([2:np_strut],1) - strut([1:np_strut-1],1));
end

function [corrected_strut, corrected_slope] = strut_slope_correction(strut, slopes, cutoff)
[~,top_index] = min(abs(strut(:,2) + cutoff));
[~,bottom_index] = min(abs(strut(:,2) - cutoff));
corrected_strut = strut(top_index:bottom_index,:);
corrected_slope = slopes(top_index - 1: bottom_index);
end

function [corrected_strut, corrected_slope] = strut_slope_correction_v2(strut, slopes, strut_height, cutoff)
[~,top_index] = min(abs(strut(:,2) + 0.5*strut_height - cutoff));
[~,bottom_index] = min(abs(-1.*strut(:,2) - 0.5*strut_height + cutoff));
corrected_strut = strut(top_index:bottom_index,:);
corrected_slope = slopes(top_index - 1: bottom_index);
end

function [B_bc] = DA_to_BC(B, cp_width, cp_height, np_strut)
y_corr = linspace(-1*cp_height,cp_height, np_strut);
x_corr = (cp_width/cp_height).*y_corr;
B_bc_y = B(:,2);
B_bc_x_corr = interp1(y_corr, x_corr, B_bc_y, "linear");
B_bc_x = B(:,1) + B_bc_x_corr;
B_bc = [B_bc_x B_bc_y];
end

function [SL] = slanted_strut_generator(HPQS, WPQS, np_strut)
y = linspace(-HPQS/2,HPQS/2, np_strut);
m = -(WPQS)/(HPQS);
x = m.*y;
SL = [x' y'];
end

function [S] = straight_strut_generator(HAPS, np_strut)
y = linspace(-HAPS/2,HAPS/2, np_strut);
x = 0.*y;
S = [x' y'];
end

function [strut_uc] = strut_uc_generator(SUC, crown_points_circular_st, cr_rc, HUCS)
crown_points_top = [crown_points_circular_st(:,2) crown_points_circular_st(:,1)];
crown_points_top_shifted(:,1) = crown_points_top(:,1);
crown_points_top_shifted(:,2) = crown_points_top(:,2);
SUC(:,1) = SUC(:,1) + cr_rc;
SUC(:,2) = SUC(:,2) - 0.5*HUCS - cr_rc;
crown_points_top_shifted_fliiped = flip(crown_points_top_shifted,1);
strut_uc = [SUC; crown_points_top_shifted_fliiped];
end

function [crown_er_top, crown_er_bottom] = crown_er_generator(crown_points_circular_er, HERS, cr_er)
crown_points_top_right = [crown_points_circular_er(:,2) crown_points_circular_er(:,1)];
crown_points_top_left = [-1.*crown_points_circular_er(:,2) crown_points_circular_er(:,1)];
crown_points_top_left = flip(crown_points_top_left, 1);
crown_er_top = [crown_points_top_left; crown_points_top_right];
crown_points_bottom(:,1) = crown_er_top(:,1);
crown_points_bottom(:,2) = -1.*crown_er_top(:,2);
crown_er_top(:,2) = crown_er_top(:,2) + 2*cr_er + HERS; 
crown_er_bottom = crown_points_bottom;
end

function [strut_rs] = rs_generator(SRSM1, SRSM2, SRSM3, crown_points_circular_rs, HRSM1, HRSM2, HRSM3, WRSM2, cr_rc)
crown_points_bottom = [(crown_points_circular_rs(:,2)) -1.*(crown_points_circular_rs(:,1))];
crown_points_bottom_shifted = crown_points_bottom;
SRSM1_shifted(:,1) = SRSM1(:,1) + cr_rc;
SRSM1_shifted(:,2) = SRSM1(:,2) + 0.5*HRSM1 + cr_rc;
SRSM2_shifted(:,1) = SRSM2(:,1) + cr_rc + 0.5*WRSM2;
SRSM2_shifted(:,2) = SRSM2(:,2) + HRSM1 + cr_rc + 0.5*HRSM2;
SRSM3_shifted(:,1) = SRSM3(:,1) + cr_rc + WRSM2;
SRSM3_shifted(:,2) = SRSM3(:,2) + HRSM1 + cr_rc + HRSM2 + 0.5*HRSM3;
crown_points_top = [-1.*crown_points_circular_rs(:,2) crown_points_circular_rs(:,1)];
crown_points_top_shifted(:,1) = crown_points_top(:,1) + 2*cr_rc + WRSM2;
crown_points_top_shifted(:,2) = crown_points_top(:,2) + 2*cr_rc + HRSM1 + HRSM2 + HRSM3;
crown_points_top_shifted_fliiped = flip(crown_points_top_shifted,1);
strut_rs = [crown_points_bottom_shifted; SRSM1_shifted; SRSM2_shifted; SRSM3_shifted; crown_points_top_shifted_fliiped];
end

function [strut_pq] = strut_pq_generator(D1PQ_corrected_strut, crown_topPQ, crown_bottomPQ, cr_ar, HPQS, WPQS)
crown_points_bottomPQ_shifted(:,1) = crown_bottomPQ(:,2) + cr_ar - (WPQS/2);
crown_points_bottomPQ_shifted(:,2) = crown_bottomPQ(:,1) - (HPQS/2) - cr_ar;
crown_points_topPQ_shifted(:,1) = crown_topPQ(:,2) - cr_ar + (WPQS/2);
crown_points_topPQ_shifted(:,2) = crown_topPQ(:,1) + (HPQS/2) + cr_ar;
crown_points_topPQ_shifted_fliiped = flip(crown_points_topPQ_shifted,1);
strut_pq = [crown_points_bottomPQ_shifted; D1PQ_corrected_strut; crown_points_topPQ_shifted_fliiped];
strut_pq(:,1) = -1*strut_pq(:,1);
end

function [strut_ap] = strut_ap_generator(S, crown_points_circular, cr_ar, HAPS)
crown_points_bottom = [-1.*(crown_points_circular(:,2)) -1.*(crown_points_circular(:,1))];
crown_points_bottom_shifted(:,1) = crown_points_bottom(:,1) + cr_ar;
crown_points_bottom_shifted(:,2) = crown_points_bottom(:,2) - (HAPS/2) - cr_ar;
crown_points_top = [crown_points_circular(:,2) crown_points_circular(:,1)];
crown_points_top_shifted(:,1) = crown_points_top(:,1) - cr_ar;
crown_points_top_shifted(:,2) = crown_points_top(:,2) + (HAPS/2) + cr_ar;
crown_points_top_shifted_fliiped = flip(crown_points_top_shifted,1);
strut_ap = [crown_points_bottom_shifted; S; crown_points_top_shifted_fliiped];
end

function [strut_da] = strut_da_generator(B, crown_points_circular, cr, cp_width, cp_height, cr_yscale)
crown_points_circular(:,1) = crown_points_circular(:,1)*cr_yscale;% multipying y points of the crown with a scale factor along y direction 
crown_points_bottom = [-1.*(crown_points_circular(:,2)) -1.*(crown_points_circular(:,1))];
crown_points_bottom_shifted(:,1) = crown_points_bottom(:,1) + cp_width/2 + cr;
crown_points_bottom_shifted(:,2) = crown_points_bottom(:,2) - cp_height/2 - cr*cr_yscale;
crown_points_top = [crown_points_circular(:,2) crown_points_circular(:,1)];
crown_points_top_shifted(:,1) = crown_points_top(:,1) - cp_width/2 - cr;
crown_points_top_shifted(:,2) = crown_points_top(:,2) + cp_height/2 + cr*cr_yscale;
crown_points_top_shifted_fliiped = flip(crown_points_top_shifted,1);
strut_da = [crown_points_bottom_shifted; B; crown_points_top_shifted_fliiped];
end

function [crown_top] = crown_top_circular(cr, np_crown)
x_mid_pts = 0:(cr/(np_crown)):cr;
y_mid_pts = sqrt(cr^2 - x_mid_pts.^2);
y_mid_pts = y_mid_pts - cr;
crown_top = [y_mid_pts' x_mid_pts'];
crown_top(end,:) = [];
end

function [crown_top] = crown_top_circular_v2(cr, np_crown)% Version where the last point of the crown is not discarded
x_mid_pts = 0:(cr/(np_crown)):cr;
y_mid_pts = sqrt(cr^2 - x_mid_pts.^2);
y_mid_pts = y_mid_pts - cr;
crown_top = [y_mid_pts' x_mid_pts'];
end

function [B] = NURBS_DA(ww,P,t,k,n,np_strut)
% This sends ww which is the vector of variables sought
w = [ww(1) ww(2) ww(3) ww(3) ww(2) ww(1)];
u=linspace(min(t),max(t),np_strut);
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

function [B] = NURBS_AB(ww,P,t,k,n, np_strut)
% This sends ww which is the vector of variables sought
w = [ww(1) ww(2) ww(3) ww(4) ww(4) ww(3) ww(2) ww(1)];
u=linspace(min(t),max(t),np_strut);
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

function [cleaned_spline_points] = removing_close_neibours(spline_points, distance_thresh)
cleaned_spline_points = [];
cleaned_spline_points(1,:) = spline_points(1,:);
last_point = spline_points(end,:);
for i=2:1:(size(spline_points,1) - 1)
    distance_prev = ((spline_points(i,1) - cleaned_spline_points(end,1))^2) + ((spline_points(i,2) - cleaned_spline_points(end,2))^2);
    distance_last = ((spline_points(i,1) - last_point(1,1))^2) + ((spline_points(i,2) - last_point(1,2))^2);
    if distance_prev > distance_thresh && distance_last > distance_thresh
        cleaned_spline_points(end+1,:) = spline_points(i,:);
    end
end
cleaned_spline_points(end+1,:) = last_point;
end

function [strut_points_1, strut_points_2] = divide_into_two(strut_points)% Divides the strut into two parts to allow easy offset
mid_point = round(0.5*size(strut_points,1));
strut_points_1 = strut_points([1:mid_point+1],:);
strut_points_2 = strut_points([mid_point:end],:);
end

function [strut_points_1, strut_points_2, strut_points_3] = divide_into_three(strut_points)% Divides the strut into two parts to allow easy offset
mid_point1 = round(size(strut_points,1)/3);
mid_point2 = round(2*size(strut_points,1)/3);
strut_points_1 = strut_points([1:mid_point1+1],:);
strut_points_2 = strut_points([mid_point1:mid_point2+1],:);
strut_points_3 = strut_points([mid_point2:end],:);
end