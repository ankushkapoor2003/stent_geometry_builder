function [Stent_main] = IR2_Generation(x, result_location)

%% Start Solidworks
py.sld_interface.start();
x(4) = round(x(4)); x(7) = round(x(7));

% Inputs
SL = x(1); % Stent length
Rc = x(2); % Crimped radius of the stent 
SRL = x(3); % Strut length (Excluding the crown part)
Nu = x(4); % Number of unit struts in each ring
St = x(5); % Strut thickness (including the connector thickness)
Sw = x(6); % Strut width (including the connector width)
Nr = x(7); % Number of rings in a stent
p1 = x(8); % Connector shape variable 1 (100% variationa allowed 
p2 = x(9); % Connector shape variable 2 (100%) variation allowed here as well

% Dependent variables
R = 2*pi*Rc/(4*Nu);% Crown radius
CL = (SL - (Nr)*(2*R + SRL))/(Nr - 1);

% stent generation parameters
np_strut = 500;
np_crown = 200;
distance_thresh_main = 0.02; % Primary for main helix

% Generatig the strut code 2 and 4
strut = strut_generator(SRL, np_strut);

% generationg the semicircular crown points for top
crown_points_circular = crown_top_circular(R, 4*np_crown);

% Code - 1: Crown top
crown_top = crown_aggregator(crown_points_circular,'standard');

% Code 3 Crown Bottom
crown_bottom = [-1.*crown_top(:,1) -1.*crown_top(:,2)];

% Code 5 Connector
connector = connector_generator(p1, p2, np_strut, CL, Sw);

main_row_config = main_config_generator(Nu);
main_feature_pos = feature_pos_generator(main_row_config, R, SRL);
end_row_config = end_config_generator(Nu);
end_feature_pos = feature_pos_generator(end_row_config, R, SRL);

% Building IR2 Stent
[stent_main, stent_sw_name] = stent_builder(strut, crown_top, crown_bottom, connector, SRL, R, CL, main_row_config, end_row_config, main_feature_pos, end_feature_pos, Nr, Nu, Rc, SL, Sw, St, distance_thresh_main);
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
py.sld_interface.ps_to_stl_ir2(fullpath_ps,fullpath_stl);
py.sld_interface.exit_sw();
% Stent built. 
end

%% Helper functions to generate geometries

function [stent_main, Stent_sw_name] = stent_builder(strut, crown_top, crown_bottom, connector, SRL, R, CL, main_row_config, end_row_config, main_feature_pos, end_feature_pos, Nr,Nu, Rc, SL, Sw, St, distance_thresh)
stent_shift_x = 0;
stent_shift_y = 0;
stent_main = [];

% Building main row through wrapping
cylinder_name = py.sld_interface.create_base_cylinder(Rc/1000, SL/1000);% Create the base crimped cylinder of crimped radius
py.sld_interface.create_base_cylinder_axis();% Create cylinder axis
[row, wrap_name] = row_builder(strut, crown_top, crown_bottom, connector, SRL, R, CL, main_row_config, main_feature_pos, Nr, stent_shift_x, stent_shift_y, Sw, St, 1, distance_thresh);
cylinder_name_2 = py.sld_interface.create_base_cylinder(Rc/1000, SL/1000);
stent_part = py.sld_interface.combine_subtract(cylinder_name_2,wrap_name);
current_row = py.sld_interface.circpattern(stent_part,"Axis1",Nu);
full_stent = current_row;
previous_row = current_row;
py.sld_interface.trimetric_zoom_fit();
stent_shift_x = stent_shift_x + 2*R;
stent_shift_y = stent_shift_y - 2*R - SRL - CL;
rel_y_shift = - 2*R - SRL - CL;
rel_theta_shift = 2*R/Rc;

for i = 2: Nr - 1
    main_row = row_builder(strut, crown_top, crown_bottom, connector, SRL, R, CL, main_row_config, main_feature_pos, Nr, stent_shift_x, stent_shift_y, Sw, St, 0, distance_thresh);
    stent_shift_x = stent_shift_x + 2*R;
    stent_shift_y = stent_shift_y - 2*R - SRL - CL;
    stent_main = [stent_main; main_row];
    py.sld_interface.zoom_fit();
    current_row = py.sld_interface.copy_move(previous_row, "Axis1", rel_y_shift/1000,rel_theta_shift);
    if i ~= 2
        full_stent = py.sld_interface.row_combine(previous_row, full_stent);
    end
    previous_row = current_row;
    py.sld_interface.trimetric_zoom_fit();
end
full_stent = py.sld_interface.row_combine(previous_row, full_stent);
cylinder_name = py.sld_interface.create_base_cylinder(Rc/1000, SL/1000);% Create the base crimped cylinder of crimped radius
[end_row, wrap_name] = row_builder(strut, crown_top, crown_bottom, connector, SRL, R, CL, end_row_config, end_feature_pos, Nr, stent_shift_x, stent_shift_y, Sw, St, 1, distance_thresh);
stent_main = [stent_main;end_row];
py.sld_interface.zoom_fit()
py.sld_interface.trimetric_zoom_fit()
cylinder_name_2 = py.sld_interface.create_base_cylinder(Rc/1000, SL/1000);
stent_end = py.sld_interface.combine_subtract(cylinder_name_2,wrap_name);
current_row = py.sld_interface.circpattern(stent_end,"Axis1",Nu);
Stent_sw_name = py.sld_interface.row_combine(full_stent, current_row);
end

function [row, wrap_name] = row_builder(strut, crown_top, crown_bottom, connector, SRL, R, CL, row_config, feature_pos, Nr, stent_shift_x, stent_shift_y, Sw, St, SW_flag, distance_thresh)
len = size(row_config);
featurex = [];
featurey = [];
sketch_name_index = 1;
sketch_name = strings(0);
sketch_name_conn_index = 0;
sketch_name_conn = strings(0);

for i = 1:len(2)
    if row_config(i) == 1
        strut_x = crown_top(:,2) + feature_pos(1,i) + stent_shift_x;
        strut_y = crown_top(:,1) + feature_pos(2,i) + stent_shift_y;
        featurex = [featurex; strut_x];
        featurey = [featurey; strut_y];

        if SW_flag == 1
            strut_SW = [strut_x strut_y];
            strut_close = removing_close_neibours(strut_SW, distance_thresh);
            sw_data_points = spline_data_points(strut_close);
            spline_id = py.sld_interface.create_spline(sw_data_points);
            last_point = strut_close(:,end)/1000;
        end

    elseif row_config(i) == 2
        strut_x = strut(:,1) + feature_pos(1,i) + stent_shift_x;
        strut_y = strut(:,2) + feature_pos(2,i) + stent_shift_y;
        featurex = [featurex; strut_x];
        featurey = [featurey; strut_y];
        if SW_flag == 1
            strut_SW = [strut_x strut_y];
            strut_close = removing_close_neibours(strut_SW, distance_thresh);
            sw_data_points = spline_data_points(strut_close);
            py.sld_interface.create_spline(sw_data_points);
            first_point = strut_close(:,end)/1000;
            last_point = strut_close(:,1)/1000; 
        end

     elseif row_config(i) == 3
        strut_x = crown_bottom(:,2) + feature_pos(1,i) + stent_shift_x;
        strut_y = crown_bottom(:,1) + feature_pos(2,i) + stent_shift_y;
        featurex = [featurex; strut_x];
        featurey = [featurey; strut_y];
        if SW_flag == 1
            strut_SW = [strut_x strut_y];
            strut_close = removing_close_neibours(strut_SW, distance_thresh);
            sw_data_points = spline_data_points(strut_close);
            py.sld_interface.create_spline(sw_data_points);
            first_point = strut_close(:,1)/1000;
            last_point = strut_close(:,end)/1000; 
        end

    elseif row_config(i) == 4
        strut_x = strut(:,1) + feature_pos(1,i) + stent_shift_x;
        strut_y = strut(:,2) + feature_pos(2,i) + stent_shift_y;
        featurex = [featurex; strut_x];
        featurey = [featurey; strut_y];
        if SW_flag == 1
            strut_SW = [strut_x strut_y];
            strut_close = removing_close_neibours(strut_SW, distance_thresh);
            extra_point = [strut_close(end,1) strut_close(end,2) + 0.2*distance_thresh];
            strut_close = [strut_close;extra_point];
            sw_data_points = spline_data_points(strut_close);
            py.sld_interface.create_spline(sw_data_points);
            first_point = strut_close(:,1)/1000;
            py.sld_interface.create_stent_segment_width(spline_id,Sw/1000);%Create sketch width
            activesketch = py.sld_interface.get_active_sketch_name();
            sketch_name(sketch_name_index) = activesketch;
            sketch_name_index = sketch_name_index + 1;
            py.sld_interface.end_current_sketch();% Unit feature ends here so exit previous sketch
        end

    elseif row_config(i) == 5
        strut_x = connector(:,1) + feature_pos(1,i) + stent_shift_x;
        strut_y = connector(:,2) + feature_pos(2,i) + stent_shift_y;
        featurex = [featurex; strut_x];
        featurey = [featurey; strut_y];

        if SW_flag == 1
            strut_SW = [strut_x strut_y];
            strut_close = removing_close_neibours(strut_SW, distance_thresh);
            sw_data_points = spline_data_points(strut_close);
            spline_id = py.sld_interface.create_spline(sw_data_points);
            py.sld_interface.create_stent_segment_width(spline_id,Sw/1000);%Create sketch width
            activesketch = py.sld_interface.get_active_sketch_name();
            sketch_name_conn_index = sketch_name_conn_index + 1;
            sketch_name_conn(sketch_name_conn_index) = activesketch;
            py.sld_interface.end_current_sketch();% Unit feature ends here so exit previous sketch
        end
    end
end

row = [featurex featurey];

if SW_flag == 1
    no_sketches = size(sketch_name,2);
    wrapped = 0;
    for k = 1:1:no_sketches
        if wrapped == 0
            [wrap_name1] = py.sld_interface.wrap_sketch_around_base_cylinder_v2(sketch_name(k), St/1000);
            if wrap_name1 == "NOWRAP"
                continue
            end
            if sketch_name_conn_index == 0
                wrap_name = wrap_name1;
                wrapped = 1;
            else
                [wrap_name2] = py.sld_interface.wrap_sketch_around_base_cylinder_v2(sketch_name_conn(k), St/1000);
                if wrap_name2 == "NOWRAP"
                    status = py.sld_interface.delete_feature(wrap_name1);
                    continue
                else
                    wrap_name = wrap_name2;
                    wrapped = 1;
                end
            end
        end
    end
end
end

function [feature_pos] = feature_pos_generator(config, R, SRL)
pos_x = [];
pos_y = [];
x_shift = 0;
y_shift = -R;
    for i = 1:size(config,2)
        if config(i) == 1
            x_shift = x_shift + R;
            pos_x(i) = x_shift;
            x_shift = x_shift + R;

            y_shift = y_shift + R;
            pos_y(i) = y_shift;
            y_shift = y_shift - R;

        elseif config(i) == 2
            pos_x(i) = x_shift;
            
            y_shift = y_shift - 0.5*SRL;
            pos_y(i) = y_shift;
            y_shift = y_shift - 0.5*SRL;

        elseif config(i) == 3
            x_shift = x_shift + R;
            pos_x(i) = x_shift;
            x_shift = x_shift + R;

            y_shift = y_shift - R;
            pos_y(i) = y_shift;
            y_shift = y_shift + R;

        elseif config(i) == 4
            pos_x(i) = x_shift;
            
            y_shift = y_shift + 0.5*SRL;
            pos_y(i) = y_shift;
            y_shift = y_shift + 0.5*SRL;            

        elseif config(i) == 5
            x_shift = x_shift - R;
            pos_x(i) = x_shift;
            x_shift = x_shift + R;

            y_shift = y_shift - R - SRL;
            pos_y(i) = y_shift;
            y_shift = y_shift + R + SRL;
        end
    end
feature_pos = [pos_x;pos_y];


end

function [crown_points] = crown_aggregator(circular_points, type)
    if type == "standard"
        circular_points(:,1) = circular_points(:,1);
        crown_points_right = [circular_points];
        unflipped_crown_points_left = [circular_points];
        crown_points_left = flip(unflipped_crown_points_left);
        crown_points_left(:,2) = -1*crown_points_left(:,2);
        crown_points = [crown_points_left; crown_points_right];
    end
end

function [end_config] = end_config_generator(Nu)
    end_config = [];
    for i = 1:Nu
        unit_config = [1 2 3 4];
        end_config = [end_config unit_config];
    end
end

function [main_config] = main_config_generator(Nu)
    main_config = [];
    for i = 1:Nu
        unit_config = [1 2 3 4 5];
        main_config = [main_config unit_config];
    end
end

function [strut] = strut_generator(SRL, np_strut)
    x_pts = zeros(1,np_strut+1);
    y_pts = 0:SRL/np_strut:SRL;
    y_pts = y_pts - 0.5*SRL;
    strut = [x_pts' y_pts'];
end


function [crown_top] = crown_top_circular(cr, np_crown)
    x_mid_pts = 0:(cr/(np_crown)):cr;
    y_mid_pts = sqrt(cr^2 - x_mid_pts.^2);
    y_mid_pts = y_mid_pts - cr;
    crown_top = [y_mid_pts' x_mid_pts'];
end

function [connector] = connector_generator(p1, p2, np_strut, CL, Sw)
    t = 0:1/np_strut:1;
    t1 = 0.1*((4*p2 + 3) - sqrt((4*p2 + 3)^2 - 40*p2));
    t2 = 0.1*((4*p2 + 3) + sqrt((4*p2 + 3)^2 - 40*p2));
    F1 = abs((t1^2)*((1-t1)^2)*(p2-t1));
    F2 = abs((t2^2)*((1-t2)^2)*(p2-t2));
    if F1 >= F2
        t0 = t1;
    elseif F2 > F1
        t0 = t2;
    end
    ft = (p1.*(t.^2).*((1-t).^2).*(p2 - t))./((t0^2)*((1-t0)^2)*(p2 - t0));
    x_main = t.*(CL - Sw);
    x_main = x_main + 0.5*Sw;
    x_start = 0:(0.5*Sw)/(0.5*np_strut):0.5*Sw;
    x_end = x_start + CL - 0.5*Sw;
    x_start = x_start(1:end-1);
    x_end = x_end(2:end);
    y_start = zeros(size(x_start));
    y_end = zeros(size(x_end));
    connector_start = [y_start' x_start'];
    connector_mid = [ft' x_main'];
    connector_end = [y_end' x_end'];
    connector = [connector_start;connector_mid; connector_end];
    connector(:,2) = -1*connector(:,2);
end

%% SW related functions

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

function [data_points] = spline_data_points(spline_pts) % Convert points in [x1 y1;x2 y2] format to [x1 yz1 z1 x2 y2 z2] format
    spline_pts(:,3) = 0;
    spline_pts_trans = spline_pts';
    data_points = reshape(spline_pts_trans, 1, []);
    data_points = data_points/1000;
end