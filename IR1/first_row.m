classdef first_row
    properties
        stent_r1_config
        feature_pos 
        D1_corrected_strut 
        D2_corrected_strut
        crown_top
        crown_bottom
        strut_half_height
        cr
    end

    methods
        function obj = first_row(stent_r1_config, feature_pos,D1_corrected_strut, D2_corrected_strut, crown_top, crown_bottom, strut_half_height, cr)
            obj.stent_r1_config = stent_r1_config;
            obj.feature_pos = feature_pos;
            obj.D1_corrected_strut = D1_corrected_strut;
            obj.D2_corrected_strut = D2_corrected_strut;
            obj.crown_top = crown_top;
            obj.crown_bottom = crown_bottom;
            obj.strut_half_height = strut_half_height;
            obj.cr = cr;
        end
    end
end