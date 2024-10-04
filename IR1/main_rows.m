classdef main_rows
    properties
        stent_r2_config
        feature_pos_r2 
        SS 
        TS 
        C
    end

    methods
        function obj = main_rows(stent_r2_config, feature_pos_r2, SS, TS, C)
            obj.stent_r2_config = stent_r2_config;
            obj.feature_pos_r2 = feature_pos_r2;
            obj.SS = SS;
            obj.TS = TS;
            obj.C = C;
        end
    end
end
