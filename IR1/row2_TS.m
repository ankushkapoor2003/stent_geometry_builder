classdef row2_TS
    properties
        ts_d2 % 221
        ts_tc % 222
        ts_d1 % 223
        ts_con_bc % 224
        r % crown radius
        gamma % row_1_strut height / TS_height
    end

    methods
        function obj = row2_TS(ts_d2, ts_tc, ts_d1, ts_con_bc, r, gamma)
            obj.ts_d2 = ts_d2;
            obj.ts_tc = ts_tc;
            obj.ts_d1 = ts_d1;
            obj.ts_con_bc = ts_con_bc;
            obj.r = r;
            obj.gamma = gamma;
        end
    end
end