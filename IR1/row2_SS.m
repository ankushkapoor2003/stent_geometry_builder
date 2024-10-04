classdef row2_SS
    properties
        con_ss_bc %211
        ss_d2 % 212
        ss_tc % 213
        ss_d1 % 214
        ss_ts_bc % 215
        tc_r % Top crown radius
        bc_r % bottom crown radius
        beta % SS_height / row_1_strut_height 
    end

    methods
        function obj = row2_SS(con_ss_bc, ss_d2, ss_tc, ss_d1, ss_ts_bc,tc_r, bc_r, beta)
            obj.con_ss_bc = con_ss_bc;
            obj.ss_d2 = ss_d2;
            obj.ss_tc = ss_tc;
            obj.ss_d1 = ss_d1;
            obj.ss_ts_bc = ss_ts_bc;
            obj.tc_r = tc_r;
            obj.bc_r = bc_r;
            obj.beta = beta;
        end
    end
end


