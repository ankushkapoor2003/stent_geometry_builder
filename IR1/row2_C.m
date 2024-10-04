classdef row2_C
    properties
        cs % 231 connector strut
        cc1 % 232
        cs1 % 233connector straight portion
        cc % 234
        cs2 % 235 connector straight portion
        cc2 % 236
        length % length of connector strut
        straight_length % length of straight portion of lateral connector 
        tc_r % Top crown radius
        bc_r % bottom crown radius
    end

    methods
        function obj = row2_C(cs, cc1, cs1, cc, cs2, cc2, length, straight_length, tc_r, bc_r)
            obj.cs = cs;
            obj.cc1 = cc1;
            obj.cs1 = cs1;
            obj.cc = cc;
            obj.cs2 = cs2;
            obj.cc2 = cc2;
            obj.length = length;
            obj.straight_length = straight_length;
            obj.tc_r = tc_r;
            obj.bc_r = bc_r;
        end
    end
end