function [pathloss] = pathloss_cal(obj_loc_A,obj_loc_B,element_area,PLE)

    distance = norm(obj_loc_A-obj_loc_B);

    pathloss = element_area * power(distance,PLE);
end