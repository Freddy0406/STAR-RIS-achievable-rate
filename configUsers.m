function [loc_TUE1,loc_TUE2,loc_RUE1,loc_RUE2] = configUsers(loc_RIS,d_0)


    loc_TUE1 = [(loc_RIS(1)-0.5*d_0) , (loc_RIS(2)+0.5*d_0)];
    loc_TUE2 = [(loc_RIS(1)+0.5*d_0) , (loc_RIS(2)+0.5*d_0)];

    loc_RUE1 = [(loc_RIS(1)-0.5*d_0) , (loc_RIS(2)-0.5*d_0)];
    loc_RUE2 = [(loc_RIS(1)+0.5*d_0) , (loc_RIS(2)-0.5*d_0)];

end