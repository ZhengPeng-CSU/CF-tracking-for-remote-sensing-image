function [params] = makeParams( video )
%MAKEPARAMS 此处显示有关此函数的摘要
%   此处显示详细说明
params.target_switch = false;
params.sz_factor = 1;
params.occ = 'false';
params.peak_thresh = 0;
params.offset_thresh = 0;
params.interp_factor = 0.025;
params.target_switch = true; 
%'2_Car_01','2_Car_02','2_Car_03','2_Car_04','2_Car_05','2_Car_06',...
%'2_Car_07','2_Car_08','2_Car_09','2_Car_10','2_Car_11','2_Car_12','2_Car_13','2_Car_14',
%'2_Plane_01','2_Plane_02','2_Plane_03','2_Ship_01','2_Train_01','2_Train_02'
switch video
    case {'2_Car_01','2_Car_3','2_Car_4'}        %%遮挡
        params.interp_factor = 0.025;
        params.occ = 'true';
        params.peak_thresh = 0.1;%0.1
        params.offset_thresh = 0.13;%0.08
%     case {'1_Car_06'}           %%
%         params.interp_factor = 0.025;
%         params.sz_factor = 1; 
%         params.target_switch = true;  
    case {'2_Train_01'}             %%
        params.interp_factor = 0.025;
        params.sz_factor = 0.4;
end

end

