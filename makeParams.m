function [params,features] = makeParams( video )
%MAKEPARAMS 此处显示有关此函数的摘要
%   此处显示详细说明
params.target_switch = false;
params.occ = 'no';
params.peak_thresh = 0;
params.offset_thresh = 0;
features.gray = false;
features.hog = false;
features.CN_Color = true;
features.HOG_CN_Color = false;
features.hog_orientations = 9;
switch video
    case {'1_Car_01'}
        params.interp_factor = 0.015;
        params.sz_factor = 1;
        params.occ = 'yes';
        params.peak_thresh = 0.1;
        params.offset_thresh = 0.08;
    case {'1_Car_10'}
        params.interp_factor = 0.005;
        params.sz_factor = 1;
        params.occ = 'yes';
        params.peak_thresh = 0.1;
        params.offset_thresh = 0.13;
    case {'1_Car_03'}
        params.interp_factor = 0.025;
        params.sz_factor = 1;
        params.target_switch = true; 
    case {'1_Car_06'}
        params.interp_factor = 0.005;
        params.sz_factor = 1; 
        params.target_switch = true;
        features.CN_Color = false;
        features.HOG_CN_Color = true;
    case {'1_Car_07','1_Car_09'}
        params.interp_factor = 0.005;
        params.sz_factor = 1;
    case {'1_Car_08','1_Car_05'}
        params.interp_factor = 0.008;
        params.sz_factor = 1;
        features.CN_Color = true;
        features.HOG_CN_Color = false;
    case {'1_Car_02','1_Car_04','1_Car_11','1_Plane_01','1_Plane_02','1_Plane_03',...
          '1_Plane_04','1_Plane_05','1_Plane_06','1_Ship_01','1_Ship_02'}
        params.interp_factor = 0.025;
        params.sz_factor = 1;       
    case {'1_Train_01'}   
        params.interp_factor = 0.025;
        params.sz_factor = 0.4;
end

end

