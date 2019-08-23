%
%  MSS_DCF
%
%  Written by ZhengPeng, 2019
%
%  This function takes care of setting up parameters, 
%  and interfacing with the online tracking benchmark.

function results=run_MSS_DCF_cn(seq, res_path, bSaveImage)

%default settings
features.gray = false;
features.hog = false;
features.CN_Color = true;
padding = 1.5;  %extra area surrounding the target
lambda = 1e-2;     %regularization
output_sigma_factor = 0.1; %spatial bandwidth (proportional to target)
interp_factor = 0.035;  %linear interpolation factor for adaptation
features.hog_orientations = 9;
cell_size = 1;
show_visualization = 0;
target_sz = seq.init_rect(1,[4,3]);
pos = seq.init_rect(1,[2,1]) + floor(target_sz/2);
img_files = seq.s_frames;
video_path = [];

%call tracker function with all the relevant parameters
[positions, time] = MSS_tracker(video_path, img_files, pos, target_sz, ...
    padding, lambda, output_sigma_factor, interp_factor, ...
    cell_size, features, show_visualization);
if bSaveImage
    imwrite(frame2im(getframe(gcf)),[res_path num2str(frame) '.jpg']); 
end
%return results to benchmark, in a workspace variable
rects = [positions(:,2) - target_sz(2)/2, positions(:,1) - target_sz(1)/2];
rects(:,3) = target_sz(2);
rects(:,4) = target_sz(1);
fps = numel(seq.s_frames) / time;

results.type = 'rect';
results.res = rects;
results.fps = fps;
%assignin('base', 'res', res);

end