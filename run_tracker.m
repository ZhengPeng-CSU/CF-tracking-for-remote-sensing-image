%
%  MSS Correlation Filters
%
%  Written by Joao F. Henriques, 2014
%%%  Adapted by ZhengPeng, 2019
function [precision, fps] = run_tracker(video,show_visualization, show_plots)

	%path to the videos (you'll be able to choose one with the GUI).
	base_path = 'F:\tracking match\data';

	%default settings
	if nargin < 1, video = 'choose'; end
	if nargin < 2, show_visualization = ~strcmp(video, 'all'); end
	if nargin < 3, show_plots = ~strcmp(video, 'all'); end
 
    video = choose_video(base_path);
    padding = 1.5;  %extra area surrounding the target
	lambda = 1e-4;  %regularization
	output_sigma_factor = 0.1;  %spatial bandwidth (proportional to target)
    cell_size = 1;
    [params,features] = makeParams(video);
    sz_factor = params.sz_factor;
	%parameters according to the paper. at this point we can override
	%parameters based on the chosen kernel or feature type

    [img_files, pos, target_sz, ground_truth, video_path,~] = load_video_info(base_path, video);
    target_sz = target_sz*sz_factor;
	%call tracker function with all the relevant parameters
	[positions,rect_results,time] = MSS_tracker(video_path, img_files, pos, target_sz, ...
		padding,lambda, output_sigma_factor, params, ...
		cell_size, features, show_visualization);
    target_sz = target_sz/sz_factor;
    fps = numel(img_files) / time
	if nargout > 0,
        precision = 'unknown';
    end
%     rects = [positions(:,2) - target_sz(2)/2, positions(:,1) - target_sz(1)/2];
%     rects(:,3) = target_sz(2);
%     rects(:,4) = target_sz(1);
    rects = rect_results;
    name = strcat(video,'.mat');
    save (name,'rects');
    
end
