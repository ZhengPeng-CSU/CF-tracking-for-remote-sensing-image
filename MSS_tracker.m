function [positions,rect_results,time] = MSS_tracker(video_path, img_files, pos, target_sz, ...
	padding, lambda, output_sigma_factor, params, cell_size, ...
	features, show_visualization)

%   Joao F. Henriques, 2014
%   revision  ZhengPeng, 2019



    interp_factor = params.interp_factor;
    sz_factor =  params.sz_factor;
    occ = params.occ;
    %if the target is large, lower the resolution, we don't need that much
	%detail
	resize_image = (sqrt(prod(target_sz)) >= 100);  %diagonal size >= threshold
	if resize_image
		pos = floor(pos / 2);
		target_sz = floor(target_sz / 2);
    end
    c_sz = target_sz/sz_factor;
    % for color features
    temp = load('w2crs');
    w2c = temp.w2crs;
	%window size, taking padding into account
	window_sz = floor(target_sz * (1 + padding));
% 	%we could choose a size that is a power of two, for better FFT
% 	%performance. in practice it is slower, due to the larger window size.
% 	window_sz = 2 .^ nextpow2(window_sz);

	
	%create regression labels, gaussian shaped, with a bandwidth
	%proportional to target size
	output_sigma = sqrt(prod(target_sz)) * output_sigma_factor / cell_size;
	y = gaussian_shaped_labels(output_sigma, floor(window_sz / cell_size));
    yf = fft2(y);
    yf_l = fft2(get_labels(y,target_sz/cell_size,'left'));
    yf_r = fft2(get_labels(y,target_sz/cell_size,'right'));
    yf_t = fft2(get_labels(y,target_sz/cell_size,'top'));
    yf_b = fft2(get_labels(y,target_sz/cell_size,'bottom'));          %多采样样本标签
    
	%store pre-computed cosine window
	cos_window = hann(size(yf,1)) * hann(size(yf,2))';

	if show_visualization  %create video interface
		update_visualization = show_video(img_files, video_path, resize_image);
    end
    
    %8个方向重新搜索检测
	offset = [-target_sz(1) 0; 0 -target_sz(2); target_sz(1) 0; 0 target_sz(2);...
        -target_sz(1) -target_sz(2);-target_sz(1) target_sz(2);target_sz(1) -target_sz(2);target_sz(1) target_sz(2)];
    offset_peak = zeros(length(offset),0);
    offset_response = zeros([size(yf) length(offset)]);
    
	%note: variables ending with 'f' are in the Fourier domain.
	time = 0;  %to calculate FPS
	positions = zeros(numel(img_files), 2);  %to calculate precision
    rotation = false; %旋转检测
	for frame = 1:numel(img_files)
		%load image
		im = imread([video_path img_files{frame}]); 
        %联合特征使用灰度图像以减小特征维度
		if (size(im,3) > 1) && (features.HOG_CN_Color)
			im = rgb2gray(im);
		end
		if resize_image
			im = imresize(im, 0.5);   
        end
        img = im;
        if size(img,3) > 1
            img = rgb2gray(img);
        end
        
		tic()
		if frame > 1
			%obtain a subwindow for detection at the position from last
			%frame, and convert to Fourier domain (its size is unchanged)
			x_patch = get_subwindow(im, pos, window_sz);
			xf = fft2(get_features(x_patch, features, cell_size, cos_window, w2c));
            response = real(ifft2(sum(wf .* xf,3)));
            %figure(3);showmax(response,floor(window_sz / cell_size));%响应图
            switch occ
                case 'yes'
                    if max(response(:)) < params.peak_thresh
                        for j=1:length(offset)
                            x_patch = get_subwindow(im, pos+offset(j,:), window_sz);
                            xf = fft2(get_features(x_patch, features, cell_size, cos_window, w2c));
                            offset_response(:,:,j) = real(ifft2(sum(wf .* xf,3)));
                            tmp = offset_response(:,:,j);
                            offset_peak(j) = max(tmp(:));
                        end
                        t = find(offset_peak == max(offset_peak(:)));
                        tmp = offset_response(:,:,t);
                        if max(tmp(:)) > params.offset_thresh
                            pos = pos+offset(t,:);
                        end
                    end
                    x_patch = get_subwindow(im, pos, window_sz);
                    xf = fft2(get_features(x_patch, features, cell_size, cos_window, w2c));
                    response = real(ifft2(sum(wf .* xf,3)));
                    %figure(3);showmax(response,floor(window_sz / cell_size));
                    [vert_delta, horiz_delta] = find(response == max(response(:)), 1);
                    if vert_delta > size(xf,1) / 2  %wrap around to negative half-space of vertical axis
                        vert_delta = vert_delta - size(xf,1);
                    end
                    if horiz_delta > size(xf,2) / 2  %same for horizontal axis
                        horiz_delta = horiz_delta - size(xf,2);
                    end
                    pos = pos + cell_size * [vert_delta - 1, horiz_delta - 1];
                    %target location is at the maximum response. 
                case 'no'
                    [vert_delta, horiz_delta] = find(response == max(response(:)), 1);
                    if vert_delta > size(xf,1) / 2  %wrap around to negative half-space of vertical axis
                        vert_delta = vert_delta - size(xf,1);
                    end
                    if horiz_delta > size(xf,2) / 2  %same for horizontal axis
                        horiz_delta = horiz_delta - size(xf,2);
                    end
                    pos = pos + cell_size * [vert_delta - 1, horiz_delta - 1];
            end
         end
            
		

		%obtain a subwindow for training at newly estimated target position

        x_patch = get_subwindow(im, pos, window_sz);
        x_patch_l = get_subwindow(im,pos+[0,-floor(target_sz(2))],window_sz);
        x_patch_r = get_subwindow(im,pos+[0,floor(target_sz(2))],window_sz);
        x_patch_t = get_subwindow(im,pos+[-floor(target_sz(1)),0],window_sz);
        x_patch_b = get_subwindow(im,pos+[floor(target_sz(2)),0],window_sz);      %多采样样本
      
        xf = fft2(get_features(x_patch, features, cell_size, cos_window, w2c));
        xf_l = fft2(get_features(x_patch_l, features, cell_size, cos_window, w2c));
        xf_r = fft2(get_features(x_patch_r, features, cell_size, cos_window, w2c));
        xf_t = fft2(get_features(x_patch_t, features, cell_size, cos_window, w2c));
        xf_b = fft2(get_features(x_patch_b, features, cell_size, cos_window, w2c));
        new_denominator = sum(xf.*conj(xf)+xf_l.*conj(xf_l)+xf_r.*conj(xf_r)+xf_t.*conj(xf_t)+xf_b.*conj(xf_b),3)+lambda;
        new_numerator = yf.*conj(xf)+yf_l.*conj(xf_l)+yf_r.*conj(xf_r)+yf_t.*conj(xf_t)+yf_b.*conj(xf_b);
        new_wf =  new_numerator./ new_denominator;
		
		if frame == 1  %first frame, train with a single image       
            wf = new_wf;
		else
			%subsequent frames, interpolate model
            wf = (1 - interp_factor) * wf + interp_factor * new_wf;   
        end
        
		%save position and timing
		positions(frame,:) = pos;
        if params.target_switch && frame > 1        %旋转检测，模板匹配的思想
            target0 = imcrop(img,[pos([2,1]) - c_sz([2,1])/2, c_sz([2,1])]);
            target90 = imcrop(img,[pos([2,1]) - c_sz([1,2])/2, c_sz([1,2])]);
            target90 = imrotate(target90,90);
            e0 = sum(sum(abs(target0 - target)));
            e90 = sum(sum(abs(target90 - target)));
            if e90 < e0
                rotation = true;
            end
        end
  
		time = time + toc();
        if rotation
            box = [pos([2,1]) - c_sz([1,2])/2, c_sz([1,2])];
        else
            box = [pos([2,1]) - c_sz([2,1])/2, c_sz([2,1])];
        end
        rect_results(frame,:) = box;
        if frame == 1
            target = imcrop(img,rect_results(frame,:));
        end
        
		%visualization
		if show_visualization
			stop = update_visualization(frame, box);
			if stop, break, end  %user pressed Esc, stop early			
			drawnow
% 			pause(0.05)  %uncomment to run slower
		end
		
	end

	if resize_image,
		positions = positions * 2;
	end
end

