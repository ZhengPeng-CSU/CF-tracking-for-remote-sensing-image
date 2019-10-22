function [positions,rect_results,time] = MSS_tracker(video_path, img_files, pos, target_sz, ...
	padding, lambda, output_sigma_factor, params, cell_size, ...
	features, show_visualization)

%   Joao F. Henriques, 2014
%   revision  ZhengPeng, 2019

    interp_factor = params.interp_factor;
    sz_factor =  params.sz_factor;
    occ = params.occ;
    target_switch = params.target_switch;
    %if the target is large, lower the resolution, we don't need that much
	%detail
	resize_image = (sqrt(prod(target_sz)) >= 100);  %diagonal size >= threshold
	if resize_image
		pos = floor(pos / 2);
		target_sz = floor(target_sz / 2);
    end
    original_sz = target_sz;
    c_sz = target_sz/sz_factor;
    % for color features
    temp = load('w2crs');
    w2c = temp.w2crs;
	%window size, taking padding into account
	window_sz = floor(target_sz * (1 + padding));
% 	%we could choose a size that is a power of two, for better FFT
% 	%performance. in practice it is slower, due to the larger window size.
% 	window_sz = 2 .^ nextpow2(window_sz);
    para.lambda = [0.01,0.1,0];
    para.nT		= 9;       %%% number of templates for the sparse representation
	para.sz_T  = target_sz;
    para.Lip	= 10;
    para.Maxit	= 8;
	if prod(para.sz_T) > 150
        para.sz_T = floor(target_sz/2);
    end
    sz_T = para.sz_T;
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
    base_interp = interp_factor;
    multisampling_pos = [-target_sz(1) 0; 0 -target_sz(2); target_sz(1) 0; 0 target_sz(2)];
    %丢失后8个方向重新搜索检测
	offset = [-target_sz(1) 0; 0 -target_sz(2); target_sz(1) 0; 0 target_sz(2);...
        -target_sz(1) -target_sz(2);-target_sz(1) target_sz(2);target_sz(1) -target_sz(2);target_sz(1) target_sz(2)];
    offset_peak = zeros(length(offset),0);
    offset_response = zeros([size(yf) length(offset)]);
    current_theta = 0;  %顺时针旋转为正，逆时针旋转为负
    search_theta = [0 1 2 3 -1 -2 -3];  %需要调整
    peak_value = 0;  %保存最大峰值
	%note: variables ending with 'f' are in the Fourier domain.
	time = 0;  %to calculate FPS
	positions = zeros(numel(img_files), 2);  %to calculate precision
	for frame = 1:numel(img_files)
		%load image
		im = imread([video_path img_files{frame}]); 
		if resize_image
			im = imresize(im, 0.5);   
        end
        img = im;
        if size(img,3) > 1
            img = rgb2gray(img);
        end
		tic()
        occ_frame = false;
		if frame > 1
			%obtain a subwindow for detection at the position from last
			%frame, and convert to Fourier domain (its size is unchanged)
			x_patch = get_subwindow(im, pos, window_sz);
			xf = fft2(get_features(x_patch, features, cell_size, cos_window, w2c));
            response = real(ifft2(sum(wf .* xf,3)/size(xf,3)));
            %figure(3);showmax(response,floor(window_sz / cell_size));%响应图
            switch occ
                case 'true'
                    if max(response(:)) < params.peak_thresh
                        occ_frame = true;
                        for j=1:length(offset)
                            x_patch = get_subwindow(im, pos+offset(j,:), window_sz);
                            xf = fft2(get_features(x_patch, features, cell_size, cos_window, w2c));
                            offset_response(:,:,j) = real(ifft2(sum(wf .* xf,3)/size(xf,3)));
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
                    response = real(ifft2(sum(wf .* xf,3)/size(xf,3)));
                    %figure(3);showmax(response,floor(window_sz / cell_size));
                    [vert_delta, horiz_delta] = find(response == max(response(:)), 1);
                    current_peak = max(response(:));
                    %动态更新率
                    if peak_value > 0
                        interp_factor = base_interp * current_peak/peak_value;
                    else
                        interp_factor = base_interp;
                    end
                    if current_peak > peak_value
                        peak_value = current_peak;
                    end

                    if vert_delta > size(xf,1) / 2  %wrap around to negative half-space of vertical axis
                        vert_delta = vert_delta - size(xf,1);
                    end
                    if horiz_delta > size(xf,2) / 2  %same for horizontal axis
                        horiz_delta = horiz_delta - size(xf,2);
                    end
                    pos = pos + cell_size * [vert_delta - 1, horiz_delta - 1];
                        
                    %target location is at the maximum response. 
                case 'false'
                    [vert_delta, horiz_delta] = find(response == max(response(:)), 1);
                    current_peak = max(response(:));
                    %动态更新率
                    if peak_value > 0
                        interp_factor = base_interp * current_peak/peak_value;
                    else
                        interp_factor = base_interp;
                    end
                    if current_peak > peak_value
                        peak_value = current_peak;
                    end
                    
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
		patch = get_subwindow(im, pos, window_sz);
		xf = fft2(get_features(patch, features, cell_size, cos_window, w2c));
		xf_0 = conj(xf) .* xf; 
        xf_n = zeros([size(xf) length(multisampling_pos)]);
        for i=1:length(multisampling_pos)
            % multi_sampling 
            patch = get_subwindow(im, pos+multisampling_pos(i,:), window_sz);
            xfn = fft2(get_features(patch, features, cell_size, cos_window, w2c));
            xf_n(:,:,:,i) = conj(xfn) .*xfn;
        end
        num = bsxfun(@times, conj(xf),yf) + bsxfun(@times, conj( xf_n(:,:,:,1)),yf_t)...
            + bsxfun(@times, conj( xf_n(:,:,:,2)),yf_l) + bsxfun(@times, conj( xf_n(:,:,:,3)),yf_b)...
            + bsxfun(@times, conj( xf_n(:,:,:,4)),yf_r); 
        den = xf_0 + sum(xf_n,4) + lambda;
        new_wf = num ./ den;
        
		if frame == 1  %first frame, train with a single image       
            wf = new_wf;
            [Tdt] = get_templates(img,pos,target_sz,para);      %%% first frame, build a dictionary 
		else
			%subsequent frames, interpolate model
            if ~occ_frame 
                wf = (1 - interp_factor) * wf + interp_factor * new_wf;  
            end
            if current_theta == 0
                [Tdt] = get_templates(img,pos,target_sz,para);%全部更新
                %[Tdt] = template_update(img,pos,target_sz,para,Tdt);
            end
        end
        
		%save position and timing
		positions(frame,:) = pos;
        if frame > 1 && target_switch        %旋转检测，稀疏的思想
            theta = current_theta + search_theta;
            for i=1:size(search_theta,2)
                if theta(i) == 0
                    theta_patch = imcrop(img,[pos([2,1]) - original_sz([2,1])/2, original_sz([2,1])]);
                elseif theta(i) == 90
                    target90 = imcrop(img,[pos([2,1]) - original_sz([1,2])/2, original_sz([1,2])]);
                    theta_patch = imrotate(target90,90);
                elseif theta(i) == -90
                    target90 = imcrop(img,[pos([2,1]) - original_sz([1,2])/2, original_sz([1,2])]);
                    theta_patch = imrotate(target90,-90);                    
                else
                    theta_patch = get_pixels(img,pos,floor(original_sz+[1,1]),theta(i)/180*pi);
                end
                %theta_e(i) = sum(sum(abs(theta_patch - target)));
                tmp = imresize(theta_patch,[sz_T(1),sz_T(2)]);
                tmp = double(reshape(tmp,prod(sz_T),1));
                [tmp,~,~] = whitening(tmp);
                Y(:,i) = tmp/norm(tmp); 
            end
            paraT.Lambda = para.lambda;
            paraT.nT = para.nT;
            paraT.Lip = para.Lip;
            paraT.Maxit = para.Maxit;
            nT = para.nT;
            Dict = Tdt.Dict;
            A = Tdt.A;
            Temp = Tdt.Temp;
            fixT = Tdt.fixT;

            for n = 1:size(search_theta,2)
                [c] = APGLASSOup(Temp'*Y(:,n),Dict,paraT);     
                D_s = (Y(:,n) - [A(:,1:nT),fixT]*[c(1:nT); c(end)]).^2;  
                r_eor(n) = sum(D_s);             %reconstruction error
                if(sum(c(1:nT))<=0)              %remove the inverse intensity patterns
                    continue;
                end
            end

            t = find(r_eor == min(r_eor));
            if length(t)>1
                t = floor(min(t));
            end
            current_theta = current_theta + search_theta(t);

            new_sz = get_sz(original_sz, current_theta);
            c_sz = new_sz/sz_factor;
            
        end
        
		time = time + toc();
        box = [pos([2,1]) - c_sz([2,1])/2, c_sz([2,1])];
        rect_results(frame,:) = box;
        
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

