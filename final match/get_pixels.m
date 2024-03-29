function [out] = get_pixels(im, pos, sz, theta )
%GET_PIXELS 此处显示有关此函数的摘要
%   此处显示详细说明
%   im:    输入的uint8类型图像
%   theta: 所求patch相对于原图，顺时针旋转的弧度数
%   pos:   要截取图像块在原图中的中心坐标
%   sz:    要截取图像块的大小
 
if isscalar(sz)
    sz = [sz, sz];
end
c = sz/2;
im = double(im);
[X, Y] = meshgrid(1:sz(2), 1:sz(1));
X = X - c(2);
Y = Y - c(1);
ct = cos(theta);
st = sin(theta);
Xi = pos(2) + ct*X - st*Y;
Yi = pos(1) + st*X + ct*Y;
method = 'linear';
if size(im,3)==3
    out(:,:,1) = interp2(im(:,:,1), Xi, Yi, method);
    out(:,:,2) = interp2(im(:,:,2), Xi, Yi, method);
    out(:,:,3) = interp2(im(:,:,3), Xi, Yi, method);
else
    out = interp2(im(:,:,1), Xi, Yi, method);
end
out = uint8(out);
end

