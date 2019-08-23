function [ ] = showmax(response,sz)
%SHOWMAX 此处显示有关此函数的摘要
%   此处显示详细说明
labels = circshift(response, floor(sz(1:2) / 2-0.5) );
surf(labels);
end

