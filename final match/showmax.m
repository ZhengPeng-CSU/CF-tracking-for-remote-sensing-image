function [ ] = showmax(response,sz)
%SHOWMAX �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
labels = circshift(response, floor(sz(1:2) / 2-0.5) );
surf(labels);
end

