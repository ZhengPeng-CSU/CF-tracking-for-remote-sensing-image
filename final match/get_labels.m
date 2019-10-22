function [out] = get_labels(labels,target_sz,dirc)
%GET_LABELS 此处显示有关此函数的摘要
%   此处显示详细说明
sz = size(labels);
new_labels = zeros(sz);
h = floor((sz(1)+1)/2);
w = floor((sz(2)+1)/2);
a = floor(target_sz(1));
b = floor(target_sz(2));
switch dirc
    case 'left'
        labels = circshift(labels, floor(sz(1:2) / 2-0.5) );
        new_labels = circshift(labels,[0,-b]);
        new_labels(:,(end-b+1):end) = 0;
        new_labels(:,1:(w-1)) = 0;
    case 'right'
        labels = circshift(labels, floor(sz(1:2) / 2-0.5) );
        new_labels = circshift(labels,[0,b]);
        new_labels(:,1:b) = 0;
        new_labels(:,(w+1):end) = 0;
    case 'top'
        labels = circshift(labels, floor(sz(1:2) / 2-0.5) );
        new_labels = circshift(labels,[-a,0]);
        new_labels((end-a+1):end,:) = 0;
        new_labels(1:(h-1),:) = 0;
    case 'bottom'
        labels = circshift(labels, floor(sz(1:2) / 2-0.5) );
        new_labels = circshift(labels,[a,0]);
        new_labels(1:a,:) = 0;
        new_labels((h+1):end,:) = 0;

end
out = circshift(new_labels, -floor(sz(1:2) / 2-0.5) );

end

