function [pos] = multipos_find(img,Tdt,all_pos,target_sz,paraT)
%MUTI_FIND 此处显示有关此函数的摘要
%   此处显示详细说明
%get candidate Y
para.Lambda = paraT.lambda;
para.nT = paraT.nT;
para.Lip = paraT.Lip;
para.Maxit = paraT.Maxit;
sz = paraT.sz_T;
nT = paraT.nT;
T = Tdt.T;
Dict = Tdt.Dict;
A = Tdt.A;
Temp = Tdt.Temp;
fixT = Tdt.fixT;
rect = floor([0,0,target_sz(2),target_sz(1)]);
Y = zeros(prod(sz),size(all_pos,1))+1;   
% 从所有峰值中裁剪候选样本
for i = 1:size(all_pos,1)    
    tmp_pos = all_pos(i,:);
    rec_pos = tmp_pos - floor(target_sz/2);
    rect(1:2)=[rec_pos(2),rec_pos(1)];
    try
        imp = imcrop(img,rect);
        if (size(imp,1) ~= rect(4)+1)||(size(imp,2) ~= rect(3)+1)
            imp((rect(4)+1-size(imp,1)):end,:) = 1;
            imp(:,(rect(3)+1-size(imp,2)):end) = 1;
        end
        imp = double(imresize(imp,[sz(1),sz(2)]));
        tmp = reshape(imp, prod(sz), 1);            %crop is a vector
        [tmp,~,~] = whitening(tmp);
        Y(:,i) = tmp/norm(tmp); 
    end   
end
n_sample = size(Y,2);
% 只保留与模板欧式距离最小的候选样本
% if n_sample > 2
%     dis = zeros(n_sample,1);
%     for i =1:n_sample
%         current = repmat(Y(:,i),1,nT);
%         dis(i) = sum(sqrt(sum((current - T).*(current - T),1)/size(Y,1)));
%     end
%     [~,I] = sort(dis);
%     Y = Y(:,I(1:2));
% end
% n_sample = size(Y,2);
r_eor = zeros(n_sample,1); % minimal error bound initialization

for n = 1:n_sample 

    [c] = APGLASSOup(Temp'*Y(:,n),Dict,para);
        
    D_s = (Y(:,n) - [A(:,1:nT),fixT]*[c(1:nT); c(end)]).^2;  
    r_eor(n) = sum(D_s);             %reconstruction error
    
    if(sum(c(1:nT))<=0) %remove the inverse intensity patterns
        continue;
    end

end

t = find(r_eor == min(r_eor));
if length(t)>1
    t =floor(mean(t));
end
pos = all_pos(t,:);

if isempty(pos)
    pos = all_pos(1,:);
end

end















