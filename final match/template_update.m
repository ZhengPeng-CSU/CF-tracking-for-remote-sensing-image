function [Tdt_update] = template_update(img,pos,target_sz,paraT,Tdt)
%TEMPLATE_UPDATE 此处显示有关此函数的摘要
%   此处显示详细说明
sz_T =  paraT.sz_T;
nT = paraT.nT;
Tdt_update = Tdt;
T = Tdt.T;
fixT = Tdt.fixT;
rect = pos - floor(target_sz/2);
rect=[rect(2),rect(1)];
rect(3) = target_sz(2);
rect(4) = target_sz(1);
try
    imp = imcrop(img,rect);
    imp = double(imresize(imp,[sz_T(1),sz_T(2)]));      %样本缩放以确保大小一致      
    tmp = reshape(imp,[sz_T(1)*sz_T(2),1]);
    [tmp,~,~] = whitening(tmp);
    crop = tmp/norm(tmp);                     %2范数归一化
    data_num = size(T,1);
    re_xm = repmat(crop,1,nT);
    current_dis = min(sqrt((sum((re_xm - T).*(re_xm - T),1))/data_num));
    if  current_dis < 0.05 && current_dis >0
        T = [crop,T(:,1:end-1)];
    end
    dim_T	= size(T,1);	
    A		= [T,eye(dim_T)]; %data matrix is composed of T, positive trivial T.
    Temp = [A fixT];
%Temaplate Matrix
    Dict = Temp'*Temp;
    Tdt_update.T = T;
    Tdt_update.Dict = Dict;
    Tdt_update.A = A;
    Tdt_update.Temp = Temp;
end  

end

