function [ Tdt ] = get_templates(img,pos,target_sz,paraT)
%GET_SAMPLES 此处显示有关此函数的摘要
%   此处显示详细说明
sz_T = paraT.sz_T;  %模板大小
nT = paraT.nT;      %模板数量
ra = 3;                         %正样本采集区域    
[rs, cs] = ndgrid(((1:ra)- (ra+1)/2),((1:ra) - (ra+1)/2));
rect = pos - floor(target_sz/2);
rect=[rect(2),rect(1)];
rs = rect(1)+rs;
cs = rect(2)+cs;
rect(3) = target_sz(2);
rect(4) = target_sz(1);
T = zeros(sz_T(1)*sz_T(2),ra*ra);
n =1;
for i =1:ra
    for j= 1:ra
        rect(1) = rs(i,j);
        rect(2) = cs(i,j);
        imp = imcrop(img,rect);
        imp = double(imresize(imp,[sz_T(1),sz_T(2)]));      %样本缩放以确保大小一致  
        tmp = reshape(imp,[sz_T(1)*sz_T(2),1]);
        [tmp,~,~] = whitening(tmp);
        T(:,n) = tmp/norm(tmp);                     %2范数归一化
        n = n+1;   
    end
end
%generate the initial templates for the 1st frame
fixT = T(:,5)/nT;
dim_T	= size(T,1);	
A		= [T,eye(dim_T)]; %data matrix is composed of T, positive trivial T.
%Temaplate Matrix
Temp = [A fixT];
Dict = Temp'*Temp;
Tdt.T = T;
Tdt.fixT = fixT;
Tdt.Dict = Dict;
Tdt.A = A;
Tdt.Temp = Temp;

end

