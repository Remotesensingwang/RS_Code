clear
% f=['D:\02FY3D\A202203090585162088\out\FY3D_MERSI_GBAL_L1_20190329_0440_1000M_MS_TOA.tif'];
% f=['D:\02FY3D\山建影像\sdjzimg.tif'];
% [a,R]=geotiffread(f); %先导入投影信息
% info=geotiffinfo(f);

% f=['E:\fy3d1920\tifout\FY3D_MERSI_GBAL_L1_20190101_0705_1000M_MS_TOA.tif']
f=['C:\Users\Wangxingtao\Downloads\TOA.tif']
data=importdata(f);
data0=data(:,:,1);


% B = transpose(data); %如果IDL中没有进行矩阵转置，用这一行 
% s=reshape(data,[2048 2000 4]);
% m=s(:,:,1);
% M=transpose(m); 
% 