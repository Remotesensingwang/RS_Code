clear
% f=['D:\02FY3D\A202203090585162088\out\FY3D_MERSI_GBAL_L1_20190329_0440_1000M_MS_TOA.tif'];
% f=['D:\02FY3D\ɽ��Ӱ��\sdjzimg.tif'];
% [a,R]=geotiffread(f); %�ȵ���ͶӰ��Ϣ
% info=geotiffinfo(f);

% f=['E:\fy3d1920\tifout\FY3D_MERSI_GBAL_L1_20190101_0705_1000M_MS_TOA.tif']
f=['C:\Users\Wangxingtao\Downloads\TOA.tif']
data=importdata(f);
data0=data(:,:,1);


% B = transpose(data); %���IDL��û�н��о���ת�ã�����һ�� 
% s=reshape(data,[2048 2000 4]);
% m=s(:,:,1);
% M=transpose(m); 
% 