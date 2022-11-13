;coding=utf-8
;*****************************************************
;FY3D数据进行云掩膜计算（云、水体、雪） 不是云的像元值为0
;cloud wuhan
;water Yang
;snow wuhan
;*****************************************************

;局部标准差计算中需要对图像与进行局部求和操作。
function get_std, banddata, filter_h, filter_w
  ;banddata 是原始的二维数组  ;fileter_h 是滑块的高度  ;filter_w 是滑块的宽度
  ; invalid_value 是数据中的无效值、空值等，根据需要其是否参与计算
  ; 函数返回一个与初始矩阵大小相同的，标准差矩阵
  w = n_elements(banddata[*, 0])
  h = n_elements(banddata[0, *])
  im = banddata
  im2 = im^2
  ones = replicate(1.0, w, h) ; 初始化一个全1数组(3*3区域内有效的像元数)
  ;  ids_invalid = WHERE(im gt invalid_value, count_ids_invalid)
  ;  IF count_ids_invalid GT 0 THEN ones[ids_invalid] = invalid_value
  kernel = replicate(1.0, filter_h, filter_w) ;设定卷积核的大小(和滑块大小相同)
  s = convol(im, kernel, /center, /edge_zero)  ;piexs*piexs区域内每个有效的像元的toa之和
  s2 = convol(im2, kernel, /center, /edge_zero) ;piexs*piexs区域内每个有效的像元的toa平方之和
  ns = convol(ones, kernel, /center, /edge_zero) ;piexs*piexs区域内有效的像元数
  kernel_mean=s/ns                               ;piexs*piexs区域内有效的像元的toa的均值
  std = sqrt(abs((s2 - s^2 / ns) / (ns - 1.0) ))
  mstd=std*kernel_mean*sqrt(filter_h*filter_w)
  return,{std:std,$
    mstd:mstd}
end



pro fy3d_cloud,FY3DFile,TOAdata,landcover,CloudData
  compile_opt idl2
  fy3d_level1b_read,FY3Dfile,Te_data,/temperature
  Data0047=TOAdata[*,*,0]  ;0.47um
  Data0055=TOAdata[*,*,1] ;0.55um
  Data0065=TOAdata[*,*,2]   ;0.65um
  Data0087=TOAdata[*,*,3] ;0.87um
  Data0138=TOAdata[*,*,4]   ;1.38um
  Data0164=TOAdata[*,*,5] ;1.64um
  Data0213=TOAdata[*,*,6] ;2.13um
  Data0055_result=get_std(Data0055,3,3)
  Data0138_result=get_std(Data0138,3,3)
  Data0108=Te_data[*,*,5]
  DIM = size(Data0047,/dimensions)
  NS = DIM[0]
  NL = DIM[1]
  CloudData = MAKE_ARRAY(NS,NL,VALUE=1,/BYTE) ;背景值为1
  ;设置非背景值为0
  CloudData[WHERE(Data0047 GT 0 AND Data0047 LT 1)] = 0B  ; 有效值为0
  ;*****************************************************去云处理*****************************************************     
  cloudpos=where(((Data0055_result.std ge 0.0025) and (Data0047 ge 0.4)) or ((Data0138_result.std ge 0.0025) and (Data0138 ge 0.015)),/null)
  CloudData[cloudpos]=50B
  ;*****************************************************去水体处理*****************************************************
  
  ;;水体
  NDWIData = (Data0055-Data0087) / (Data0055+Data0087)
  CloudData[WHERE((NDWIData GT 0.0) and (Data0213 lt 0.16))] = 60B
  ;CloudData[WHERE((NDWIData GT 0.0) and (Data0213 lt 0.16))] = 60B
  ;CloudData[WHERE((NDWIData GT 0.0) and (landcover ne 0))] = 70B
;  ndvi=(Data0087-Data0065)/(Data0087+Data0065)
;  waterpos=where((ndvi lt 0.1) and (Data0213 lt 0.08),/null)
;  CloudData[waterpos]=60B
  ;*****************************************************去雪处理*****************************************************  


  ndsi=(Data0087-Data0164)/(Data0087+Data0164)
  ;snowpos=where((ndsi gt 0.1) and (Data0108 lt 285),/null)
  ;snowpos=where((ndsi gt 0.1),/null)
  snowpos=where((ndsi gt 0.1) and (Data0108 lt 145),/null)
  CloudData[snowpos]=70B
  ;result_tiff_name_cloud='F:\FYdata\02watermaster\00\'+'0018.tif'
  ;write_tiff,result_tiff_name_cloud,Data0108,planarconfig=2,compression=1,/float
end