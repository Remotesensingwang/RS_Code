;coding=utf-8
;*****************************************************
;对敦煌定标场进行TOA计算（以敦煌定标场中心点为中心，计算0.1°*0.1°范围的均值）
;cloud wuhan
;water Yang
;snow Yang
;*****************************************************

;;*****************************************************读取数据集数据****************************************************
function get_hdf5_data,hd_name,filename
  file_id = H5F_OPEN(hd_name)
  dataset_id=H5D_OPEN(file_id,filename)
  data=H5D_READ(dataset_id)
  return,data
  h5d_close,dataset_id
  h5d_close,file_id
end

;*****************************************************读取数据集标签属性值****************************************************
;hd_name=文件路径名称，filename=数据集具体标签名称，attr_name=具体标签的属性名称
function get_hdf5_attr_data,hd_name,filename,attr_name
  file_id = H5F_OPEN(hd_name)
  dataset_id=H5D_OPEN(file_id,filename)
  attr_id=H5A_OPEN_Name(dataset_id,attr_name)
  data=H5A_READ(attr_id) ;获取属性值
  return,data
  h5d_close,dataset_id
  h5d_close,file_id
end




;*****************************************************查找距离站点最近的坐标值下标****************************************************
function Spatial_matching,extract_lon,extract_lat,lon,lat
  x=(lon-extract_lon)
  y=(lat-extract_lat)
  distance=sqrt(x^2+y^2)
  min_dis=min(distance)
  pos=where(distance eq min_dis)
  return,pos
end

;*****************************************************计算1-19波段的TOA数据+24波段的亮温*****************************************************
function get_TOAdata,file,szdata
  file_id=H5F_OPEN(file)
  earthsun_distance_ratio_id= H5A_OPEN_Name(file_id,'EarthSun Distance Ratio')
  earthsun_distance_ratio=H5A_READ(earthsun_distance_ratio_id);获取日地距离
  refSB_band1000m_data=get_hdf5_data(file,'/Data/EV_1KM_RefSB');获取5-19波段的DN值
  refSB_band250m_data=get_hdf5_data(file,'/Data/EV_250_Aggr.1KM_RefSB');获取1-4波段的DN值
  refSB_band_data=[[[refSB_band250m_data]],[[refSB_band1000m_data]]] ;获取1-19波段的DN值
  ;获得每个波段的slpe和Intercept属性值（1-4波段）（5-19波段）
  refSB_band250m_slope=get_hdf5_attr_data(file,'/Data/EV_250_Aggr.1KM_RefSB','Slope')
  refSB_band250m_Intercept=get_hdf5_attr_data(file,'/Data/EV_250_Aggr.1KM_RefSB','Intercept')
  refSB_band1000m_slope=get_hdf5_attr_data(file,'/Data/EV_1KM_RefSB','Slope')
  refSB_band1000m_Intercept=get_hdf5_attr_data(file,'/Data/EV_1KM_RefSB','Intercept')

  refSB_band_slope=[refSB_band250m_slope,refSB_band1000m_slope]
  refSB_band_Intercept=[refSB_band250m_Intercept,refSB_band1000m_Intercept]
  ;获取定标值（1-19波段）
  caldata=get_hdf5_data(file,'/Calibration/VIS_Cal_Coeff')
  cal_0=caldata[0,*]
  cal_1=caldata[1,*]
  cal_2=caldata[2,*]

  band_data_size=size(refSB_band_data)

  ;计算24波段（热红外）
  t_nir_data=get_hdf5_data(file,'/Data/EV_250_Aggr.1KM_Emissive');获取24波段的DN值
  t_nirband24_data=t_nir_data[*,*,0]
  
  ;获取热红外波段属性值slope和Intercept
  t_nir_slope=get_hdf5_attr_data(file,'/Data/EV_250_Aggr.1KM_Emissive','Slope')
  t_nir_Intercept=get_hdf5_attr_data(file,'/Data/EV_250_Aggr.1KM_Emissive','Intercept') 
  
  rad=t_nirband24_data*t_nir_slope[0]+t_nir_Intercept[0]
  
 ; e1=2*!pi*imsl_constant(!Planck)*(!C^2)*(10.8^5)
  c1=1.191066e-5
  c2=1.438833
  te=(c2*933.364)/(alog(1+c1*933.364^3/rad)) 
  ;(FY3MERSI地表温度反演和专题制图的MATLAB 实现)
  
  Tbb=te*1.00133-0.0734
  dd=max(Tbb)
  ;print,TBB
  ;存储读取的1-19波段的Ref与TOA数据
  band_data_ref=fltarr(band_data_size[1],band_data_size[2],band_data_size[3])
  TOA_data=fltarr(band_data_size[1],band_data_size[2],band_data_size[3])

  for layer_i=0,band_data_size[3]-1 do begin
    band_data_dn=((refSB_band_data[*,*,layer_i] ge 0) and (refSB_band_data[*,*,layer_i] le 4095))*refSB_band_data[*,*,layer_i]*refSB_band_slope[layer_i]+refSB_band_Intercept[layer_i]
    band_data_ref[*,*,layer_i]=cal_2[layer_i]*band_data_dn^2.0+cal_1[layer_i]*band_data_dn+cal_0[layer_i]
    TOA_data[*,*,layer_i]=(earthsun_distance_ratio[0]^2*band_data_ref[*,*,layer_i])/cos(szdata*!dtor)*0.01
  endfor
  return,{TOA_data:TOA_data,$
    t:Tbb}
end

;*****************************************************
;局部标准差计算中需要对图像与进行局部求和操作。
;https://blog.csdn.net/weixin_33853827/article/details/94181585?spm=1001.2101.3001.6650.2&utm_medium=distribute.pc_relevant.none-task-blog-2%7Edefault%7ECTRLIST%7ERate-2-94181585-blog-109385441.t0_edu_mlt&depth_1-utm_source=distribute.pc_relevant.none-task-blog-2%7Edefault%7ECTRLIST%7ERate-2-94181585-blog-109385441.t0_edu_mlt&utm_relevant_index=3
;*****************************************************
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

pro dunhuang_removecsw_TOA
  compile_opt idl2
  input_directory='H:\00data\FY3D\2019-1'
  out_directory='H:\00data\FY3D\tifout\snowwatercloud\Yang\2019\'
  ;文件日期 角度 匹配站点范围各个波段的toa均值
  openw,lun,'H:\00data\toa\FY3D\snowwatercloud\Yang(tbb)\2019\dh_toa2019_fy3d.txt',/get_lun,/append,width=500

  ;敦煌定标场中心坐标
  dh_lon=94.32083333333334
  dh_lat=40.1375;

  file_list_hdf=file_search(input_directory,'*_1000M_MS.HDF',count=file_n_hdf)

  ;*****************************************************文件批处理 *****************************************************
  for file_i_hdf=0,file_n_hdf-1 do begin
    starttime1=systime(1)

    ;获取文件的时间
    datetime=strmid(file_basename(file_list_hdf[file_i_hdf],'.hdf'),19,8)+strmid(file_basename(file_list_hdf[file_i_hdf],'.hdf'),28,4)

    ;获取GEO文件的经纬度及四个角度数据
    basefile_i_geo=file_basename(file_list_hdf[file_i_hdf])
    strput, basefile_i_geo, "GEO1K_MS.HDF",33 ;字符串替换
    file_i_geo= input_directory+'\'+basefile_i_geo
    Latdata=get_hdf5_data(file_i_geo,'/Geolocation/Latitude')
    Longdata=get_hdf5_data(file_i_geo,'/Geolocation/Longitude')
    ;pos=Spatial_matching(dh_lon,dh_lat,Longdata,Latdata) ;获取距离站点最近的经纬度下标
    szdata=get_hdf5_data(file_i_geo,'/Geolocation/SolarZenith')*0.01;太阳天顶角
    sadata=get_hdf5_data(file_i_geo,'/Geolocation/SolarAzimuth')*0.01;太阳方位角
    vzdata=get_hdf5_data(file_i_geo,'/Geolocation/SensorZenith')*0.01;观测天顶角
    vadata=get_hdf5_data(file_i_geo,'/Geolocation/SensorAzimuth')*0.01;观测方位角
    landcover=get_hdf5_data(file_i_geo,'/Geolocation/LandCover')   ;陆地覆盖类型
    ;print,STRCOMPRESS(string(file_basename(file_i_geo)))+'文件的经纬度为:'+STRCOMPRESS(string(Londata[pos]))+STRCOMPRESS(string(Latdata[pos]))

    ;获取每个文件的1-19波段的TOA数据
    data=get_TOAdata(file_list_hdf[file_i_hdf],szdata)
    TOAdata=data.TOA_data
    TOAdata_size=size(TOAdata)

    ;*****************************************************去云处理*****************************************************
    toa_blue=TOAdata[*,*,0]  ;0.47um
    toa_green=TOAdata[*,*,1] ;0.55um
    toa_red=TOAdata[*,*,2]   ;0.65um
    toa_nir87=TOAdata[*,*,3] ;0.87um
    toa_nir=TOAdata[*,*,4]   ;1.38um
    toa_nir164=TOAdata[*,*,5] ;1.64um
    toa_nir213=TOAdata[*,*,6] ;2.13um
    ;toa_blue_result=get_std(toa_blue,3,3)
    toa_green_result=get_std(toa_green,3,3)
    toa_nir_result=get_std(toa_nir,3,3)
    cloudpos=where(((toa_green_result.std ge 0.0025) and (toa_blue ge 0.4)) or ((toa_nir_result.std ge 0.0025) and (toa_nir ge 0.015)),/null)
    ;cloudpos=where((toa_blue gt 0.4) or ((toa_blue_result.std gt 0.0075) and (toa_blue_result.mstd gt 0.0025)) or (toa_nir gt 0.025) or (toa_nir_result.std gt 0.003) ,/null)
    ;cloudpos=where((toa_blue gt 0.4) or ((toa_blue_result.std gt 0.0075) and (toa_blue_result.mstd gt 0.0025)),/null)

    ;*****************************************************去水体、雪处理*****************************************************
    ndvi=(toa_nir87-toa_red)/(toa_nir87+toa_red)
    waterpos=where((ndvi lt 0.1) and (toa_nir213 lt 0.08),/null)
    
    tbb=data.t
    
    ndsi=(toa_nir87-toa_nir164)/(toa_nir87+toa_nir164)
    snowpos=where((ndsi gt 0.1) and (tbb lt 285),/null)
    
    ;screenpos=[cloudpos,waterpos,snowpos]

    for layer_i=0,TOAdata_size[3]-1 do begin
      data=TOAdata[*,*,layer_i]
      data[cloudpos]=!VALUES.F_NAN
      data[waterpos]=!VALUES.F_NAN
      data[snowpos]=!VALUES.F_NAN
      TOAdata[*,*,layer_i]=data
    endfor

    ;*****************************************************以站点为中心取0.1°*0.1范围内的数据下标*****************************************************
    x=0.1
    lon_min=dh_lon-x
    lon_max=dh_lon+x
    lat_min=dh_lat-x
    lat_max=dh_lat+x

    pos=where(((Longdata ge lon_min) and (Longdata le lon_max)) and ((Latdata ge lat_min)and(Latdata le lat_max)),/null)

    ;判断pos最近点是否为！null，否则失败
    if pos eq !null then begin
      ;print,file_basename(file_list_hdf[file_i_hdf])+'pos失败'+string(systime(1)-starttime1)+string(file_n_hdf-file_i_hdf-1)
      print,file_basename(file_list_hdf[file_i_hdf])+'失败'
      continue
    endif

    ;获取站点范围四个角度的均值
    point_degree=[mean(szdata[pos]),mean(sadata[pos]),mean(vzdata[pos]),mean(vadata[pos])]


    ;*****************************************************获取以站点为中心0.1°*0.1°范围内各波段表观反射率的均值 *****************************************************
    pixs_TOAdata_mean=[]
    for band=0,TOAdata_size[3]-1 do begin
      ;pixs_TOAdata[*,*,band]=get_spmatching_data(TOAdata[*,*,band],3,pos_col[0],pos_line[0])
      pixs_TOAdata=TOAdata[*,*,band]
      pixs_TOAdata_mean=[pixs_TOAdata_mean,mean(pixs_TOAdata[pos],/nan)]
    endfor
    ;print,pixs_TOAdata_mean

    ;文件日期 角度 匹配站点范围各个波段的toa均值  ,逗号分隔
    data=[string(datetime),string(n_elements(pos)),string(point_degree),string(pixs_TOAdata_mean)]
    printf,lun,strcompress(data,/remove_all);,format='(25(a,:,","))'

    result_tiff_name=out_directory+strmid(file_basename(file_list_hdf[file_i_hdf],'.hdf'),19,13)+'_TOA_removeall.tif'
    ;write_tiff,result_tiff_name,TOAdata,planarconfig=2,/float   ;planarconfig=2(BSQ) 说明导入的数据是（列，行，通道数）这也是IDL的常用的，用envi打开格式为（2048 x 2000 x 4）,matlab打开格式为（2000，2048，4）
    print,file_basename(file_list_hdf[file_i_hdf])+STRCOMPRESS(string(Longdata[median(pos)]))+STRCOMPRESS(string(Latdata[median(pos)]))+string(systime(1)-starttime1)+string(file_n_hdf-file_i_hdf-1)

    szdata=!null
    sadata=!null
    vzdata=!null
    vadata=!null
    point_degree=!null
    TOAdata=!null
    pixs_TOAdata=!null
    pixs_TOAdata_mean=!null

  endfor
  free_lun,lun
  print,'所有文件提取完成'
end

