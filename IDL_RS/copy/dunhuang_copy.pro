;读取数据集数据
function get_hdf5_data,hd_name,filename
  file_id = H5F_OPEN(hd_name)
  dataset_id=H5D_OPEN(file_id,filename)
  data=H5D_READ(dataset_id)
  return,data
  h5d_close,dataset_id
  h5d_close,file_id
end

;读取数据集标签属性值
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

;查找距离站点最近的坐标值下标
function Spatial_matching,extract_lon,extract_lat,lon,lat
  x=(lon-extract_lon)
  y=(lat-extract_lat)
  distance=sqrt(x^2+y^2)
  min_dis=min(distance)
  pos=where(distance eq min_dis)
  return,pos
end

;计算1-19波段的TOA数据值
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

  ;存储读取的1-19波段的Ref与TOA数据
  band_data_ref=fltarr(band_data_size[1],band_data_size[2],band_data_size[3])
  TOA_data=fltarr(band_data_size[1],band_data_size[2],band_data_size[3])

  for layer_i=0,band_data_size[3]-1 do begin
    band_data_dn=((refSB_band_data[*,*,layer_i] ge 0) and (refSB_band_data[*,*,layer_i] le 4095))*refSB_band_data[*,*,layer_i]*refSB_band_slope[layer_i]+refSB_band_Intercept[layer_i]
    band_data_ref[*,*,layer_i]=cal_2[layer_i]*band_data_dn^2.0+cal_1[layer_i]*band_data_dn+cal_0[layer_i]
    ;s=cos(szdata[*,*,file_i_hdf])
    ;TOA_data[*,*,layer_i]=band_data_ref[*,*,layer_i]/cos(szdata)
    TOA_data[*,*,layer_i]=(earthsun_distance_ratio[0]^2*band_data_ref[*,*,layer_i])/cos(szdata*!dtor)*0.01
    ;print,min(TOA_data[*,*,layer_i])
    ;print,max(TOA_data[*,*,layer_i])
  endfor
  return,TOA_data
end

;以站点为中心取pixs*pixs个数据（pixs为奇数）
function get_spmatching_data,arr,pixs,xloc,yloc
  iw=intarr(pixs)+1  ;3   3/2=1
  m=indgen(pixs)-pixs/2 ;3 3/2=1
  mx=m#iw
  my=iw#m
  arrsize=size(arr)
  if (xloc gt pixs/2)&&(yloc ge pixs/2)&& $
    (xloc le arrsize[1]-pixs/2)&&(yloc ge pixs/2)&& $
    (xloc gt pixs/2)&&(yloc le arrsize[2]-pixs/2)&& $
    (xloc le arrsize[1]-pixs/2)&&(yloc le arrsize[2]-pixs/2) then begin
    data=arr[xloc+mx,yloc+my]
  endif else begin
    data=!VALUES.F_NAN
  endelse
  return,data
end

pro dunhuang_copy
  compile_opt idl2
  input_directory='H:\data\FY3D\temp'
  out_directory='H:\data\FY3D\tifout\2021\'
  ;dir_test=file_test(out_directory,/directory)
  ;文件日期 角度 匹配站点范围各个波段的toa均值
  openw,lun,'F:\FYdata\fy3d\toadata\dh_dingbiao2020day1111.txt',/get_lun,/append,width=500
;  if dir_test eq 0 then begin
;    file_mkdir,out_directory
;  endif
  
  ;敦煌定标场中心坐标
  dh_lon=94.32083333333334
  dh_lat=40.1375;
  
  file_list_hdf=file_search(input_directory,'*_1000M_MS.HDF',count=file_n_hdf)
  
  for file_i_hdf=0,file_n_hdf-1 do begin
    starttime1=systime(1)
 
    ;获取文件的时间
    datetime=strmid(file_basename(file_list_hdf[file_i_hdf],'.hdf'),19,8)+strmid(file_basename(file_list_hdf[file_i_hdf],'.hdf'),28,4)
    
    ;获取GEO文件的经纬度及四个角度数据
    basefile_i_geo=file_basename(file_list_hdf[file_i_hdf])
    strput, basefile_i_geo, "GEO1K_MS.HDF",33 ;字符串替换
    file_i_geo= input_directory+'\'+basefile_i_geo
    Latdata=get_hdf5_data(file_i_geo,'/Geolocation/Latitude')
    Londata=get_hdf5_data(file_i_geo,'/Geolocation/Longitude')
    pos=Spatial_matching(dh_lon,dh_lat,Londata,Latdata) ;获取距离站点最近的经纬度下标
    szdata=get_hdf5_data(file_i_geo,'/Geolocation/SolarZenith')*0.01;太阳天顶角
    sadata=get_hdf5_data(file_i_geo,'/Geolocation/SolarAzimuth')*0.01;太阳方位角
    vzdata=get_hdf5_data(file_i_geo,'/Geolocation/SensorZenith')*0.01;观测天顶角
    vadata=get_hdf5_data(file_i_geo,'/Geolocation/SensorAzimuth')*0.01;观测方位角
    ;print,STRCOMPRESS(string(file_basename(file_i_geo)))+'文件的经纬度为:'+STRCOMPRESS(string(Londata[pos]))+STRCOMPRESS(string(Latdata[pos]))
    
    point_degree=[szdata[pos],sadata[pos],vzdata[pos],vadata[pos]];获取距离站点最近四个角度
    
    TOAdata=get_TOAdata(file_list_hdf[file_i_hdf],szdata) ;获取每个文件的1-19波段的TOA数据
    
    ;help,TOAdata
    szdatasize=size(szdata)
    
    ;获取pos具体所对应的列，行号,可应用于数据的提取（截取）
    londata_col=szdatasize[1]
    pos_col=pos mod londata_col ;pos的列（类型为数组）
    pos_line=pos / londata_col    ;pos的行（类型为数组）
    TOAdata_size=size(TOAdata)
   
    pixs_TOAdata=MAKE_ARRAY(3,3,TOAdata_size[3],/float)
    
    pixs_TOAdata_mean=[]
    
    for band=0,TOAdata_size[3]-1 do begin
      pixs_TOAdata[*,*,band]=get_spmatching_data(TOAdata[*,*,band],3,pos_col[0],pos_line[0])
      pixs_TOAdata_mean=[pixs_TOAdata_mean,mean(pixs_TOAdata[*,*,band])]
    endfor
    ;print,pixs_TOAdata_mean
    if ~FINITE(pixs_TOAdata_mean[0]) then begin ;判段第一个数组是否为NAN，若不是返回1，即为真
      print,file_basename(file_list_hdf[file_i_hdf])+'失败'+string(systime(1)-starttime1)+string(file_n_hdf-file_i_hdf-1)
      continue
    endif
      
    ;文件日期 角度 匹配站点范围各个波段的toa均值  ,逗号分隔
    data=[string(datetime),string(point_degree),string(pixs_TOAdata_mean)]   
    printf,lun,strcompress(data,/remove_all);,format='(24(a,:,","))'
    
    ;result_tiff_name=out_directory+file_basename(file_list_hdf[file_i_hdf],'.hdf')+'_TOA.tif'
    result_tiff_name=out_directory+strmid(file_basename(file_list_hdf[file_i_hdf],'.hdf'),19,13)+'_TOA.tif'
    write_tiff,result_tiff_name,TOAdata,planarconfig=2,/float   ;planarconfig=2(BSQ) 说明导入的数据是（列，行，通道数）这也是IDL的常用的，用envi打开格式为（2048 x 2000 x 4）,matlab打开格式为（2000，2048，4）
  
    szdata=!null
    sadata=!null
    vzdata=!null
    vadata=!null
    point_degree=!null
    TOAdata=!null
    pixs_TOAdata=!null
    pixs_TOAdata_mean=!null
    
    print,file_basename(file_list_hdf[file_i_hdf])+string(systime(1)-starttime1)+string(file_n_hdf-file_i_hdf-1)
  endfor
  free_lun,lun
  print,'所有文件提取完成'
end