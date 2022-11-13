;coding=utf-8
;*****************************************************
;对敦煌定标场进行TOA计算（以敦煌定标场中心点为中心，计算0.1°*0.1°范围的均值）
;*****************************************************

pro fy3d_calculate_toa
  compile_opt idl2
  input_directory='H:\00data\FY3D\FY3D_dunhuang\2021-yang'
  ;input_directory='F:\FYdata\fy3d\2021'
  out_directory='H:\00data\FY3D\FY3D_dunhuang\tifout\snowwatercloud\Yang\2021\'
  ;文件日期 角度 匹配站点范围各个波段的toa均值
  ;openw,lun,'H:\00data\toa\FY3D\snowwatercloud\Yang(tbb)\2019\dh_toa2019_fy3d111.txt',/get_lun,/append,width=500

  ;敦煌定标场中心坐标
  dh_lon=94.32083333333334
  dh_lat=40.1375
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
    
    ;*****************************************************计算TOA(1-19波段)*****************************************************   
    fy3d_level1b_read,file_list_hdf[file_i_hdf],szdata=szdata,toadata,/reflectance
    
    ;*****************************************************去云处理************************************************************* 
    fy3d_cloud,file_list_hdf[file_i_hdf],toadata,clouddata
    
    
    toadata_size=size(toadata)
    cloudpos=where(clouddata ne 0)
    
    for layer_i=0,toadata_size[3]-1 do begin
      data=toadata[*,*,layer_i]
      data[cloudpos]=!VALUES.F_NAN
      ;data[cloudpos]=-100
      toadata[*,*,layer_i]=data
    endfor    
    result_tiff_name_cloud='H:\00data\FY3D\FY3D_dunhuang\tifout\cloud\'+strmid(file_basename(file_list_hdf[file_i_hdf],'.hdf'),19,13)+'_cloud.tif'
    write_tiff,result_tiff_name_cloud,clouddata,planarconfig=2,compression=1,/float
    
;    ;重采样为6KM即（6*6数组）
;    single_band_size=size(toadata[*,*,0])
;    single_col=single_band_size[1]
;    single_row=single_band_size[2]
;    data=toadata[*,*,0]
;    data_6km=MAKE_ARRAY((single_col-2)/6,(single_row-2)/6,value=!VALUES.F_NAN,/FLOAT)
;    k=0
;    
;    ;按行存储，for循环j要在前面
;    for j=0,single_row-3,6 do begin  
;      for i=0,single_col-3,6 do begin
;        value=data[i:i+5,j:j+5]
;        NotNaN_pos=WHERE(FINITE(value),/null) ;查找不是NaN的索引下标
;        if NotNaN_pos eq !null then begin
;          data_6km[k]=!VALUES.F_NAN
;        endif else begin
;          NotNaN_value=value[NotNaN_pos]
;          sort_value=NotNaN_value[sort(NotNaN_value)]
;          count_value=n_elements(sort_value)
;          if count_value gt 4 then begin
;            dark_pixels=ceil(count_value*0.2)
;            light_pixels=ceil(count_value*0.5)
;            data_6km[k]=mean(sort_value[dark_pixels:count_value-light_pixels-1])
;          endif else begin
;            data_6km[k]=!VALUES.F_NAN
;          endelse
;        endelse
;        k=k+1            
;      endfor
;    endfor


    ;*****************************************************获取以站点为中心0.1°*0.1°范围内各波段表观反射率的均值 *****************************************************
    pixs_TOAdata_mean=[]
    for band=0,toadata_size[3]-1 do begin
      ;pixs_TOAdata[*,*,band]=get_spmatching_data(TOAdata[*,*,band],3,pos_col[0],pos_line[0])
      pixs_TOAdata=toadata[*,*,band]
      pixs_TOAdata_mean=[pixs_TOAdata_mean,mean(pixs_TOAdata[pos],/nan)]
      
    endfor
    ;print,pixs_TOAdata_mean

    ;文件日期 角度 匹配站点范围各个波段的toa均值  ,逗号分隔
    data=[string(datetime),string(n_elements(pos)),string(point_degree),string(pixs_TOAdata_mean)]
    ;printf,lun,strcompress(data,/remove_all);,format='(25(a,:,","))'

    result_tiff_name=out_directory+strmid(file_basename(file_list_hdf[file_i_hdf],'.hdf'),19,13)+'_TOA_1km.tif'
    write_tiff,result_tiff_name,TOAdata,planarconfig=2,compression=1,/float   ;planarconfig=2(BSQ) 说明导入的数据是（列，行，通道数）这也是IDL的常用的，用envi打开格式为（2048 x 2000 x 4）,matlab打开格式为（2000，2048，4）
    print,file_basename(file_list_hdf[file_i_hdf])+STRCOMPRESS(string(Longdata[median(pos)]))+STRCOMPRESS(string(Latdata[median(pos)]))+string(systime(1)-starttime1)+string(file_n_hdf-file_i_hdf-1)

    szdata=!null
    sadata=!null
    vzdata=!null
    vadata=!null
    point_degree=!null
    toadata=!null
    clouddata=!null
    pixs_TOAdata=!null
    pixs_TOAdata_mean=!null

  endfor
  ;free_lun,lun
  print,'所有文件提取完成'
end