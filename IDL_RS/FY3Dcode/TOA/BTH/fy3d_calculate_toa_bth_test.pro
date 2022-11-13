;coding=utf-8
;*****************************************************
;对BTH进行TOA计算（以敦煌定标场中心点为中心，计算0.1°*0.1°范围的均值）
;*****************************************************

pro fy3d_calculate_toa_bth_test
  compile_opt idl2
  input_directory='H:\bth2021\3'
  ;out_directory='H:\00data\FY3D\FY3D_dunhuang\tifout\snowwatercloud\Yang\2021\'

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
    szdata=get_hdf5_data(file_i_geo,'/Geolocation/SolarZenith')*0.01;太阳天顶角
;    sadata=get_hdf5_data(file_i_geo,'/Geolocation/SolarAzimuth')*0.01;太阳方位角
;    vzdata=get_hdf5_data(file_i_geo,'/Geolocation/SensorZenith')*0.01;观测天顶角
;    vadata=get_hdf5_data(file_i_geo,'/Geolocation/SensorAzimuth')*0.01;观测方位角
    landcover=get_hdf5_data(file_i_geo,'/Geolocation/LandCover')   ;陆地覆盖类型

    ;*****************************************************计算TOA(1-19波段)*****************************************************
    fy3d_level1b_read,file_list_hdf[file_i_hdf],szdata=szdata,toadata,/reflectance
    
    data1=[[[toadata]],[[Longdata]],[[Latdata]]]

    result_tiff_name_clouds='H:\dtcloud\01toabth\1\pro\'+strmid(file_basename(file_list_hdf[file_i_hdf],'.hdf'),19,13)+'_TOA.tif'
    ;write_tiff,result_tiff_name_clouds,data1,planarconfig=2,compression=1,/float
    
    ;*****************************************************去云处理*************************************************************
    fy3d_cloud_test,file_list_hdf[file_i_hdf],toadata,landcover,clouddata

    result_tiff_name_cloud='H:\dtcloud\02cloudmaster\1\pro\'+strmid(file_basename(file_list_hdf[file_i_hdf],'.hdf'),19,13)+'_cloud047ge030.tif'
    ;write_tiff,result_tiff_name_cloud,clouddata,planarconfig=2,compression=1,/float

    toadata_size=size(toadata)
    cloudpos=where(clouddata ne 0)

    for layer_i=0,toadata_size[3]-1 do begin
      data=toadata[*,*,layer_i]
      data[cloudpos]=!VALUES.F_NAN
      ;data[cloudpos]=-100
      toadata[*,*,layer_i]=data
    endfor
    data2=[[[toadata]],[[Longdata]],[[Latdata]]]
    result_tiff_name='H:\dtcloud\03cloudrbth\3\pro\'+strmid(file_basename(file_list_hdf[file_i_hdf],'.hdf'),19,13)+'_TOA_removecloud047ge024.tif'
    write_tiff,result_tiff_name,data2,planarconfig=2,compression=1,/float   ;planarconfig=2(BSQ) 说明导入的数据是（列，行，通道数）这也是IDL的常用的，用envi打开格式为（2048 x 2000 x 4）,matlab打开格式为（2000，2048，4）
   
    print,file_basename(file_list_hdf[file_i_hdf])+string(systime(1)-starttime1)+string(file_n_hdf-file_i_hdf-1)

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
  print,'所有文件提取完成'
end