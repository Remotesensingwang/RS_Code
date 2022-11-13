;coding=utf-8
;*****************************************************
;对BTH(京津冀)进行TOA计算,并重采样为6KM
;*****************************************************

pro fy3d_calculate_toa_BTH
  compile_opt idl2
  input_directory='D:\02FY3D\A202203090585162088'
  out_directory='D:\02FY3D\tif\'
  ;文件日期 角度 匹配站点范围各个波段的toa均值
  ;openw,lun,'H:\00data\toa\FY3D\snowwatercloud\Yang(tbb)\2019\dh_toa2019_fy3d111.txt',/get_lun,/append,width=500


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
;    Latdata=get_hdf5_data(file_i_geo,'/Geolocation/Latitude')
;    Longdata=get_hdf5_data(file_i_geo,'/Geolocation/Longitude')
    szdata=get_hdf5_data(file_i_geo,'/Geolocation/SolarZenith')*0.01;太阳天顶角



    ;*****************************************************计算TOA(1-19波段)*****************************************************
    fy3d_level1b_read,file_list_hdf[file_i_hdf],szdata=szdata,toadata,/reflectance

    ;*****************************************************去云处理*************************************************************
    fy3d_cloud,file_list_hdf[file_i_hdf],toadata,clouddata
    toadata_size=size(toadata)
    cloudpos=where(clouddata ne 0)


    ;去云
    for layer_i=0,toadata_size[3]-1 do begin
      data=toadata[*,*,layer_i]
      data[cloudpos]=!VALUES.F_NAN
      toadata[*,*,layer_i]=data       
    endfor
    
    
    ;重采样为6KM即（6*6数组)
    toadata_6km=MAKE_ARRAY((toadata_size[1]-2)/6,(toadata_size[2]-2)/6,toadata_size[3],value=!VALUES.F_NAN,/FLOAT)
    ;data_6km=MAKE_ARRAY((toadata_size[1]-2)/6,(toadata_size[2]-2)/6,value=!VALUES.F_NAN,/FLOAT)    
    
;    k=0L
;    help,k
    Data0065=toadata[*,*,2]
    toadata_6km_col=(toadata_size[1]-2)/6
    for layer_i=0,toadata_size[3]-1 do begin
      ;按行存储，for循环j要在前面
      time2=systime(1)
      data_layer=toadata[*,*,layer_i]
      k=0L
      for j=0,toadata_size[2]-3,6 do begin
        for i=0,toadata_size[1]-3,6 do begin
          ;获取K具体所对应的列，行号,可应用于数据的提取（截取）
         
          k_col=k mod toadata_6km_col ;k的列（类型为数组）
          k_line=k / toadata_6km_col    ;k的行（类型为数组）
          value0065=Data0065[i:i+5,j:j+5]
          NotNaN_pos=WHERE(FINITE(value0065),/null) ;查找不是NaN的索引下标
          if NotNaN_pos eq !null then begin
            ;toadata_6km[k_col,k_line,*]=!VALUES.F_NAN
            k=k+1
            continue
          endif else begin
            NotNaN_value=value0065[NotNaN_pos]
            sort_value=NotNaN_value[sort(NotNaN_value)]
            count_value=n_elements(sort_value)
            if count_value gt 4 then begin
              dark_pixels=round(count_value*0.2)
              light_pixels=round(count_value*0.5)
              value=sort_value[dark_pixels:count_value-light_pixels-1]
              ;data_6km[k]=mean(value)
              pos=where((value0065 ge value[0]) and (value0065 le value[-1]))


;              data_layer=toadata[*,*,layer_i]
              value_layer=data_layer[i:i+5,j:j+5]
              ;value_sort=value_layer[pos]
              toadata_6km[k_col,k_line,layer_i]=mean(value_layer[pos])
              ;print,'111111111'


            endif else begin
              ;toadata_6km[k_col,k_line,*]=!VALUES.F_NAN
              k=k+1
              continue
            endelse
          endelse
          k=k+1
        endfor
      endfor
      print,string(systime(1)-time2)+string(layer_i)
    endfor
;    ;按行存储，for循环j要在前面
;    for j=0,toadata_size[2]-3,6 do begin
;      for i=0,toadata_size[1]-3,6 do begin
;        ;获取K具体所对应的列，行号,可应用于数据的提取（截取）
;        toadata_6km_col=(toadata_size[1]-2)/6
;        k_col=k mod toadata_6km_col ;k的列（类型为数组）
;        k_line=k / toadata_6km_col    ;k的行（类型为数组）
;        value0065=Data0065[i:i+5,j:j+5]
;        NotNaN_pos=WHERE(FINITE(value0065),/null) ;查找不是NaN的索引下标
;        if NotNaN_pos eq !null then begin
;          toadata_6km[k_col,k_line,*]=!VALUES.F_NAN
;          ;continue
;        endif else begin
;          NotNaN_value=value0065[NotNaN_pos]
;          sort_value=NotNaN_value[sort(NotNaN_value)]
;          count_value=n_elements(sort_value)
;          if count_value gt 4 then begin
;            dark_pixels=round(count_value*0.2)
;            light_pixels=round(count_value*0.5)
;            value=sort_value[dark_pixels:count_value-light_pixels-1]
;            ;data_6km[k]=mean(value)
;            pos=where((value0065 ge value[0]) and (value0065 le value[-1]))
;            
;            for layer_i=0,toadata_size[3]-1 do begin
;              data_layer=toadata[*,*,layer_i]
;              value_layer=data_layer[i:i+5,j:j+5]
;              value_sort=value_layer[pos]
;              toadata_6km[k_col,k_line,layer_i]=mean(value_sort)
;              ;print,'111111111'    
;            endfor
;            
;          endif else begin
;            toadata_6km[k_col,k_line,*]=!VALUES.F_NAN
;          endelse
;        endelse
;        k=k+1
;      endfor
;    endfor
    ;toadata_6km[*,*,layer_i]=data_6km
    
    
    result_tiff_name=out_directory+strmid(file_basename(file_list_hdf[file_i_hdf],'.hdf'),19,13)+'_TOA_6km.tif'
    write_tiff,result_tiff_name,toadata_6km,planarconfig=2,/float   ;planarconfig=2(BSQ) 说明导入的数据是（列，行，通道数）这也是IDL的常用的，用envi打开格式为（2048 x 2000 x 4）,matlab打开格式为（2000，2048，4）
    print,file_basename(file_list_hdf[file_i_hdf])+string(systime(1)-starttime1)+string(file_n_hdf-file_i_hdf-1)

    szdata=!null
    toadata=!null
    clouddata=!null
    toadata_6km=!null
    data_6km=!null
  endfor
  ;free_lun,lun
  print,'所有文件提取完成'
end