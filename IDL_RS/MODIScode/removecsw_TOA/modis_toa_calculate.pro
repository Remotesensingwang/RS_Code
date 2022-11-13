;coding=utf-8
;*****************************************************
;对敦煌定标场进行TOA计算（以敦煌定标场中心点为中心，计算0.1°*0.1°范围的均值）
;*****************************************************


;获取数据集内容
function get_dataset,filename,dataset_name
  ; 获取文件id
  file_id = hdf_sd_start(filename)
  ;获取数据集index
  dataset_index = hdf_sd_nametoindex(file_id,dataset_name)
  ;获取数据集id
  dataset_id=hdf_sd_select(file_id,dataset_index)
  ;获取数据集的内容
  hdf_sd_getdata,dataset_id,data
  ;关闭文件
  hdf_sd_end, file_id  ; 传入文件id
  return,data
end

; 获取数据集属性内容
function get_attr, filename, dataset_name, attr_name
  ; 获取文件id
  file_id = hdf_sd_start(filename, /read)  ; 传入文件路径, 打开方式为read
  ; 获取数据集index
  dataset_index = hdf_sd_nametoindex(file_id, dataset_name)  ; 传入参数: 文件id, 数据集名称
  ; 获取数据集id
  dataset_id = hdf_sd_select(file_id, dataset_index)  ; 传入参数: 文件id, 数据集索引indedx
  ; 获取属性index
  attr_index = hdf_sd_attrfind(dataset_id, attr_name)
  ; 获取属性内容
  hdf_sd_attrinfo, dataset_id, attr_index, data=attr_data
  ; 关闭文件
  hdf_sd_end, file_id  ; 传入文件id
  ; 返回属性内容
  return, attr_data
end

pro modis_toa_calculate
  compile_opt idl2
  input_directory='H:\data\MODIS\2021'
  out_directory='H:\data\MODIS\tifout\2021\'

  ;文件日期 角度 匹配站点范围各个波段的toa均值
  openw,lun,'H:\data\toa\MODIS\dh_dingbiao2021mean_modis.txt',/get_lun,/append,width=500

  ;敦煌定标场中心坐标
  dh_lon=94.32083333333334
  dh_lat=40.1375;
  ;filename='E:\0000\MYD021KM.A2019001.0600.061.2019001192529.hdf'
  file_list_hdf=file_search(input_directory,'*.HDF',count=file_n_hdf)
  for file_i_hdf=0,file_n_hdf-1 do begin
    starttime1=systime(1)

    ;获取文件的时间、经纬度、四个角度
    datetime=strmid(file_basename(file_list_hdf[file_i_hdf],'.hdf'),10,7)+strmid(file_basename(file_list_hdf[file_i_hdf],'.hdf'),18,4)
    latdata=get_dataset(file_list_hdf[file_i_hdf],'Latitude')
    Longdata=get_dataset(file_list_hdf[file_i_hdf],'Longitude')
    ;printf,lun,strcompress(latdata,/remove_all)
    szdata=get_dataset(file_list_hdf[file_i_hdf],'SolarZenith')*0.01
    sadata=get_dataset(file_list_hdf[file_i_hdf],'SolarAzimuth')*0.01
    vzdata=get_dataset(file_list_hdf[file_i_hdf],'SensorZenith')*0.01
    vadata=get_dataset(file_list_hdf[file_i_hdf],'SensorAzimuth')*0.01


    ;获取1-19波段,26波段的DN值
    EV_250_RefSB_data=get_dataset(file_list_hdf[file_i_hdf],'EV_250_Aggr1km_RefSB') ;1-2波段
    EV_500_RefSB_data=get_dataset(file_list_hdf[file_i_hdf],'EV_500_Aggr1km_RefSB') ;3-7波段
    EV_1KM_RefSB_data=get_dataset(file_list_hdf[file_i_hdf],'EV_1KM_RefSB') ;8-19波段，26波段
    DN_band_data=[[[EV_250_RefSB_data]],[[EV_500_RefSB_data]],[[EV_1KM_RefSB_data]]]

    ;获取1-19波段,26波段的reflectance_scales值
    EV_250_RefSB_scales=get_attr(file_list_hdf[file_i_hdf],'EV_250_Aggr1km_RefSB','reflectance_scales')
    EV_500_RefSB_scales=get_attr(file_list_hdf[file_i_hdf],'EV_500_Aggr1km_RefSB','reflectance_scales')
    EV_1KM_RefSB_scales=get_attr(file_list_hdf[file_i_hdf],'EV_1KM_RefSB','reflectance_scales')
    scales_band_data=[EV_250_RefSB_scales,EV_500_RefSB_scales,EV_1KM_RefSB_scales]

    ;获取1-19波段,26波段的reflectance_offsets值
    EV_250_RefSB_offsets=get_attr(file_list_hdf[file_i_hdf],'EV_250_Aggr1km_RefSB','reflectance_offsets')
    EV_500_RefSB_offsets=get_attr(file_list_hdf[file_i_hdf],'EV_500_Aggr1km_RefSB','reflectance_offsets')
    EV_1KM_RefSB_offsets=get_attr(file_list_hdf[file_i_hdf],'EV_1KM_RefSB','reflectance_offsets')
    offsets_band_data=[EV_250_RefSB_offsets,EV_500_RefSB_offsets,EV_1KM_RefSB_offsets]


    ;获取SI数据的列、行号（有时行号为2030，有时为2040）
    DN_band_data_size=size(DN_band_data)

    ;经纬度、四个角度数据重采样 插值方法为双线性插值法（interp）
    x_size=DN_band_data_size[1]  ;列
    y_size=DN_band_data_size[2]  ;行
    Latdata=congrid(Latdata,x_size,y_size,/interp)
    Longdata=congrid(Longdata,x_size,y_size,/interp)
    szdata=congrid(szdata,x_size,y_size,/interp)
    sadata=congrid(sadata,x_size,y_size,/interp)
    vzdata=congrid(vzdata,x_size,y_size,/interp)
    vadata=congrid(vadata,x_size,y_size,/interp)

    ;以站点为中心取取0.1°*0.1范围内的数据下标
    x=0.1  
    lon_min=dh_lon-x
    lon_max=dh_lon+x
    lat_min=dh_lat-x
    lat_max=dh_lat+x

    pos=where(((Longdata ge lon_min) and (Longdata le lon_max)) and ((Latdata ge lat_min)and(Latdata le lat_max)),/null)

    ;判断pos最近点是否为！null，否则失败
    if pos eq !null then begin
      ;print,file_basename(file_list_hdf[file_i_hdf])+'pos失败'+string(systime(1)-starttime1)+string(file_n_hdf-file_i_hdf-1)
      print,file_basename(file_list_hdf[file_i_hdf])+'pos失败'
      continue
    endif

    ;print,STRCOMPRESS(string(file_basename(file_list_hdf[file_i_hdf])))+'文件的经纬度为:'+STRCOMPRESS(string(Latdata[median(pos)]))+STRCOMPRESS(string(Latdata[median(pos]))

    ;获取站点范围四个角度的均值
    point_degree=[mean(szdata[pos]),mean(sadata[pos]),mean(vzdata[pos]),mean(vadata[pos])]

    ;存储读取的1-19波段,26波段的TOA数据
    ;计算TOA   R=reflectance_scale*(DN-reflectance _offset)
    TOA_data=fltarr(DN_band_data_size[1],DN_band_data_size[2],DN_band_data_size[3])
    for layer_i=0,DN_band_data_size[3]-1 do begin
      TOA_data[*,*,layer_i]=(scales_band_data[layer_i]*(((DN_band_data[*,*,layer_i] ge 0) and (DN_band_data[*,*,layer_i] le 32767))*DN_band_data[*,*,layer_i]-offsets_band_data[layer_i]))/cos(szdata*!dtor)
      ;     print,min(TOA_data[*,*,layer_i])
      ;     print,max(TOA_data[*,*,layer_i])
    endfor

    ;获取以站点为中心0.1°*0.1°范围内各波段表观反射率的均值
    toadata_size=size(TOA_data)
    pixs_TOAdata_mean=[]
    for band=0,toadata_size[3]-1 do begin
      pixs_TOAdata=TOA_data[*,*,band]
      ;temp=pixs_TOAdata[pos]
      pixs_TOAdata_mean=[pixs_TOAdata_mean,mean(pixs_TOAdata[pos])]
    endfor
    ;print,pixs_TOAdata_mean

    ;文件日期 角度 匹配站点范围各个波段的toa均值  ,逗号分隔
    data=[string(datetime),string(n_elements(pos)),string(point_degree),string(pixs_TOAdata_mean)]
    printf,lun,strcompress(data,/remove_all);,format='(25(a,:,","))'

    result_tiff_name=out_directory+datetime+'_TOA.tif'
    write_tiff,result_tiff_name,TOA_data,planarconfig=2,/float   ;planarconfig=2(BSQ) 说明导入的数据是（列，行，通道数）这也是IDL的常用的，
    print,file_basename(file_list_hdf[file_i_hdf])+STRCOMPRESS(string(Longdata[median(pos)]))+STRCOMPRESS(string(Latdata[median(pos)]))+string(systime(1)-starttime1)+string(file_n_hdf-file_i_hdf-1)

    DN_band_data=!null
    scales_band_data=!null
    offsets_band_data=!null

    Latdata=!null
    Longdata=!null
    szdata=!null
    sadata=!null
    vzdata=!null
    vadata=!null

    point_degree=!null
    toadata=!null
    pixs_TOAdata=!null
    pixs_TOAdata_mean=!null

  endfor
  free_lun,lun
  print,'所有文件提取完成'
end