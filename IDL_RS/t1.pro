
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
pro t1

;  print,cc
;  p = scatterplot(a,a,xtitle='FY-',ytitle='iquam SST(K)',title='btd45c_cloud',$
;    symbol = '+', SYM_COLOR = 'blue',xrange = [260,300],yrange = [-2,3],font_size=10 )
;   

;    cc=findgen(3,4,3)+1
;    ;print,cc
;    help,cc
;    pos=5
;    
;    ;获取pos具体所对应的列，行号,可应用于数据的提取（截取）
;    londata_col=3
;    pos_col=pos mod londata_col ;pos的列（类型为数组）
;    pos_line=pos / londata_col    ;pos的行（类型为数组）
;    
;    ss='20210101_0625_TOA.tif'
;    datetime=strmid(file_basename(ss,'.tif'),0,8)+strmid(file_basename(ss,'.tif'),9,4)
;    print,datetime
    
;  file='H:\data\MODIS\2021\MYD021KM.A2021001.0720.061.2021001223639.hdf'
;  result_tiff_name='H:\data\MODIS\tifout\MYD021KM.A2021001.0720.061.2021001223639_0645_cos.tif'
;  MODIS_LEVEL1B_READ,File,1,Data0645,/REFLECTANCE
;  szdata=get_dataset(file,'SolarZenith')*0.01
;  ;经纬度、四个角度数据重采样 插值方法为双线性插值法（interp）
;  DN_band_data_size=size(Data0645)
;  x_size=DN_band_data_size[1]  ;列
;  y_size=DN_band_data_size[2]  ;行
;  szdata=congrid(szdata,x_size,y_size,/interp)
;  data=Data0645/cos(szdata*!dtor)
;  write_tiff,result_tiff_name,data,planarconfig=2,/float   ;planarconfig=2(BSQ) 说明导入的数据是（列，行，通道数）这也是IDL的常用的，

  juldaytime =JULDAY(1,2,2021,2,19)
  basetime=JULDAY(12,31,2020,00,00,00)
  day_of_year=juldaytime - basetime
  print,day_of_year,format='(f0.6)'
  aa=2.0966782406903803 - 2.0965277780778706
  print,'1111'
end
;87.289502   2.0966782406903803  2.0965277780778706    0.00015044212   0.020833333
