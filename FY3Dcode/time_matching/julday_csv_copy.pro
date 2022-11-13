;coding=utf-8
;HDF文件读取，并获取日期+文件名
function getfiletime,input_directory
  file_list_hdf=file_search(input_directory,'*_1000M_MS.HDF',count=file_n_hdf)
  datetime=[]
  filename=[]
  for file_i_hdf=0,file_n_hdf-1 do begin
    file_name=file_basename(file_list_hdf[file_i_hdf])
    filename=[filename,file_name]
    date_time=strmid(file_basename(file_list_hdf[file_i_hdf],'.hdf'),19,13)
    ;time=strmid(file_basename(file_list_hdf[file_i_hdf],'.hdf'),28,4)
    datetime=[datetime,date_time]
  endfor
  return,{datetime:datetime,$
    filename:filename}
end

pro julday_csv_copy
  input_directory='D:\02FY3D\A202203090585162088'
  file_datetime_structs=getfiletime(input_directory)

  filenames=file_datetime_structs.filename ;获取文件名

  basetime=JULDAY(12,31,2018,00,00,00)
  Day_of_Year=[]
  filedatetime=file_datetime_structs.datetime ;获取日期

  ;计算儒略日，并与2018年12月30日00:00做比较，求出差值
  for date_i=0,n_elements(filedatetime)-1 do begin
    year=strmid(filedatetime[date_i],0,4)
    month=strmid(filedatetime[date_i],4,2)
    day=strmid(filedatetime[date_i],6,2)
    hour=strmid(filedatetime[date_i],9,2)
    minute=strmid(filedatetime[date_i],11,2)
    juldaytime =JULDAY(month,day,year,hour,minute)
    Day_of_Year=[Day_of_Year,juldaytime- basetime]
  endfor

  ;help,Day_of_Year
  ;help,filenames

  ;CSV文件导出
  out_directory='D:\02FY3D\AERONETData\'
  outfilename=out_directory+'file_juldaytimes'+'.csv'
  WRITE_CSV,outfilename,filenames,Day_of_Year,HEADER=['filename','Day_of_Year']
  Day_of_Year=!null
  filenames=!null
end