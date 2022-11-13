pro test
   cc=findgen(14,13)
   ;print,cc
   print,cc.TYPENAME
   ;重采样为6KM即（6*6数组）
;   single_band_size=size(cc)
;   single_col=single_band_size[1]
;   single_row=single_band_size[2]
;   data=cc
;   for i=0,single_col-3,6 do begin
;     for j=0,single_row-2,6 do begin
;       value=data[i:i+5,j:j+5]
;       print,'123'
;     endfor
;   endfor
  data=findgen(20,14)+1
  data1=findgen(20,14)+10000
  s=[[[data]],[[data1]]]
  help,s
  
  print,data1
  help,data1
  print,'111'
end