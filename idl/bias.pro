PRO bias
  bins=10L
  all=8192L
  bis1=fltarr(bins*all)
  bis2=fltarr(bins*all)
	bfl1=fltarr(all);average for ll
	bfl2=fltarr(all)
	bff1=fltarr(bins);average for flux
	bff2=fltarr(bins)
  g=fltarr(bins)
  lff=fltarr(bins)
	lfl=lindgen(all)+1L
  temp=0L
	path='/huawei/osv1/FAN/kappadata/sphere_0050/'

  filename=path+'/result/info_bias1.dat'
  openR,lun,filename,/get_lun
  readu,lun,temp
  IF temp NE bins*all THEN print,'Wrong header!'
  readu,lun,bis1
  readu,lun,temp
  IF temp NE bins*all THEN print,'Wrong header!'
	free_lun,lun

  filename=path+'result/info_bias2.dat'
  openR,lun,filename,/get_lun
  readu,lun,temp
  IF temp NE bins*all THEN print,'Wrong header!'
  readu,lun,bis2
  readu,lun,temp
  IF temp NE bins*all THEN print,'Wrong header!'
	free_lun,lun

  filename=path+'slp.dat'
  openR,lun,filename,/get_lun
  readu,lun,temp
  IF temp NE bins THEN print,'Wrong header'
  readu,lun,g
  readu,lun,temp
  IF temp NE bins THEN print,'Wrong header'
	free_lun,lun

  FOR j=0,all-1 DO BEGIN
		index=indgen(bins)
    bfl1[j] = total(bis1[index*all+j])/bins
    bfl2[j] = total(bis2[index*all+j])/bins
  ENDFOR
	FOR i=0,bins-1 DO BEGIN
    index=indgen(all)
    bff1[i] = total(bis1[i*all+index])/all
    bff2[i] = total(bis2[i*all+index])/all
  ENDFOR

  lff[0]=600.0
  lff[1]=1500.0
  lff[2]=3000.0
  lff[3]=5000.0
  lff[4]=7000.0
  lff[5]=9000.0
  lff[6]=12500.0
  lff[7]=17500.0
  lff[8]=35000.0
  lff[9]=100000.0

  x=fltarr(2)
  y=fltarr(2)
  x[0] = 1.0e2
  x[1] = 1.0e5
  y[0] = 1.0
  y[1] = 1.0

	;average for flux
	entry_device = !d.name
  set_plot,'ps'
  device,filename=path+'result/info_bias1.ps'
	tvlct,[0,0,0,255],[0,0,255,0],[0,255,0,0]
  plot,lff,smooth(g,4),xlog=1,xtitle='s',ytitle='g&b',charsize=0.8,yrange=[-2,5]
  oplot,lff,smooth(bff1,4),color=1
  oplot,lff,smooth(bff2,4),color=1,linestyle=1
	oplot,x,y,linestyle=2
	al_legend,['g','b1','b2'],linestyle=[1,2,3],box=0,charsize=0.8
  device,/close_file
	set_plot,entry_device

	;average for all
	entry_device = !d.name
  set_plot,'ps'
  device,filename=path+'result/info_bias2.ps'
	tvlct,[0,0,0,255],[0,0,255,0],[0,255,0,0]
  plot,lfl,bfl1,color=0,xlog=1,xtitle='l',ytitle='b',charsize=0.8
  oplot,lfl,bfl2,color=1
	al_legend,['b1','b2'],linestyle=0,color=['black','red'],box=0,charsize=0.8
  device,/close_file
	set_plot,entry_device

END
