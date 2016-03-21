PRO powerspectrum
  all=8192L
  cmb=fltarr(all+1)
	cg=fltarr(all+1)
	ci=fltarr(all+1)
  ll=lindgen(all)+1
  temp=0L
	path='/huawei/osv1/FAN/kappadata/sphere_0050/'

  filename=path+'alm_cl/h_darkmatter_cl.dat'  ;the power spectrum of dark matter particles
  openR,lun,filename,/get_lun
  readu,lun,temp
  IF temp NE all+1 THEN print,'Wrong header'
  readu,lun,cmb
  readu,lun,temp
  IF temp NE all+1 THEN print,'Wrong header'
	free_lun,lun

  filename=path+'alm_cl/h_source_cl.dat'  ;the power spectrum of source galaxies
  openR,lun,filename,/get_lun
  readu,lun,temp
  IF temp NE all+1 THEN print,'Wrong header'
  readu,lun,cg
  readu,lun,temp
  IF temp NE all+1 THEN print,'Wrong header'
	free_lun,lun

	;filename=path+'alm_cl/h_images_cl.dat'  ;the power spectrum of image galaxies
  ;openR,lun,filename,/get_lun
  ;readu,lun,temp
  ;IF temp NE all+1 THEN print,'Wrong header'
  ;readu,lun,ci
  ;readu,lun,temp
  ;IF temp NE all+1 THEN print,'Wrong header'
  ;free_lun,lun

  set_plot,'ps'
  device,filename=path+'result/powerspectrum.ps',xsize=12,ysize=12
	tvlct,[0,0,0,255],[0,0,255,0],[0,255,0,0]
  plot,ll,cmb[1:all],color=0,thick=0.5,xlog=1,ylog=1,xtitle='l',charsize=0.8
  oplot,ll,cg[1:all],color=2,thick=0.5
	;cgplot,ll,ci[1:all],color='green',thick=0.5,/overplot
	al_legend,['cmb','cg','ci'],color=['black','red','green'],box=0,charsize=0.8
  device,/close

END
