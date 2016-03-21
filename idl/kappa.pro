PRO kappa
  all=8192L
  cmb=fltarr(all+1)
  kap1=fltarr(all)
  kap2=fltarr(all)
  ll=lindgen(all)+1
  temp=0L
	path='/huawei/osv1/FAN/kappadata/sphere/'

  filename=path+'alm_cl/h_darkmatter_cl.dat'	;the power spectrum of dark matter
  openR,lun,filename,/get_lun
  readu,lun,temp
  IF temp NE all+1 THEN print,'Wrong header'
  readu,lun,cmb
  readu,lun,temp
  IF temp NE all+1 THEN print,'Wrong header'

  filename=path+'result/info_kappa1.dat'	;the kappa power spectrum recovered from magnification
  openR,lun,filename,/get_lun
  readu,lun,temp
  IF temp NE all THEN print,'Wrong header!'
  readu,lun,kap1
  readu,lun,temp
  IF temp NE all THEN print,'Wrong header!'

  filename=path+'result/info_kappa2.dat'      ;the kappa power spectrum recovered from magnification
  openR,lun,filename,/get_lun
  readu,lun,temp
  IF temp NE all THEN print,'Wrong header!'
  readu,lun,kap2
  readu,lun,temp
  IF temp NE all THEN print,'Wrong header!'

	filename=path+'halofit/pskap_revised_L500_1.142.dat'
	nlines=FILE_LINES(filename)
	data=dblarr(5,nlines-1)
	openr,lun,filename,/get_lun
	skip_lun,lun,1,/lines
	readf,lun,data,format='(5(e13.5))'
	free_lun,lun

  set_plot,'ps'
  device,filename=path+'result/powerspectrum_kappa.ps'
  cgplot,ll,cmb[1:all],linestyle=0,xlog=1,ylog=1,xtitle='l',ytitle='Pk(k)',ystyle=1,yrange=[1e-9,1e1],charsize=0.8
  cgplot,ll,kap1,color='red',linestyle=1,/overplot
  cgplot,ll,kap2,color='blue',linestyle=1,/overplot
	cgplot,data[0,*],data[0,*]*(data[0,*]+1)*data[3,*]/(2*!PI),linestyle=2,/overplot
	al_legend,['cmb','kap1','kap2','halofit'],color=['black','red','blue','black'],linestyle=[0,1,1,2],box=0,charsize=0.8
  device,/close

END
