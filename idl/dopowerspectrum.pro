PRO dopowerspectrum
	Plane = 0L
	Bins = 0L
	Order = 0L
	read,Plane,prompt='Plane: '
	read,Bins,prompt='Bins: '
	read,Order,prompt='Order: '
	Nside = long(2^Order)
	Npix=long64(nside2npix(Nside))
	nn=0LL
	path=strcompress('/huawei/osv1/FAN/kappadata/sphere_00'+string(Plane)+'/',/remove)

	;Read power spectrum of dark matter
	print,'Read dark matter...'
	dpix=fltarr(Npix)
	filename=strcompress(path+'deltat_'+string(Order)+'/h_darkmatter_deltat.dat',/remove)
	openr,lun,filename,/get_lun
	readu,lun,nn
	if nn ne Npix then print,'nn='+string(nn)
	readu,lun,dpix
	readu,lun,nn
	if nn ne Npix then print,'nn='+string(nn)
	free_lun,lun
	print,'Completed!'

	print,'Computing Cmb...'
	ianafast,dpix,cl_out,/ring,/silent,/double
	cl=float(cl_out)

	print,'Write cl to file...'
	n_ele=size(cl,/n_elements)
	ll=indgen(n_ele)
	cl=cl*ll*(ll+1)/(2*!PI)
	filename=strcompress(path+'alm_cl_'+string(Order)+'/h_darkmatter_cl.dat',/remove)
	openw,lun,filename,/get_lun
	writeu,lun,n_ele
	writeu,lun,cl
	writeu,lun,n_ele
	free_lun,lun
	print,'Done!'

	;Read power spectrum of source galaxies
	print,'Read source data...'
	dpix=fltarr(Npix)
	filename=strcompress(path+'deltat_'+string(Order)+'/h_source_deltat.dat',/remove)
	openr,lun,filename,/get_lun
	readu,lun,nn
	if nn ne Npix then print,'nn='+string(nn)
	readu,lun,dpix
	readu,lun,nn
	if nn ne Npix then print,'nn='+string(nn)
	free_lun,lun
	print,'Completed!'

	print,'Computing cg...'
	ianafast,dpix,cl_out,/ring,/silent,/double
	cl=float(cl_out)

	print,'Write cg to file...'
	n_ele=size(cl,/n_elements)
	ll=indgen(n_ele)
	cl=cl*ll*(ll+1)/(2*!PI)
	filename=strcompress(path+'alm_cl_'+string(Order)+'/h_source_cl.dat',/remove)
	openw,lun,filename,/get_lun
	writeu,lun,n_ele
	writeu,lun,cl
	writeu,lun,n_ele
	free_lun,lun
	print,'Done!'

	;Read power spectrum of image galaxies
	print,'Read image data...'
	dpix_all = fltarr(Npix)
	filename=strcompress(path+'deltat_'+string(Order)+'/h_image_deltat.dat',/remove)
	openr,lun,filename,/get_lun
	readu,lun,nn
	if nn ne Npix then print,'nn='+string(nn)
	readu,lun,dpix_all
	readu,lun,nn
	if nn ne Npix then print,'nn='+string(nn)
	free_lun,lun
	dpix = fltarr(Npix,Bins)
	dpix_tmp = fltarr(Npix)
	for ii=0,Bins-1 do begin
  	filename=strcompress(path+'deltat_'+string(Order)+'/h_image_deltat'+string(ii)+'.dat',/remove)
  	print,ii
  	print,filename
  	openr,lun,filename,/get_lun
  	readu,lun,nn
  	if nn ne Npix then print,'0nn='+string(nn)
  	readu,lun,dpix_tmp
  	readu,lun,nn
  	if nn ne Npix then print,'1nn='+string(nn)
  	free_lun,lun
		dpix[*,ii] = dpix_tmp
	endfor
	print,'Completed!'

	print,'Computing ci, alm and write to file...'
	print,'all'
	filename = strcompress(path+'fits_'+string(Order)+'/h_image.fits',/remove)
	ianafast,dpix_all,cl_out,alm1_out=filename,/ring,/silent,/double
	fits2alm,index,alm_out,filename
	cl=float(cl_out)
	alm=float(alm_out)
	filename=strcompress(path+'alm_cl_'+string(Order)+'/h_image_cl.dat',/remove)
	n_ele=size(cl,/n_elements)
	ll=indgen(n_ele)
	cl=cl*ll*(ll+1)/(2*!PI)
	openw,lun,filename,/get_lun
	writeu,lun,n_ele
	writeu,lun,cl
	writeu,lun,n_ele
	free_lun,lun
	filename=strcompress(path+'alm_cl_'+string(Order)+'/h_image_alm.dat',/remove)
	n_ele=size(index,/n_elements)
	openw,lun,filename,/get_lun
	writeu,lun,n_ele
	writeu,lun,index
	writeu,lun,alm
	writeu,lun,n_ele
	free_lun,lun
	for ii=0,Bins-1 do begin
  	print,'bin'+string(ii)+string(ii)
  	pix1=dpix[*,ii]
  	filename = strcompress(path+'fits_'+string(Order)+'/h_image'+string(ii)+string(ii)+'.fits',/remove)
  	ianafast,pix1,cl_out,alm1_out=filename,/ring,/silent,/double
  	fits2alm,index,alm_out,filename
  	cl=float(cl_out)
  	alm=float(alm_out)
  	filename=strcompress(path+'alm_cl_'+string(Order)+'/h_image'+string(ii)+string(ii)+'_cl.dat',/remove)
  	n_ele=size(cl,/n_elements)
  	ll=indgen(n_ele)
  	cl=cl*ll*(ll+1)/(2*!PI)
  	openw,lun,filename,/get_lun
  	writeu,lun,n_ele
  	writeu,lun,cl
  	writeu,lun,n_ele
  	free_lun,lun
  	filename=strcompress(path+'alm_cl_'+string(Order)+'/h_image'+string(ii)+string(ii)+'_alm.dat',/remove)
  	n_ele=size(index,/n_elements)
  	openw,lun,filename,/get_lun
  	writeu,lun,n_ele
  	writeu,lun,index
  	writeu,lun,alm
  	writeu,lun,n_ele
  	free_lun,lun
  	for jj=ii+1,Bins-1 do begin
    	print,'bin'+string(ii)+string(jj)
    	pix2=dpix[*,jj]
    	ianafast,pix1,cl_out,map2_in=pix2,/ring,/silent
    	cl=float(cl_out)
    	filename=strcompress(path+'alm_cl_'+string(Order)+'/h_image'+string(ii)+string(jj)+'_cl.dat',/remove)
    	n_ele=size(cl,/n_elements)
    	ll=indgen(n_ele)
    	cl=cl*ll*(ll+1)/(2*!PI)
    	openw,lun,filename,/get_lun
    	writeu,lun,n_ele
    	writeu,lun,cl
    	writeu,lun,n_ele
    	free_lun,lun
  	endfor
	endfor
	print,'Done!'
end
