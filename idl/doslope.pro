FORWARD_FUNCTION doslopefunc

PRO doslope
;Read the images from file...
path='/huawei/osv1/FAN/kappadata/sphere_0050/'
openr,lun,path+'images_flux.bin',/get_lun
TotNumGals=0LL
readu,lun,TotNumGals
flux=fltarr(TotNumGals)
readu,lun,flux
readu,lun,TotNumGals
free_lun,lun

print,'Flux:',min(flux),max(flux)

;Minus the residuals used for division of flux
flux_res = floor(min(flux)/100.0)*100.0
flux -= flux_res

;These info used for images division of flux
bins=10L
flxp = fltarr(bins+1)
flxp[0] = 200.0
flxp[1] = 1000.0
flxp[2] = 2000.0
flxp[3] = 4000.0
flxp[4] = 6000.0
flxp[5] = 8000.0
flxp[6] = 10000.0
flxp[7] = 15000.0
flxp[8] = 20000.0
flxp[9] = 50000.0
flxp[10] = max(flux)+1

;This is some information about dividing the flux bins
;flux_div_point
fdp = [0, 1.0e4, 1.0e6, 1.0e7]	;for 50
;fdp = [0, 1.0e5, 1.0e7, 1.0e8] ;for 24
n_fdp = size(fdp, /n_elements)

;flux_bin_width
fbw = [1.0e3, 1.0e4, 1.0e5, 1.0e6] ;for 50
;fbw = [1.0e3, 1.0e5, 1.0e6, 1.0e7] ;for 24

;flux_bin_need_to_average
fbna = fltarr(n_fdp)
fbna[0] = 1.0
for i=1,n_fdp-1 do begin
	fbna[i] = fbna[0]*fbw[i]/fbw[0]
endfor

;flux_bin_between_div_point
fbbdp = lonarr(n_fdp)
for i=0,n_fdp-1 do begin
  if i eq n_fdp-1 then begin
    fbbdp[i] = long((max(flux)-fdp[i])/fbw[i])
		if fdp[i]+fbw[i]*fbbdp[i] ne max(flux) then fbbdp[i]++
  endif else begin
    fbbdp[i] = long((fdp[i+1]-fdp[i])/fbw[i])
  endelse
endfor

;sum of the flux_bin_between_div_point
s_fbbdp = lonarr(n_fdp)
s_fbbdp[0] = fbbdp[0]
for i=1,n_fdp-1 do s_fbbdp[i] = s_fbbdp[i-1]+fbbdp[i]

fbin = long(total(fbbdp))
print,'fbin:',fbin

dex=intarr(TotNumGals)
for i=0,TotNumGals-1 do begin
  if flux[i] lt fdp[1] then begin
    dex[i] = floor(flux[i]/fbw[0])
  endif else if flux[i] lt fdp[2] then begin
    dex[i] = s_fbbdp[0]+floor((flux[i]-fdp[1])/fbw[1])
  endif else if flux[i] lt fdp[3] then begin
    dex[i] = s_fbbdp[1]+floor((flux[i]-fdp[2])/fbw[2])
  endif else begin
    dex[i] = s_fbbdp[2]+floor((flux[i]-fdp[3])/fbw[3])
  endelse
endfor

num=lon64arr(fbin)
num[dex]++

ctt = 0LL
for i=0,TotNumGals-1 do begin
	if flux[i] lt 200.0 then ctt++
endfor
print,'Number details in flux bin:'
print,'0 ~ 200.0',ctt
print,fdp[0],' ~ ',fdp[1],'divided by ',fbbdp[0],':',total(num[0:s_fbbdp[0]-1])
print,num[0:s_fbbdp[0]-1]
print,fdp[1],' ~ ',fdp[2],'divided by ',fbbdp[1],':',total(num[s_fbbdp[0]:s_fbbdp[1]-1])
print,num[s_fbbdp[0]:s_fbbdp[1]-1]
print,fdp[2],' ~ ',fdp[3],'divided by ',fbbdp[2],':',total(num[s_fbbdp[1]:s_fbbdp[2]-1])
print,num[s_fbbdp[1]:s_fbbdp[2]-1]
print,fdp[3],' ~ ',max(flux),'divided by ',fbbdp[3],':',total(num[s_fbbdp[2]:s_fbbdp[3]-1])
print,num[s_fbbdp[2]:s_fbbdp[3]-1]

cnt=fltarr(fbin);number of images in each flux division
for i=0,fbin-1 do begin
  if i lt s_fbbdp[0] then begin
    cnt[i] = num[i]/fbna[0]
  endif else if i lt s_fbbdp[1] then begin
    cnt[i] = num[i]/fbna[1]
  endif else if i lt s_fbbdp[2] then begin
    cnt[i] = num[i]/fbna[2]
  endif else begin
    cnt[i] = num[i]/fbna[3]
  endelse
endfor

flx=fltarr(fbin);the average flux of each flux division
for i=0,fbin-1 do begin
  if i lt s_fbbdp[0] then begin
    flx[i] = flux_res+fdp[0]+i*fbw[0]+0.5*fbw[0]
  endif else if i lt s_fbbdp[1] then begin
    flx[i] = flux_res+fdp[1]+(i-s_fbbdp[0])*fbw[1]+0.5*fbw[1]
  endif else if i lt s_fbbdp[2] then begin
    flx[i] = flux_res+fdp[2]+(i-s_fbbdp[1])*fbw[2]+0.5*fbw[2]
  endif else begin
    flx[i] = flux_res+fdp[3]+(i-s_fbbdp[2])*fbw[3]+0.5*fbw[3]
  endelse
endfor

set_plot,'ps'
device,filename=path+'flx_cnt.ps'
cgplot,flx,cnt,xstyle=1,xrange=[0,1.0e5],xtitle='flx',ytitle='cnt',charsize=0.8
device,/close

for i=0,fbin-2 do begin
  if num[i] eq 0 then num[i] = (num[i-1]+num[i+1])/2.0
  if cnt[i] eq 0 then cnt[i] = (cnt[i-1]+cnt[i+1])/2.0
endfor
cntl=alog(cnt);log of the cnt
flxl=alog(flx);log of the flx

set_plot,'ps'
device,filename=path+'flxl_cntl.ps'
cgplot,flxl,cntl,xtitle='Log(flux)',ytitle='Log(num)',charsize=0.8
device,/close

;calculate the slope of the num-flx relation
n_slp = 5
x=fltarr(n_slp)
y=fltarr(n_slp)
slpx=fltarr(fbin-n_slp+1)
slpy=fltarr(fbin-n_slp+1)
num_slp = 0 ; the number for slpx and slpy
do_slp = 1  ; check whether we should computer the slope continuely
i = 0       ; loop for the max slpx and slpy
while i lt fbin-n_slp+1 do begin
  j = 0     ; loop for the max of x and y
  k = 0     ; number of the x and y
  if k eq 0 and num_slp gt 0 then begin
    if flxl[i]-x[0] lt 0.1 then begin
      i++
      continue
    endif
  endif
  while k lt n_slp do begin
    if k eq 0 then begin
      x[k] = flxl[i+j]
      y[k] = cntl[i+j]
      slpx[num_slp] += flxl[i+j]
      j++
      k++
    endif else begin
      if i+j eq fbin then begin
        if k lt 2 then do_slp = 0
        break
      endif
      if flxl[i+j]-x[k-1] ge 0.1 then begin
        x[k] = flxl[i+j]
        y[k] = cntl[i+j]
        slpx[num_slp] += flxl[i+j]
        j++
        k++
      endif else begin
        j++
      endelse
    endelse
  endwhile
  if do_slp eq 0 then break
  slpx[num_slp] /= k
  slpy[num_slp] = doslopefunc(k,x,y)
  num_slp++
  i++
endwhile

;calculate the parameter g of flux
numg = 0L
for i=0,num_slp-1 do begin
  if slpx[i] gt 14.5 then begin
    if slpy[i] gt slpy[numg-1] then continue
  endif
  slpx[j] = slpx[i]
  slpy[j] = slpy[i]
  numg++
endfor
g = 2*(-slpy[0:numg-1]-1-1)
ff = exp(slpx[0:numg-1])

set_plot,'ps'
device,filename=path+'flxl_slp.ps'
cgplot,flxl,cntl,ystyle=1,yrange=[-20,30],xtitle='Log(flux)',charsize=0.8
cgplot,slpx[0:numg-1],slpy[0:numg-1],psym=-7,/overplot
cgplot,slpx[0:numg-1],g,psym=-5,/overplot
al_legend,['num','slope','g'],psym=[0,-7,-5],box=0,charsize=0.8
device,/close

set_plot,'ps'
device,filename=path+'flx_slp.ps'
cgplot,ff,g,psym=-2,/xlog,ystyle=1,yrange=[-4,30],xtitle='flux',ytitle='g',charsize=0.8
device,/close

openw,lun,path+'gdata.dat',/get_lun
writeu,lun,j
writeu,lun,ff
writeu,lun,g
writeu,lun,j
free_lun,lun

;calculate the parameter g of flux in each flux bin
nump = intarr(bins)
s_g = fltarr(bins)
for ii=0,bins-1 do begin
	for jj=0,numg-1 do begin
		if ff[jj] gt flxp[ii] and ff[jj] le flxp[ii+1] then begin
			nump[ii]++
			s_g[ii] += g[jj]
		endif
	endfor
endfor
slpr = fltarr(bins)
need_ex = 0
for ii=0,bins-1 do begin
	if nump[ii] ne 0 then begin
		slpr[ii] = s_g[ii]/nump[ii]
	endif else begin
		for jj=0,fbin-1 do begin
			if ff[jj] le flxp[ii] then begin
				nump[ii] = 1
				s_g[ii] = g[jj]
			endif
			if ff[jj] gt flxp[ii+1] then begin
        nump[ii]++
        s_g[ii] += g[jj]
				if nump[ii] eq 2 then begin
					break
        endif else begin
					need_ex = 1
				endelse
      endif
		endfor
		if need_ex eq 0 then begin
			slpr[ii] = s_g[ii]/nump[ii]
		endif else begin
			slpr[ii] = (s_g[ii]-g[jj])-(2*g[jj]-s_g[ii])
		endelse
	endelse
endfor
print,'slpr:',slpr

openw,lun,path+'slp.dat',/get_lun
writeu,lun,bins
writeu,lun,slpr
writeu,lun,bins
free_lun,lun

end
