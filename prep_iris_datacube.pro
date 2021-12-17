;+
;
; NAME:
; prep_iris_datacube
;
; PURPOSE:
; Given an user-supplied name and path to an IRIS slitjaw Level2 fits file this procedure opens the 
; datacube and header, automatically masks off regions of missing data as NAN, and uses an image alignment procedure
; to correct for any spatial drift of the field of view over time. The user can crop the
; datacube to a spatial subregion, or crop in time, via the optional subregion and timerange keywords. 
;
; CALLING SEQUENCE:
; prep_iris_datacube,fitsfile [,datacube,hdr,savedir=savedir, lim=lim,  $
;                     subregion=subregion,create_imagealign_plots=create_imagealign_plots, $
;                     fit_shift=fit_shift,timerange=timerange]
;
; INPUTS:
;   fitsfile = full path and filename of an IRIS Level2 fits file containing the slitjaw datacube and header 
;
; OUTPUTS:
;   datacube = the processed datacube
;   
;   hdr = the header metadata 
;
; OPTIONAL KEYWORDS
;
;
; savedir = specifies directory where user requires output plots to be saved. If not set
; then the output plot is saved in the directory of the fitsfile.
;
; timerange =  range of time step indices units
;   e.g. timerange=[50,220] will extract time steps from 50 to 220 inclusive
;   so datacube [x,y,t] will be shortened to datacube[*,*,50:220]
;   Default is to use all time steps in datacube
;
; lim = multiplier. If there are more than lim*nx*nt bad pixels in
;   a single image column or row (x,t plane or y,t plane)
;   then that column is defined as bad throughout the whole datacube and set to NAN.
;   Useful for automatically setting margins and central slit position as missing data.
;   Default value is 0.25.
;
; Subregion = 4-element vector defining a sub region. The datacube will be cropped
;       to this subregion. Format is [xleft, ybottom, xright, ytop] in pixels.
;
; fit_shift = keyword controlling how the image translations over time are calculated, for
;   datacube alignment. If set (/fit_shift), then the per-time-step cumulative global image shifts
;   in x and y are fitted to a straight line. If not set, the default is to smooth the cumulative shifts over time.
;   The alignment corrects for region drift due to e.g. solar rotation, using FLCT
;
; OPTIONAL INPUTS
;   (NONE)
;
; OPTIONAL OUTPUTS
;   datacube = the output datacube, aligned and optionally cropped
;   
;   hdr = the modified header
;
; PROCEDURE:
; See https://arxiv.org/abs/2107.13635
;
; USE & PERMISSIONS
; If you use this code in your own work, please cite https://arxiv.org/abs/2107.13635
; Any problems/queries, or suggestions for improvements, please email Huw Morgan, hmorgan@aber.ac.uk
;
; ACKNOWLEDGMENTS:
;  This code was developed with the financial support of:
;  STFC and Coleg Cymraeg Cenedlaethol Studentship to Aberystwyth University (Humphries)
;  STFC Consolidated grant to Aberystwyth University (Morgan)
;
; MODIFICATION HISTORY:
; Created at Aberystwyth University 2021 - Huw Morgan hmorgan@aber.ac.uk
;
;
;-


pro prep_iris_datacube,fitsfile,datacube,hdr,tai,savedir=savedir, lim=lim,  $
        subregion=subregion,create_imagealign_plots=create_imagealign_plots, $
        fit_shift=fit_shift,timerange=timerange

  if ~keyword_set(savedir) then savedir=file_dirname(fitsfile);savedir is used to save optional plots

  ;Read in IRIS Level 2 SJI file
  read_iris_l2,fitsfile,hdr,datacube
  datacube=float(datacube>0);convert to float, set negative values to zero
  
  if keyword_set(timerange) then begin
    print,'Restricting time range of datacube to ',timerange
    datacube=datacube[*,*,timerange[0]:timerange[1]]
    hdr=hdr[timerange[0]:timerange[1]]
  endif
  
  ;query size of datacube
  sizearr,datacube,nx,ny,nt

  ;ad-hoc removal of slit
  ;AFTER new FOV: datacube[318:321,*,0]=!values.f_nan
  ;datacube[360:364,*,0]=!values.f_nan ;BEFORE new FOV, 2014-
  sltpx=hdr[0].sltpx1ix
  datacube[sltpx-1:sltpx+1,*,*]=!values.f_nan;2014-02-20

  if keyword_set(subregion) then begin

    if n_elements(subregion) ne 4 then begin
      print,'Subregion keyword should have 4 elements giving [xleft, ybottom, xright, ytop] in pixels'
      print,'Datacube not cropped in space'
    endif else begin

      xra=subregion[[0,2]]
      yra=subregion[[1,3]]
      xra[0]=xra[0]>0
      xra[1]=xra[1]<(nx-1)
      yra[0]=yra[0]>0
      yra[1]=yra[1]<(ny-1)

      datacube=datacube[xra[0]:xra[1],yra[0]:yra[1]]

      hdr.naxis1=xra[1]-xra[0]+1
      hdr.naxis2=yra[1]-yra[0]+1
      hdr.fovx=hdr.naxis1*hdr.cdelt1
      hdr.fovy=hdr.naxis2*hdr.cdelt2
      hdr.crpix1=hdr.crpix1-xra[0]
      hdr.crpix2=hdr.crpix2-yra[0]

    endelse
  endif

  ;set 'bad' regions of image to NANs. Subsequently treated as missing data.
  maskbad=datacube eq 0 or datacube gt 150*10.
  if n_elements(lim) eq 0 then lim=0.25
  tx=total(total(long(datacube lt 1),3),2)
  indx=where(tx gt ny*nt*lim) ;0.25 default
  ;stop
  maskbad[indx,*,*]=1
  ty=total(total(long(datacube lt 1),3),1)
  indy=where(ty gt nx*nt*lim) ;0.25 default
  maskbad[*,indy,*]=1
  maskbad=float(maskbad)
  maskbad=smooth(maskbad, 3, /edge_trunc)
  maskbad=maskbad gt 0.01
  indnan=where(maskbad,cntnan)
  if cntnan gt 0 then datacube[indnan]=!values.f_nan

  ;do sub-pixel alignment on datacube time series
  sizearr,datacube,nx,ny,nt
  
  ;Use fft_corr_track_huw to calculate image-by-image shifts
  xy=fltarr(2,nt)
  print,''
  print,'Calculating image shifts at each time step...'
  for i=1,nt-1 do begin
    if i mod 40 eq 0 then print,i,' out of ',nt-1
    xy[*,i]=fft_corr_track_huw(datacube[*,*,i-1],datacube[*,*,i])
  endfor
  
  ;xs and ys are the raw x and y cumulative shifts
  xs=total(reform(xy[0,*]),/cum)
  ys=total(reform(xy[1,*]),/cum)
  
  ;
  tai=anytim2tai(hdr.startobs)+dindgen(nt)*hdr.cdelt3
  midtai=mean(tai)
  tref=tai-midtai
  if keyword_set(fit_shift) then begin
    ;fit the cumulative shifts to a linear function of time
    px=ladfit(tref,xs)
    py=ladfit(tref,ys)
    xy=[[px[0]+px[1]*tref],[py[0]+py[1]*tref]]
    xy=rotate(xy,4)
  endif else begin
    ;smooth the cumulative shifts
    xy=[[smooth(median(xs,15),5,/edge_trunc)], $
      [smooth(median(ys,15),5,/edge_trunc)]]
    xy=rotate(xy,4)
  endelse

  ;highly recommend taking a look and checking image alignment
  ;when using for first time on new dataset.
  if keyword_set(create_imagealign_plots) then begin
    psfile=savedir+'/iris_align.ps'
    set_ps_huw,psfile,7,8
    !p.multi=[0,1,2]
    plot,tref,xs,/psym,title='X shift',xtitle='Time (sec from median)',ytitle='Shift (pixel)'
    oplot,tref,reform(xy[0,*]),col=cgcolor('red')
    plot,tref,ys,/psym,title='Y shift',xtitle='Time (sec from median)',ytitle='Shift (pixel)'
    oplot,tref,reform(xy[1,*]),col=cgcolor('red')
    setx,/cle,/clo
    loadct,0,/sile
    spawn,'open '+psfile
  endif
  
  print,'Applying image alignment...'
  datacube2=datacube
  hdr2=hdr
  for i=1,nt-1 do begin
    if i mod 40 eq 0 then print,i,' out of ',nt-1
    datacube2[*,*,i]=image_align(datacube[*,*,i-1],datacube[*,*,i],xy=xy[*,i],/silent,/use_given,cubic=0)
    hdr2[i].crpix1=hdr2[i-1].crpix1
    hdr2[i].crpix2=hdr2[i-1].crpix2
    hdr2[i].xcen=hdr2[i-1].xcen
    hdr2[i].ycen=hdr2[i-1].ycen
    hdr2[i].crval1=hdr2[i-1].crval1
    hdr2[i].crval2=hdr2[i-1].crval2
  endfor
  datacube=datacube2
  hdr=hdr2

end
