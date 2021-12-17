

pro call_iris_master

;specify top level directory where IRIS data is kept
topdatadir='/Users/hum2/data/llyr_bright_points'

if n_elements(date) eq 0 then begin;just for convenience when testing.
  overwrite=0
  date='2016/11/04';User can change this date to their existing datasets for testing
  ;subregion=[0,-320,70,-250];and change this
  create_imagealign_plots=1
  fit_shift=0
  timerange=[130,454]
  time_var_sigma=1
endif

;data stored in this top level directory under subdirectories labelled yyyymmdd (e.g. 20161105 for 5th November 2016)
date8=anytim2cal(date,form=8,/date)
datadir=topdatadir+'/'+date8
if ~keyword_set(savedir) then savedir=datadir;if user has not specified differently, save output in data directory

;search for slitjaw level 2 fits files
searchstring=datadir+'/iris_l2_'+date8+'_*_SJI_*.fits'
fitsfiles=file_search(searchstring,count=nfitsfiles)

if nfitsfiles eq 0 then begin
  print,'DETECTION_IRIS_MASTER: No IRIS level 2 FITS files found'
  print,'Search string used is ',searchstring
  print,'Returning'
  return
endif

;identify wavelength of these data from the fits filename
wl=lonarr(nfitsfiles)
for ifile=0,nfitsfiles-1 do begin
  f=file_basename(fitsfiles[ifile],'.fits')
  f=strsplit(f,'_',/extract)
  ind=where(strmatch(f,'SJI'))+1
  wl[ifile]=long(f[ind])
endfor

;if user has not set wlreq (wavelength required) keyword, then loop through all existing wavelengths
if ~keyword_set(wlreq) then wlreq=wl

nwlreq=n_elements(wlreq)
for iwlreq=0,nwlreq-1 do begin;loop through files

  indwl=where(wl eq wlreq[iwlreq],cntwl);does file exist for this WL?
  if cntwl eq 0 then begin
    print,'No FITS file for wavelength ',wlreq[iwlreq]
    continue
  endif

  fitsfilenow=fitsfiles[indwl]
  savename=datadir+'/'+file_basename(fitsfilenow,'.fits')+'_bp.dat'
  
  prep_iris_datacube,fitsfilenow,datacube,hdr,savedir=savedir, lim=lim,  $
    subregion=subregion,create_imagealign_plots=create_imagealign_plots, $
    fit_shift=fit_shift,timerange=timerange

  det=bright_point_detect(datacube,1/hdr.exptime,flow=flow,fhigh=fhigh, $
              time_var_sigma=time_var_sigma, tlow=tlow, thigh=thigh)
              
  if ~is_struct(det) then begin
    print,'No success with detections for ',fitsfilenow
    continue
  endif
  
  
  save,datacube,hdr,det,filename=savename

endfor

end
