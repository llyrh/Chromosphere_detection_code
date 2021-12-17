;+
;
; NAME:
; bright_point_filter
;
; PURPOSE:
; Applies bandpass filtering across all dimensions of a datacube of dimensions [spatial,spatial,temporal] e.g. [x,y,t].
; Calculates the variance expected from Poisson noise for a filtered datacube.
; The method used is that desribed by https://arxiv.org/abs/2107.13635. Please cite this work if you
; use the method in your own work.
; The datacube should be given without radiometric calibration, so in suitable units for estimating
; Poisson statistics. 
; One limitation of this program is that datacube be large enough in all dimensions to enable
; the bandpass filtering.
;
; CALLING SEQUENCE:
; datacube_out=bright_point_filter(datacube_in [,flow=flow,fhigh=fhigh, $
;                            sigma=sigma,time_var_sigma=time_var_sigma])
;
; INPUTS:
;   datacube_in = A datacube of dimensions [spatial,spatial,temporal] e.g. [Nx,Ny,Nt], in units suitable
;                 for Poisson statistics (i.e. the detector counts e.g. DN)
;
; OUTPUTS:
;  datacube_out = The filtered datacube
;
; OPTIONAL KEYWORDS
;
; flow= sets the lower frequency for bandpass filtering, see IDL's digital_filter. Default 0.09
;
; fhigh= sets the lower frequency for bandpass filtering, see IDL's digital_filter. Default 0.2
;
; time_var_sigma= keyword, if set, is to create a time-varying estimate of Poisson significance. Default is to
;                   estimate the significance over the whole datacube time range
;
; OPTIONAL INPUTS
;   (NONE)
;
; OPTIONAL OUTPUTS
; 
; sigma = 2D or 3D array giving the expected variance from a datacube comprised purely from bandpass-filtered Poisson noise.
;         Default is a 2D array. If user sets keyword /time_var_sigma, then sigma is a 3D array (time-varying variance).
;
; PROCEDURE:
; See https://arxiv.org/abs/2107.13635. Since publishing this paper, the capability for time-varying variance has been added 
; to this routine.
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

function bright_point_filter,datacube_in,flow=flow,fhigh=fhigh, $
                            sigma=sigma,time_var_sigma=time_var_sigma


  if n_elements(flow) eq 0 then flow=0.09      ;low frequency limit
  if n_elements(fhigh) eq 0 then fhigh=0.2       ;high frequency limit
  xy_nterms=10;number of spatial kernel terms
  t_nterms=10  ;number of temporal kernel terms
  filter_power=50        ;
  
  ;query size of datacube
  sizearr, datacube_in, nx, ny, nt
  if n_elements(nt) eq 0 then begin
    print,'Please provide 3-dimensional datacube'
    return,-1
  endif
  
  indnan=where(~finite(datacube_in),cntnan)
  
  print,'Creating filter coefficients'
  coeffx=digital_filter(flow, fhigh, filter_power, xy_nterms)
  coeffy=digital_filter(flow, fhigh, filter_power, xy_nterms)
  coefft=digital_filter(flow, fhigh, filter_power,t_nterms)

  coeffy=reform(coeffy,1,n_elements(coeffy),1)
  coefft=reform(coefft,1,1,n_elements(coefft))
  
  ;check that datacube is larger than the filter kernels
  ncx=n_elements(coeffx)
  ncy=n_elements(coeffy)
  nct=n_elements(coefft)
  if nx lt ncx or ny lt ncy or nt lt nct then begin
    print,'Please provide datacube of size greater than the filter coefficients'
    return,-1
  endif

  ;Estimating thresholds for significant bright regions
  ;based on amplitude of Poisson noise
  ;Estimate standard deviation of filtered Poisson noise (with means based on data)
  ;Then apply filter to determine the noise amplitude in filtered datacube
  
  print,'Calculating filtered noise standard deviations'
  ;numerical estimate of standard deviation of filtered noise
  nlevels=5;estimate sigma for 5 different signal levels ranging from min to max of data
  nmod=ncx*3;width of small datacube used to create model noise
  mnlev=interpol(minmax(datacube_in,/nan),nlevels);mean signal levels to create noise models for
  sigmamod=fltarr(nlevels);to store measurements of sigma at each mean signal level
  for ilevel=0,nlevels-1 do begin;loop through mean signal levels
    noisemodel=randomn(seed,[nmod,nmod,nmod],poisson=mnlev[ilevel]);create Poisson noise model based on current mean value
    fsd=convol(noisemodel,  coeffx,/edge_trunc,/nan);and filter as we will apply to data
    fsd=convol(fsd,coeffy,/edge_trunc,/nan)
    fsd=convol(fsd,coefft,/edge_trunc,/nan)
    sigmamod[ilevel]=stddev(fsd);store standard deviation of filtered noise model
  endfor
  factor=mean(sigmamod/sqrt(mnlev));sigmamod always in proportion to square root of input signal
  
  ;calculate mean of data across time - either across all time, or locally
  if ~keyword_set(time_var_sigma) then begin
     mn=mean(datacube_in,dim=3,/nan)
  endif else begin
    ;mean calculated locally in time
    sigtwdth=7
    kert=gaussian_function(sigtwdth,/norm)
    kert=reform(kert,1,1,n_elements(kert))
    mn=convol(datacube_in,kert,/edge_trunc,/nan);time-local mean
  endelse
  sigma=sqrt(mn)*factor

  print,'Filtering datacube...'
  ;Filter input datacube
  datacube_out=convol(datacube_in,coeffx,/edge_trunc,/nan)
  datacube_out=convol(datacube_out,coeffy,/edge_trunc,/nan)
  datacube_out=convol(datacube_out,coefft,/edge_trunc,/nan)
  
  if cntnan gt 0 then datacube_out[indnan]=!values.f_nan
  
  return,datacube_out

end