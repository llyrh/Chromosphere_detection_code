;+
;
; NAME:
; image_align
;
; PURPOSE:
; Given x,y pixel translation of input image im1 relative to reference input image im0, this function applies
; interpolation in order to translate im1 so it is optimally aligned with im0.
; Method described well by
; Fisher & Welsh http://articles.adsabs.harvard.edu/pdf/2008ASPC..383..373F (and references within)
;
; CALLING SEQUENCE:
; imageout=image_align(im0,im1 [,xy=xy,silent=silent, $
;        whole_pixel=whole_pixel,log=log,subregion=subregion, $
;        percent=percent,use_given_xy=use_given_xy,cubic=cubic])
;
; INPUTS:
;   im0 = 2-dimensional image
;
;   im1 = 2-dimensional image. Can be different size to im0
;
; OUTPUTS:
;   imageout = output image, which is im1 shifted to align with im0
;  
;
; OPTIONAL KEYWORDS
;
; xy : returns the x,y pixel shift, or if use_given_xy keyword is supplied by user,
;     uses this supplied xy value for the shift (does not calculate a new shift)
;
; use_given_xy  : if set, uses the user-supplied xy to translate IM1, rather than calculate the translation
;    
; whole_pixel : Passed to fft_corr_track_huw. Calculate whole pixel translations only (default sub-pixel)
;
; subregion : Passed to fft_corr_track_huw. User can define a subregion of the image in order 
;             to calculate the shift (does not crop the return image)
;             Subregion=[xleft,ybottom,xright,ytop] in pixels
;
; percent  : Passed to fft_corr_track_huw. Clips image values to a given percentile
;           e.g. percent=[4,98] limits values to 4th percentile minium and 2nd percentile maximum
;
; log : takes logarithm base 10 of image to calculate the shift (but does not change return image values)
; 
; silent : suppresses messages
;
; OPTIONAL INPUTS
;   (NONE)
;
; OPTIONAL OUTPUTS
;   (None)
;
; PROCEDURE:
;   See Fisher & Welsh http://articles.adsabs.harvard.edu/pdf/2008ASPC..383..373F
;
; USE & PERMISSIONS
; If you use this code in your own work, please cite
; Fisher & Welsh http://articles.adsabs.harvard.edu/pdf/2008ASPC..383..373F (and references within)
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


function image_align,im0,im1,xy=xy,silent=silent, $
        whole_pixel=whole_pixel,log=log,subregion=subregion, $
        percent=percent,use_given_xy=use_given_xy,cubic=cubic

if ~keyword_set(use_given_xy) or n_elements(xy) eq 0 then begin
  if keyword_set(log) then $
  xy=fft_corr_track_huw(alog10(im0+mean(im0,/nan)),alog10(im1+mean(im1,/nan)), $
                          whole_pixel=whole_pixel,subregion=subregion,percent=percent) $
  else $ 
  xy=fft_corr_track_huw(im0,im1,whole_pixel=whole_pixel,subregion=subregion,percent=percent)

  if ~keyword_set(silent) then print,xy
endif
 
sizearr,im0,nx0,ny0
x0=findgen(nx0)
y0=findgen(ny0)

sizearr,im1,nx1,ny1
x1=findgen(nx1)
y1=findgen(ny1)

ix=interpol(x1,x1-xy[0],x0)
iy=interpol(y1,y1-xy[1],y0)
im1out=interpolate(im1,ix,iy,missing=!values.f_nan,/grid,cubic=cubic)

return,im1out

end