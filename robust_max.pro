;+
;
; NAME:
; robust_max
;
; PURPOSE:
; Calculates the robust maximum of an input array. I believe it uses linear interpolation between closest ranks method.
; Wrote intiuitively, but really need to visit the literature and write properly
; e.g. https://en.wikipedia.org/wiki/Percentile#Weighted_percentile
;
; CALLING SEQUENCE:
; min=robust_max(array[,percent,min=min,dimension=dimension])
;
; INPUTS:
;   array=input array, numerical
;
; OPTIONAL INPUT:
;   percent=the percentile minimum to use, default is 1%
;
; OUTPUTS:
;   function returns robust minimum
;
;
; OPTIONAL KEYWORD
;   min = returns robust minimum of the array at same percentile
;
;   dimension = calculate the robust max over this dimension of the array, see dimension keyword for IDL 
;               functions such as mean, median, min, max etc. Note this feature is currently only implemented
;               for arrays of 3 dimensions or less.
;   
; PROCEDURE:
;   Sorts array members into ascending order, then interpolates to the value
;   where number of array members is equal to the percentile
;
; USE & PERMISSIONS
;  If you reuse in your own code, please include acknowledgment to Huw Morgan (see below)
;
; ACKNOWLEDGMENTS:
;  This code was developed with the financial support of:
;  STFC Consolidated grant to Aberystwyth University (Morgan)
;
; MODIFICATION HISTORY:
; Created at Aberystwyth University 07/2019 - Huw Morgan hmorgan@aber.ac.uk
;
;
;-

function robust_max,y0,per,min=min,dimension=dimension

if n_params() lt 2 then per=1

if ~keyword_set(dimension) then begin
  indok=where(finite(y0),n)
  if n eq 0 then return, !values.f_nan
  indsort=sort(y0[indok])
  indmin=per*float(n-1)/100.
  indmax=(100-per)*float(n-1)/100.
  min=interpol(y0[indok[indsort]],findgen(n),indmin)
  max=interpol(y0[indok[indsort]],findgen(n),indmax)

endif else begin

  sizearr,y0,nx,ny,nz,ndim=ndim
  ind=sort_nd(y0,dimension)
  case dimension of
    1:begin
      indmin=per*float(nx-1)/100.
      indmax=(100-per)*float(nx-1)/100.
      if arg_present(min) then begin
        min=interpolate(y0[ind],indmin,findgen(ny),findgen(nz),/grid)
        min=reform(min)
      endif
      case ndim of
        2:begin
          max=interpolate(y0[ind],indmax,findgen(ny),/grid)
          max=reform(max)
        end
        3:begin
          max=interpolate(y0[ind],indmax,findgen(ny),findgen(nz),/grid)
          max=reform(max)
        end
      endcase
    end
    2:begin
      indmin=per*float(ny-1)/100.
      indmax=(100-per)*float(ny-1)/100.
      if arg_present(min) then begin
        min=interpolate(y0[ind],findgen(nx),indmin,findgen(nz),/grid)
        min=reform(min)
      endif
      max=interpolate(y0[ind],findgen(nx),indmax,findgen(nz),/grid)
      max=reform(max)
      case ndim of
        2:begin
          max=interpolate(y0[ind],indmax,findgen(ny),/grid)
          max=reform(max)
        end
        3:begin
          max=interpolate(y0[ind],indmax,findgen(ny),findgen(nz),/grid)
          max=reform(max)
        end
      endcase
    end
    3:begin
      indmin=per*float(nz-1)/100.
      indmax=(100-per)*float(nz-1)/100.
      if arg_present(min) then min=interpolate(y0[ind],findgen(nx),findgen(ny),indmin,/grid)
      max=interpolate(y0[ind],findgen(nx),findgen(ny),indmax,/grid)
    end
  endcase

endelse

return,max

end