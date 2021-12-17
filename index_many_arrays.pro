;+
;
; NAME:
; index_many_arrays
;
; PURPOSE:
; Convenient wrapper. Given index ind, and array variables p0,p1.... will change p0,p1 into indexed arrays.
;   E.g. Simple usage
;   a=[2,4,6,8,10]
;   b=[1,3,5,7,9,11,13,15]
;   index=[1,2,4]
;   index_many_arrays,index,a,b
;   print,a
;   > 4,6,10
;   print,b
;   > 3,5,9
;   Equivalent to a=a[index], b=b[index]
;
;   E.g. More complex use setting index as a string
;   IDL> a=fltarr(30,30)
;   IDL> index='[10:20,5:6]';index is a string in this case
;   IDL> index_many_arrays,index,a
;   IDL> help,a
;   A               FLOAT     = Array[11, 2]
;
; CALLING SEQUENCE:
;   index_many_arrays,ind,p0,p1,p2,p3,p4,p5, $
;                 p6,p7,p8,p9,p10,p11,p12,p13,p14,p15
;
; INPUTS:
;  ind = index e.g. [1,3,4,5,7] or string e.g. '[2:5]'
;
;   p0 -> p15 = arrays on which the indexing will be applied
;
; OUTPUTS:
;   (NONE)
;
; OPTIONAL KEYWORDS
;   (NONE)
;
; OPTIONAL INPUTS
;   (NONE)
;
; OPTIONAL OUTPUTS
;   (NONE)
;
; PROCEDURE:
;
; EXAMPLE:
;
; USE & PERMISSIONS
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

pro index_many_arrays,ind,p0,p1,p2,p3,p4,p5, $
	p6,p7,p8,p9,p10,p11,p12,p13,p14,p15

if size(ind,/type) eq 7 then str=ind else str='[ind]'
for i=0,n_params()-2 do begin
  pname=strcompress('p'+string(i),/remove_all)
  r0=execute('n=n_elements('+pname+')')
  if n gt 0 then r=execute(strcompress(pname+'='+pname+str,/remove_all))
endfor
end