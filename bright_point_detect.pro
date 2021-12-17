;+
;
; NAME:
; bright_point_detect
;
; PURPOSE:
; Identifies and characterises brightenings within a datacube of dimensions [spatial,spatial,temporal] e.g. [x,y,t].
; Brightenings are small regions that brighten temporarily compared to the background.
; The method used is that desribed by https://arxiv.org/abs/2107.13635. Please cite this work if you 
; use the method in your own work.
; The datacube should be given without radiometric calibration, so in suitable units for estimating
; Poisson statistics. User also supplies 'cal_factor', which this program then uses to calibrate each image in 
; the datacube in order to characterise the intensities of each detected brightening.
; The function returns an array of structures, each containing information on each bright point.
; One limitation of this program is that datacube be large enough in all dimensions to enable
; the bandpass filtering applied in the bright_point_filter function.
;
; CALLING SEQUENCE:
; detection_structures=bright_point_detect(datacube_in,cal_factor[,flow=flow,fhigh=fhigh, $
;              time_var_sigma=time_var_sigma, tlow=tlow, thigh=thigh])
;
; INPUTS:
;   datacube_in = A datacube of dimensions [spatial,spatial,temporal] e.g. [Nx,Ny,Nt], in units suitable
;                 for Poisson statistics (i.e. the detector counts e.g. DN)
;                
;   cal_factor = A single variable or vector of calibration factors for each image. After applying filtering
;                 and calculating significance levels using Poisson statistics on datacube, the datacube is then
;                 copied, and calibrated, using cal_factor. If cal_factor is a single variable, then each image
;                 in the datacube is multiplied by cal_factor (i.e. each image has the same calibration factor).
;                 If cal_factor is a vector, then it must have the same number of elements as the number of time
;                 steps in the datacube (=Nt). In this case, each image is multipled by the corresponding element
;                 of cal_factor. cal_factor is often just a factor of the exposure time, thus converting DN to DN/s.
;                 In this case, cal_factor should be set equal to 1/exposure_time.
;
; OUTPUTS:
;   detection_structures = an array of structures, each containing information on each bright point. The structure
;               contains variables:
; BG              FLOAT           3.98933 ;estimated median background intensity in region local to brightening
; TOTB            FLOAT           678.976 ;the total brightness summed over the brightening (minus the background estimate)
; MAXB            FLOAT           20.6044 ;the maximum brightness recorded in the brightening volume
; NVOX            LONG               161  ;number of voxels contained in the brightening
; NT              LONG                 8  ;number of time steps spanned by the brightening
; NFRAG           FLOAT           1.00000 ;over all time steps, the maximum recorded number of isolated 
;                                           regions comprising the brightening, see paper
; MEANSPEED       FLOAT          0.630572 ;the mean speed of the brightening's centroid, in pixel per time step
; MEANXPOS        FLOAT           59.9531 ;the overall mean x-position of the brightening's centroid, in pixels from image left
; MEANYPOS        FLOAT           22.1430 ;the overall mean y-position of the brightening's centroid, in pixels from image bottom
; MEANTPOS        FLOAT           1.97326 ;the overall mean t-position of the brightening's centroid, in time steps from datacube start
; XPOS            POINTER   <PtrHeapVar2> ;the x-position of the brightening's centroid as a function of its lifetime
; YPOS            POINTER   <PtrHeapVar3> ;the y-position of the brightening's centroid as a function of its lifetime
; TPOS            POINTER   <PtrHeapVar4> ;the t-position of the brightening's centroid as a function of its lifetime
; BRT             POINTER   <PtrHeapVar5> ;the total brightness at each time step over it's lifetime 
;                                           (summed over all pixels at a given time step, minus background estimate)
; AREA            POINTER   <PtrHeapVar6> ;the area of the brightening as a function of time over it's lifetime
; FRAG            POINTER   <PtrHeapVar7> ;at each time step, the number of isolated spatial regions comprising the brightening
; A               POINTER   <PtrHeapVar8> ;A, B, and C are approximate estimates of the brightening's shape and orientation calculated
;                                             at each time step, if we assume the brightening to have a 2D Gaussian distribution
;                                             in brightness. See section 'Two-dimensional Gaussian Function' in 
;                                             e.g. https://en.wikipedia.org/wiki/Gaussian_function. We haven't done any proper
;                                             analysis or tests of these, so please take care. They obviously cannot be valid for 
;                                             any detection time step that has frag > 1.
; B               POINTER   <PtrHeapVar9>
; C               POINTER   <PtrHeapVar10>
;
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
; thigh= Sets the threshold for initial detection. Default 9. The variance expected from Poisson noise (after bandpass filtering)
;       is multiplied by tlow, and any regions above this threshold are initially considered as brightening candidate events.
; 
; tlow= Sets the threshold for brightening characterisation. Default 7. The initial brightenings are identified using Thigh. 
;       We then loop through each brightening and set a new threshold which is less stringent, thus increasing the 
;       brightening volume. See paper. 
; 
; minvol= The minimum number of voxels for a detection to be considered. Default 25.
; 
; mintstep= The minimum number of time steps for a detection to be considered. Default 5.
; 
; local_range= The number of pixels surrounding each initial detection, defines a local region over which
;               the background intensity is estimated, plus over which the tlow threshold is applied.
;
; OPTIONAL INPUTS
;   (NONE)
;
; OPTIONAL OUTPUTS
;   (None)
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

function bright_point_detect, datacube_in,cal_factor,flow=flow,fhigh=fhigh, $
              time_var_sigma=time_var_sigma, tlow=tlow, thigh=thigh, $
              minvol=minvol,mintstep=mintstep,local_range=local_range

  thresh=n_elements(thigh) eq 0?9:thigh
  thresh2=n_elements(tlow) eq 0?7:tlow
  minvol=keyword_set(minvol)?minvol:25
  mintstep=keyword_set(mintstep)?mintstep:5
  local_range=keyword_set(local_range)?local_range:7

  sizearr, datacube_in, nx, ny, nt,ndim=ndim
  if ndim ne 3 or nx lt 60 or ny lt 60 or nt lt 60 then begin
    print,'Please provide 3-dimensional datacube of sufficient size to apply filtering and detection'
    return,-1
  endif
  
  ncal=n_elements(cal_factor)
  case 1 of
    ncal eq 0:begin
      print,'Please provide calibration factors for images'
      return,-1
    end
    ncal eq 1:begin
      ;if user provides just one exposure time, assume all images have same exposure time
      cal_factor2=fltarr(nt)+cal_factor
    end
    ncal gt 1 and ncal ne nt:begin
      print,'Please provide vector of calibration factors of same number of time steps as datacube'
      return,-1 
    end
    else:cal_factor2=cal_factor
  endcase
  
  ;correct for exposure time (this is done after filter since the sigma thresholds are based on Poisson statistics)
  datacube=fltarr(nx,ny,nt)
  for i=0, nt-1 do datacube[*,*,i]=datacube_in[*,*,i]*cal_factor2[i]
  
  datacube_filt=bright_point_filter(datacube_in,sigma=sigma,flow=flow,fhigh=fhigh, $
    time_var_sigma=time_var_sigma);note we filter datacube_in, not corrected for exposure time (Poisson stats)
  
  if n_elements(datacube_filt) eq 1 then begin
    print,'Problem with datacube, returning (detect_iris_datacube)'
    return,-1
  endif
  
  ;calculate ratio of filtered datacube to sigma 
  indnan=where(~finite(datacube_filt),cntnan)
  datacube_filt=datacube_filt-median(datacube_filt)
  sizearr,sigma,nxsig,nysig,ntsig,ndim=ndimsig
  if ~(ndimsig eq 3) then sigma=rebin(sigma,nx,ny,nt)
  ratio=(datacube_filt>0)/sigma
  if cntnan gt 0 then ratio[indnan]=!values.f_nan

  ;mask of significantly bright regions
  mask=ratio gt thresh
  mask=pad_3d(mask,1)

  ;filter out small-volume or very large-volume elements
  lr=label_region(mask,/ulong)
  h=histogram(lr,min=1,max=max(lr),loc=xh,rev=ri)
  indvol=where(h lt minvol or h gt robust_max(h,1),cntvol)
  for i=0,cntvol-1 do begin
    ind=get_rev_ind(ri,indvol[i],n)
    mask[ind]=0
  endfor
  lr=label_region(mask,/ulong)

  ;filter out very short or very long regions (in time)
  h=histogram(lr,min=1,max=max(lr),loc=xh,rev=ri)
  nh=n_elements(h)
  for i=0,nh-1 do begin
    ind=get_rev_ind(ri,i,n)
    one2n,ind,mask,ix,iy,it
    if max(it)-min(it) lt mintstep then begin
      mask[ind]=0
    endif ;else print, max(it)-min(it)
  endfor

  lr=label_region(mask,/ulong)
  mask=unpad_3d(mask,1)
  lr=unpad_3d(lr,1)

  ;now begin BP analysis
  h=histogram(lr,min=1,max=max(lr),loc=xh,rev=ri)
  nh=n_elements(h)

  masksmo=smooth(float(mask),5,/edge_trunc) gt 0.01
  maskmaster=lonarr(nx,ny,nt)
  brtmaster=fltarr(nx,ny,nt)

  

  d={ $
    bg:0., $
    totb:0., $
    maxb:0., $
    nvox:0l, $
    nt:0l, $
    nfrag:0., $
    meanspeed:0., $
    meanxpos:-1., $
    meanypos:-1., $
    meantpos:-1., $
    xpos:ptr_new(), $
    ypos:ptr_new(), $
    tpos:ptr_new(), $
    brt:ptr_new(), $
    area:ptr_new(), $
    frag:ptr_new(), $
    a:ptr_new(), $
    b:ptr_new(), $
    c:ptr_new() $
  }
  d=replicate(d,nh)

  iprint=long(nh/20.)
  for i=0l,nh-1 do begin

    if i mod iprint eq 0 then print,i,' out of ',nh-1

    ind=get_rev_ind(ri,i,n)
    one2n,ind,mask,ix,iy,it

    ix0=(min(ix)-local_range)>0
    ix1=(max(ix)+local_range)<(nx-1)
    iy0=(min(iy)-local_range)>0
    iy1=(max(iy)+local_range)<(ny-1)
    it0=(min(it)-local_range)>0
    it1=(max(it)+local_range)<(nt-1)

    ix=ix-ix0
    iy=iy-iy0
    it=it-it0
    imnow=datacube[ix0:ix1,iy0:iy1,it0:it1];original image values
    maskbp=lr[ix0:ix1,iy0:iy1,it0:it1] eq xh[i];only this current bright point (high threshold)
    rationow=ratio[ix0:ix1,iy0:iy1,it0:it1];ratio
    sizearr,imnow,nxnow,nynow,ntnow

    maskbp2=rationow gt thresh2;mask for low threshold
    l2=unpad_3d(label_region(pad_3d(maskbp2,1)),1)
    indbp=where(maskbp)
    maskbp2=l2 eq l2[indbp[0]];only keep the region including original BP detection

    ;BG estimate
    indnobp=where(~maskbp2,cntnobp)
    bgestimate=median(imnow[indnobp])
    indbp2=where(maskbp2)

    brt=(((imnow-bgestimate)>0)*maskbp);brightness in core of BP (high threshold)
    brt2=(((imnow-bgestimate)>0)*maskbp2);brightness in extended BP (low threshold)

    tot=total(total(brt2,1),1)
    indbp=where(tot gt 0,nbp)

    if nbp eq 0 then continue

    l=lonarr(nxnow,nynow,nbp)
    for it=0,nbp-1 do begin
      itim=it+min(indbp)
      mnow=maskbp2[*,*,itim]
      lnow=unpad_2d(label_region(pad_2d(mnow,1)),1)
      l[*,*,it]=lnow
    endfor
    nb=max(max(l,dim=1),dim=1)
    nbmx=max(nb)

    maskmaster[ix0:ix1,iy0:iy1,it0:it1]=maskmaster[ix0:ix1,iy0:iy1,it0:it1]>(maskbp2*nbmx)
    brtmaster[ix0:ix1,iy0:iy1,it0:it1]=brtmaster[ix0:ix1,iy0:iy1,it0:it1]>brt2

    mxt=fltarr(nbp)
    mxt[*]=!values.f_nan
    a=fltarr(nbp)
    a[*]=!values.f_nan
    b=fltarr(nbp)
    b[*]=!values.f_nan
    c=fltarr(nbp)
    c[*]=!values.f_nan
    ixmx=fltarr(nbp)
    ixmx[*]=!values.f_nan
    iymx=fltarr(nbp)
    iymx[*]=!values.f_nan
    iix=rebin(findgen(nxnow),nxnow,nynow)
    iiy=rebin(reform(findgen(nynow),1,nynow),nxnow,nynow)
    itt=lonarr(nbp)
    totbrtmain=fltarr(nbp)
    area=fltarr(nbp)
    for it=0,nbp-1 do begin
      itim=it+min(indbp)

      brtnow=brt[*,*,itim]
      brtnow2=brt2[*,*,itim]

      totbrt=total(brtnow)
      ixnow=total(brtnow*iix)/totbrt
      iynow=total(brtnow*iiy)/totbrt
      mxt[it]=max(brtnow,/nan);,ind)
      ;one2n,ind,[nxnow,nynow],ixnow,iynow,/dim

      ixmx[it]=ixnow+ix0
      iymx[it]=iynow+iy0

      a[it]=total(brtnow*((iix-ixnow)^2))/totbrt
      b[it]=total(brtnow*((iiy-iynow)^2))/totbrt
      c[it]=total(brtnow*((iix-ixnow)*(iiy-iynow)))/totbrt

      itt[it]=itim+it0

      totbrtmain[it]=total(brtnow2)
      area[it]=total(maskbp2[*,*,itim])


    endfor;it
    
    xmean=total(rebin(iix,nxnow,nynow,nbp)*brt[*,*,indbp])/total(brt[*,*,indbp])+ix0
    ymean=total(rebin(iiy,nxnow,nynow,nbp)*brt[*,*,indbp])/total(brt[*,*,indbp])+iy0
    tmean=total(rebin(reform(indbp,1,1,nbp),nxnow,nynow,nbp)*brt[*,*,indbp])/total(brt[*,*,indbp])+it0
    
    ind=where(finite(ixmx),cnt)
    index_many_arrays,ind,ixmx,iymx,a,b,c,mxt,itt,totbrtmain,area,nb
    ind=sort(itt)
    index_many_arrays,ind,ixmx,iymx,a,b,c,mxt,itt,totbrtmain,area,nb


    d[i].nvox=n ;number of voxels above thresh_low??
    d[i].nt=nbp ;number of frames
    d[i].bg=bgestimate
    d[i].totb=total(brt2)
    d[i].maxb=max(brt)
    d[i].nfrag=max(float(nb))
    d[i].meanxpos=xmean
    d[i].meanypos=ymean
    d[i].meantpos=tmean
    dx=finite_difference_1d(ixmx);pixels
    dy=finite_difference_1d(iymx)
    dt=1;time step
    d[i].meanspeed=mean(sqrt(dx^2+dy^2)/dt)
    d[i].xpos=ptr_new(ixmx)
    d[i].ypos=ptr_new(iymx)
    d[i].tpos=ptr_new(itt)
    d[i].brt=ptr_new(totbrtmain)
    d[i].area=ptr_new(area)
    d[i].frag=ptr_new(nb)
    d[i].a=ptr_new(a)
    d[i].b=ptr_new(b)
    d[i].c=ptr_new(c)

  endfor

  index=where(d.nvox gt 0,cnt)
  if cnt eq 0 then begin
    print,'No detections'
    return,-1
  endif
  d=d[index]
  
  return,d

end
