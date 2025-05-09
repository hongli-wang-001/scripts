;*************************************************
; polyg_8_lbar.ncl
;
; Concepts illustrated:
;   - Drawing a scatter plot on a map
;   - Changing the marker color and size in a map plot
;   - Plotting station locations using markers
;   - Manually creating a labelbar
;   - Adding text to a plot
;   - Generating dummy data using "random_uniform"
;   - Binning data
;*************************************************
;
; This example creates some dummy station data,
; and then plots each value by coloring and sizing it
; depending on which range of values it falls in.
;
;*************************************************
;
; These files are loaded by default in NCL V6.2.0 and newer
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

;----------------------------------------------------------------------
; Procedure for adding a labelbar at a given NDC location, given
; the levels and colors to use.
;----------------------------------------------------------------------
undef("add_labelbar")
procedure add_labelbar(wks,plot,levels,colors)
local lbres, labels
begin
  nlevels = dimsizes(levels)
;----------------------------------------------------------------------
; Draw a labelbar
;----------------------------------------------------------------------
  lbres                    = True
  lbres@vpWidthF           = 0.80             ; width
  lbres@vpHeightF          = 0.10             ; height
  lbres@lbPerimOn          = False            ; Turn off perimeter.
  lbres@lbOrientation      = "Horizontal"     ; Default is vertical.
  lbres@lbLabelAlignment   = "InteriorEdges"  ; Default is "BoxCenters"
  lbres@lbFillColors       = colors           ; Colors for boxes.
  lbres@lbMonoFillPattern  = True             ; Fill them all solid.
  lbres@lbLabelFontHeightF = 0.012            ; label font height

  labels = sprintf("%4.2f",levels)

  lbid = gsn_create_labelbar(wks,nlevels+1,labels,lbres)
   
;
; Now, create some annotation resources indicating how we want to
; attach the labelbar to the plot. Here, we are using the top right
; corner of the labelbar as the point which we are going to position
; it, and then we use amParallelPosF and amOrthogonalPosF to indicate
; where we want to place it.
;
; amParallelPosF/amOrthogonalPosF                                       ;
;    0.0/ 0.0  -  annotation in dead center of plot
;    0.5/ 0.5  -  annotation at bottom right of plot
;    0.5/-0.5  -  annotation at top right of plot
;   -0.5/-0.5  -  annotation at top left of plot
;   -0.5/ 0.5  -  annotation at bottom left of plot
;                                                                 
  amres                  = True
  amres@amJust           = "TopCenter"
  amres@amParallelPosF   =  0.0   ; keep labelbar centered
  amres@amOrthogonalPosF =  0.6   ; move down and outside of plot

;
; Give both annotation id and labelbar id unique names.
;
; Attaching them to plot with unique names ensures that
; labelbar "lives" outside this procedure.
;
  tmpid1 = "anno"+unique_string("id")
  tmpid2 = "lbar"+unique_string("id")

  plot@$tmpid1$ = gsn_add_annotation(plot,lbid,amres)
  plot@$tmpid2$ = lbid
end

;----------------------------------------------------------------------
; Main code.
;----------------------------------------------------------------------
begin
;-------Options--------

  ;levels = fspan(1,41,21) ; (/0.,5.,10.,15.,20.,23.,26./) ; bin settings (bin0 = < 0., 
  levels = (/0.1,1,2,4,6,8,10,12,14,20,25,30,40,60,100,200./)                                        ; bin1 = 0.:4.999, etc.)
  nlevels = dimsizes(levels)
;  colors  = span_color_rgba("NCV_jet",nlevels+1)
;   colors  = span_color_rgba("perc2_9lev",nlevels+1) 
   colors  = span_color_rgba("precip2_17lev",nlevels+1)
;---------------------------
;    Variable: Station_ID
;    Variable: Observation_Class
;    Variable: Observation_Type
;    Variable: Observation_Subtype
;    Variable: Latitude
;    Variable: Longitude
;    Variable: Station_Elevation
;    Variable: Pressure
;    Variable: Height
;    Variable: Time
;    Variable: Prep_QC_Mark
;    Variable: Prep_Use_Flag
;    Variable: Nonlinear_QC_Rel_Wgt
;    Variable: Analysis_Use_Flag
;    Variable: Errinv_Input
;    Variable: Errinv_Adjust
;    Variable: Errinv_Final
;    Variable: Observation
;    Variable: Obs_Minus_Forecast_adjusted
;    Variable: Obs_Minus_Forecast_unadjusted
;  fname  = "diag.nc"
;  a      = addfile(fname,"r")
;  lat    = a->Latitude 
;  lon    = a->Longitude 
;  R      = a->Observation
;https://www.ncl.ucar.edu/Applications/Scripts/fao56_1.ncl
   diri = "./"
   fili = systemfunc("cd "+diri+" ; ls pe*pm10_01.nc4")
   nfili= dimsizes(fili)
   print("nfili="+nfili)
   a    = addfiles(diri+fili,"r") 
   lat    = a[:]->Latitude
   lon    = a[:]->Longitude
   obs    = a[:]->Observation ;
   R      = a[:]->Obs_Minus_Forecast_unadjusted ;Observation
   aflag  = a[:]->Analysis_Use_Flag
   errinvf = a[:]->Errinv_Final
   errinvi = a[:]->Errinv_Input
   R = obs;-R
   printMinMax (R, False)
   ;R = mask(R,(aflag.eq.1), True)
   errinvf = mask(errinvf,(aflag.eq.1), True)
   errinvi = mask(errinvi,(aflag.eq.1), True)
   printMinMax (R, False)
  ;print(fili)
  ;print(lat)
  ;print(lon)
  ;print(R)
  npts= dimsizes(lat)
  print(npts)
  npts= dimsizes(lon)
  print(npts)
  npts= dimsizes(R)
  print(npts)
  N = num(.not.ismissing(R))
  print("Non-Missing: "+N)
  print("avgo= "+avg(R))
  std_f = stddev(R)
  print("stdo= "+std_f)
  print("avgerri= "+1/avg(errinvi))
  print("avgerrf= "+1/avg(errinvf))
;  npts = 100                                    ; Number of points.
;  lat  = random_uniform( 25., 50.,npts)         ; Create some dummy latitude
;  lon  = random_uniform(235.,290.,npts)         ; and longitude data that
                                                ; will contain the position of
                                                ; our markers.

;  R   = random_uniform(-1.2,35.,npts)   ; This is dummy data for determining
                                        ; how to color the markers.

;------------------------------ 
; Create X and Y arrays to hold the points for each range and initialize
; them to missing values.  We want to use num_distinct_markers
; different colors, so we need num_distinct_markers sets of X and
; Y points.
;
  num_distinct_markers = nlevels+1        ; number of distinct markers
  lat_new = new((/num_distinct_markers,npts/),double,-999)
  lon_new = new((/num_distinct_markers,npts/),double,-999)
;
; Group the points according to which range they fall in. At the
; same time, create the label that we will use later in the labelbar
;       
  do i = 0, num_distinct_markers-1
    if (i.eq.0) then
      indexes = ind(R.lt.levels(0))
    end if
    if (i.eq.num_distinct_markers-1) then
      indexes = ind(R.ge.max(levels))
    end if
    if (i.gt.0.and.i.lt.num_distinct_markers-1) then       
      indexes = ind(R.ge.levels(i-1).and.R.lt.levels(i))
    end if
;
; Now that we have the set of indexes whose values fall within 
; the given range, take the corresponding lat/lon values and store
; them, so later we can color this set of markers with the appropriate
; color.
;
    if (.not.any(ismissing(indexes))) then 
      npts_range = dimsizes(indexes)   ; # of points in this range.
      lat_new(i,0:npts_range-1) = lat(indexes)
      lon_new(i,0:npts_range-1) = lon(indexes)
    end if
    delete(indexes)            ; Necessary b/c "indexes" may be a different
                               ; size next time.
  end do

;===========================================================================
; Begin plotting section.
;
  wks = gsn_open_wks("png","obs_allpm10_2d")              ; send graphics to PNG file

;---Set up some map resources.
  mpres                       = True
  mpres@gsnMaximize           = True             ; Maximize plot in frame.
  mpres@gsnDraw               = False            ; Will draw later
  mpres@gsnFrame              = False            ; Don't advance the frame
  mpres@pmTickMarkDisplayMode = "Always"

;---Zoom in on United States.
  mpres@mpMinLatF             = 25.
  mpres@mpMaxLatF             = 55.
  mpres@mpMinLonF             = 235.
  mpres@mpMaxLonF             = 295.
  mpres@mpFillOn              = False ;True                ; turn on map fill
  mpres@mpLandFillColor       = "Background"        ; set to white  
  mpres@mpOutlineOn           = True                ; turn on map outline
  mpres@mpOutlineBoundarySets = "USStates"          ; turn on state boundaries
  mpres@tiMainString          = "OBS: PM10 (ug/Kg)" ;"Dummy station data colored and~C~sized according to range of values"
  map = gsn_csm_map(wks,mpres)

;--Create logical variables to hold the marker resources.
  gsres               = True
  gsres@gsMarkerIndex = 16          ; Use filled dots for markers.
;
; Loop through each grouping of markers, and draw them one set at
; a time, assigning the proper color and size with gsn_marker.
;
  base_size = 0.006
  pmid = new(num_distinct_markers,graphic)
  do i = 0, num_distinct_markers-1
    if (.not.ismissing(lat_new(i,0)))
      gsres@gsMarkerColor      = colors(i,:)
      gsres@gsMarkerSizeF      = base_size ;base_size * (i+1)/3.
      gsres@gsMarkerThicknessF = 0.7 ;0.7*(i+1)
      pmid(i) = gsn_add_polymarker(wks,map,lon_new(i,:),lat_new(i,:),gsres)
    end if
  end do

;---Draw labelbar and advance frame.
  add_labelbar(wks,map,levels,colors)
  draw(map)
  frame(wks)
end     
