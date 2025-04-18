;----------------------------------------------------------------------
; station_2.ncl
;----------------------------------------------------------------------
;
; Concepts illustrated:
;   - Drawing markers on a map indicating the locations of station data
;   - Generating dummy data using "random_uniform"
;   - Drawing markers of different sizes and colors on a map
;   - Drawing a custom legend outside of a map plot
;   - Attaching a custom labelbar to a plot
;
; This example creates some dummy station data, and then plots each
; value by coloring and sizing it depending on which range of values 
; it falls in.
;
;
; This creates one with a labelbar.
;----------------------------------------------------------------------
;
; These files are loaded by default in NCL V6.2.0 and newer
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"


;----------------------------------------------------------------------
; This function attaches a labelbar to the given plot.
;----------------------------------------------------------------------
function attach_labelbar(wks,map,arr[*]:numeric,colors[*])
local lbres, vph, vpw, nboxes
begin
  getvalues map
    "vpHeightF" : vph
    "vpWidthF"  : vpw
  end getvalues

  nboxes = dimsizes(colors)
   
  lbres                    = True          ; labelbar only resources
  lbres@lbAutoManage       = False          ; Necessary to control sizes
  lbres@lbFillColors       = colors
  lbres@vpWidthF           = 0.7 * vpw     ; labelbar width
  lbres@vpHeightF          = 0.2 * vph     ; labelbar height
  lbres@lbMonoFillPattern  = True          ; Solid fill pattern
  lbres@lbLabelFontHeightF = 0.01          ; font height. default is small
  lbres@lbOrientation      = "horizontal"
  lbres@lbPerimOn          = False
  lbres@lbLabelAlignment   = "InteriorEdges"
  ;lbid = gsn_create_labelbar(wks,nboxes,""+toint(arr),lbres)
  lbid = gsn_create_labelbar(wks,nboxes,""+arr,lbres)
  amres                  = True
  amres@amJust           = "TopCenter"
  amres@amParallelPosF   =  0.0    ; Center
  amres@amOrthogonalPosF =  0.6    ; Bottom
  annoid = gsn_add_annotation(map,lbid,amres)
  return(annoid)
end

;----------------------------------------------------------------------
; Main code
;----------------------------------------------------------------------
begin
;---Set some needed arrays
  ;arr = (/0.,5.,10.,15.,20.,23.,26./)    ; bin settings (bin0 = < 0., 
  arr = fspan(-10,10,21)*0.03                                       ; bin1 = 0.:4.999, etc.)
  narr = dimsizes(arr)
  labels = new(narr+1,string)   ; Labels for legend.

;  npts = 100                                    ; Number of points.
;  lat  = random_uniform( 25., 50.,npts)         ; Create some dummy latitude
;  lon  = random_uniform(235.,290.,npts)         ; and longitude data that
                                                ; will contain the position of
                                                ; our markers.
;  R   = random_uniform(-1.0,10.,npts)   ; This is dummy data for determining
                                        ; how to color the markers.
  fname  = "fort.1001"
  nrows = numAsciiRow(fname)
  ncols = numAsciiCol(fname)

  ;print("'" + fname + "' has " + nrows + " rows and " + ncols + \
  ;      " columns of data.")
  ;data = new ( (/nrows, ncols/), float )

  dataobs = asciiread("fort.1001",(/nrows, ncols/),"float")
  printVarSummary(dataobs)
  lat  = dataobs(:,0)
  lon  = dataobs(:,1)
  prs  = dataobs(:,2)
 ; print(lat(0:2))
 ; print(lon(0:2))
 ; print(max(prs))
 ; print(min(prs))
  fso    = dataobs(:,5) 
  print(avg(fso))
  print(stddev(fso))
  ;print(max(fso))
  ;print(min(fso))
  xscl=toint(stddev(fso))
  xpos=where(fso.le.0.0,-1,-999)
  xpos@_FillValue = -999
  print(sum(xpos))
  xneg=where(fso.gt.0.0,1,-999)
  xneg@_FillValue = -999
  print(sum(xneg))
;  R = where(prs.gt.800,fso,-999)
; R = where(prs.lt.100,fso,-999)
R = where(prs.ge.475.and.prs.ge.475,fso,-999)
  R@_FillValue = -999
  print(avg(R))
  xpos=where(R.lt.0.0,-1,-999)
  xpos@_FillValue = -999
  print(sum(xpos))
  xneg=where(R.gt.0.0,1,-999)
  xneg@_FillValue = -999
  print(sum(xneg))
  print(stddev(R))
;
; Create X and Y arrays to hold the points for each range and initialize
; them to missing values.  We want to use num_distinct_markers
; different colors, so we need num_distinct_markers sets of X and
; Y points.
;
  num_distinct_markers = dimsizes(arr)+1        ; number of distinct markers
  lat_new = new((/num_distinct_markers,dimsizes(R)/),float,-999)
  lon_new = new((/num_distinct_markers,dimsizes(R)/),float,-999)

;---Group the points according to which range they fall in.
  do i = 0, num_distinct_markers-1
    if (i.eq.0) then
      indexes = ind(R.lt.arr(0))
    end if
    if (i.eq.num_distinct_markers-1) then
      indexes = ind(R.ge.max(arr))
    end if
    if (i.gt.0.and.i.lt.num_distinct_markers-1) then       
      indexes = ind(R.ge.arr(i-1).and.R.lt.arr(i))
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

;----------------------------------------------------------------------
; Begin plotting section.
;----------------------------------------------------------------------
  wks = gsn_open_wks("x11","fso_u")   ; send graphics to PNG file

;---Set up some map resources.
  mpres              = True
  mpres@gsnMaximize  = True             ; Maximize plot in frame.
  mpres@gsnFrame     = False            ; Don't advance the frame
  mpres@gsnDraw      = False            ; Don't advance the frame

;---Zoom in on United States.
  mpres@mpMinLatF    = -90.
  mpres@mpMaxLatF    = 90.
  mpres@mpMinLonF    = 0.
  mpres@mpMaxLonF    = 359.
  mpres@mpCenterLonF = 0 

  mpres@mpFillColors = (/"transparent","transparent","gray","transparent"/)

  ;mpres@tiMainString = "Dummy station data colored according to range of values"
  map = gsn_csm_map(wks,mpres)

  gsres               = True
  gsres@gsMarkerIndex = 16          ; Use filled dots for markers.
;---Get nice spacing through color map for marker colors
  getvalues wks
    "wkColorMapLen" : clen     ; number of colors in color map
  end getvalues

  nstep = (clen-2)/narr
  colors = ispan(2,clen-1,nstep)

;---Loop through each "bin" and attach the markers to the map.
  do i = 0, num_distinct_markers-1
    if (.not.ismissing(lat_new(i,0)))
      gsres@gsMarkerColor      = colors(i)
      dumstr = unique_string("marker")
      map@$dumstr$ = gsn_add_polymarker(wks,map,lon_new(i,:),lat_new(i,:),gsres)
    end if
  end do

;----------------------------------------------------------------------
; First version of this plot has a legend at the bottom
;----------------------------------------------------------------------
  ;draw(map)                  ; Drawing the map will draw the markers
  ;draw_legend(wks,lat_new,lon_new,arr,colors) ; Draw a legend at the bottom

  ;frame(wks)    ; Advance the frame 

;----------------------------------------------------------------------
; Second version of this plot has a labelbar added.
;----------------------------------------------------------------------
  lbid = attach_labelbar(wks,map,arr,colors)   ; Attach a labelbar

  draw(map)     ; Drawing the map will draw everything
  frame(wks)    ; Advance the frame.
end     
