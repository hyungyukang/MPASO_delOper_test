;======================================================================
; mpas_1.ncl
;
; Concepts illustrated:
;   - Plotting MPAS data
;   - Plotting unstructured data
;   - Plotting data with missing lat/lon coordinates
;   - Using cnFillPalette to assign a color palette to contours
;======================================================================
;
; These files are loaded by default in NCL V6.2.0 and newer
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin
;---Read data from MPAS Grid
    f  = addfile("output.nc","r")
    lonCell = f->lonCell
    latCell = f->latCell

;---Convert to degrees from radians
    RAD2DEG = 180.d0/(4.0*atan(1.0))
    lonCell = lonCell*RAD2DEG
    latCell = latCell*RAD2DEG

;---Get time dimension
;   fdims = getfiledimsizes(f)
;   nt = fdims(2)
;   print(nt)

;****************************
;****************************

;---Start the graphics
    ifile = "fig"
    wks = gsn_open_wks("png",ifile)

    sp1 = f->num
    sp2 = f->anl
    sp3 = f->anl

    sp3 = log10(abs(sp1 - sp2))

;---levels for contours
    nlevels   = 251
;   levels    = fspan(min(sp3),max(sp3),nlevels)
    levels    = fspan(-5,0,nlevels)
;   levels    = fspan(-1.0,1.0,nlevels)
;   levels    = fspan(-4,0.2,nlevels)
;   levels    = 10^levels

    res                      = True              ; Plot mods desired.
    res@gsnMaximize          = True              ; Maximize plot

;   res@mpDataBaseVersion     = "MediumRes"
;   res@mpLandFillColor       = "wheat2"
;   res@mpOceanFillColor      = "transparent"   ; no fill
;   res@mpGridAndLimbOn       = False
    res@mpOutlineOn           = False
;   res@mpFillDrawOrder       = "PostDraw"

    res@cnFillMode           = "RasterFill"      ; turn raster on      
;   res@cnFillMode           = "AreaFill"      ; turn raster on      
    res@cnFillOn             = True              ; color plot desired
    res@cnLinesOn            = True              ; turn off contour lines
    res@cnLineLabelsOn       = False             ; turn off contour labels

    res@cnLevelSelectionMode  = "ExplicitLevels"
    res@cnLevels              = levels
;   res@cnFillPalette         = "BlGrYeOrReVi200"
    res@cnFillPalette         = "BkBlAqGrYeOrReViWh200"
;   res@cnFillPalette         = "WhiteBlueGreenYellowRed"
;   res@cnFillPalette         = "hotres"
;   res@cnFillPalette         = "GMT_ocean"
;   res@cnFillPalette         = "BlueWhiteOrangeRed"
    res@cnLinesOn            = False             ; turn off contour lines

;   res@cnLevelSelectionMode = "ManualLevels"
;   res@cnMinLevelValF       = 0
;   res@cnMaxLevelValF       = 0.5
;   res@cnLevelSpacingF      = 0.005
    res@lbBoxLinesOn         = False

    res@lbLabelStrings        = sprintf("%.2f",levels)  ; format decimal places of labels
    res@lbLabelStride         = 25                ; Over 250 levels, so only plot every 25th level
  
    res@pmLabelBarWidthF      = 0.75              ; make labelbar slightly longer
    res@pmLabelBarHeightF     = 0.08              ;   and narrower
    res@pmLabelBarOrthogonalPosF = 0.1            ; move labelbar down

    res@gsnAddCyclic         = False

    res@sfXArray             = lonCell
    res@sfYArray             = latCell
    res@gsnAddCyclic         = False

    res@vpXF = 0.0 
    res@vpYF = 0.6 
    res@vpWidthF = 1.0 
    res@vpHeightF = 0.6 

    plot = gsn_csm_contour_map(wks,sp3,res)  

    sfile = ifile+".png"
    afile = "trim_"+ifile+".png"

    system("convert "+" "+sfile+" -trim   "+" "+afile)

    exit

end
