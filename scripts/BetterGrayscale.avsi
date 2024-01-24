# Create proper grayscale using oklab color space: https://bottosson.github.io/posts/oklab/
# needs masktools2

#LoadPlugin("masktools2.dll")

function narrow_rgb_to_xyz(clip rgb){
    #values from my excel sheet
    x=mt_lutxyz(rgb.ExtractR(),rgb.ExtractG(),rgb.ExtractB(),expr=mt_polish("0.412410846*x+0.357584568*y+0.180453804*z"),Y=3,U=0,V=0,use_expr=2)
    y=mt_lutxyz(rgb.ExtractR(),rgb.ExtractG(),rgb.ExtractB(),expr=mt_polish("0.212649343*x+0.715169136*y+0.072181522*z"),Y=3,U=0,V=0,use_expr=2)
    z=mt_lutxyz(rgb.ExtractR(),rgb.ExtractG(),rgb.ExtractB(),expr=mt_polish("0.019331758*x+0.119194856*y+0.950390034*z"),Y=3,U=0,V=0,use_expr=2)
    MergeRGB(x,y,z)
}
function xyz_to_narrow_rgb(clip xyz){
    #values from my excel sheet
    r=mt_lutxyz(xyz.ExtractR(),xyz.ExtractG(),xyz.ExtractB(),expr=mt_polish("3.240812399*x-1.537308446*y-0.498586523*z"),Y=3,U=0,V=0,use_expr=2)
    g=mt_lutxyz(xyz.ExtractR(),xyz.ExtractG(),xyz.ExtractB(),expr=mt_polish("-0.969243017*x+1.875966303*y+0.041555031*z"),Y=3,U=0,V=0,use_expr=2)
    b=mt_lutxyz(xyz.ExtractR(),xyz.ExtractG(),xyz.ExtractB(),expr=mt_polish("0.055638398*x-0.204007461*y+1.05712957*z"),Y=3,U=0,V=0,use_expr=2)
    MergeRGB(r,g,b)
}

function wide_rgb_to_xyz(clip rgb){
    #values from my excel sheet
    x=mt_lutxyz(rgb.ExtractR(),rgb.ExtractG(),rgb.ExtractB(),expr=mt_polish("0.636103731*x+0.145005373*y+0.169340114*z"),Y=3,U=0,V=0,use_expr=2)
    y=mt_lutxyz(rgb.ExtractR(),rgb.ExtractG(),rgb.ExtractB(),expr=mt_polish("0.262402691*x+0.678462002*y+0.059135307*z"),Y=3,U=0,V=0,use_expr=2)
    z=mt_lutxyz(rgb.ExtractR(),rgb.ExtractG(),rgb.ExtractB(),expr=mt_polish("0.000049510*x+0.028316941*y+1.060550197*z"),Y=3,U=0,V=0,use_expr=2)
    MergeRGB(x,y,z)
}
function xyz_to_wide_rgb(clip xyz){
    #values from my excel sheet
    r=mt_lutxyz(xyz.ExtractR(),xyz.ExtractG(),xyz.ExtractB(),expr=mt_polish("1.719288549*x-0.356830388*y-0.254625614*z"),Y=3,U=0,V=0,use_expr=2)
    g=mt_lutxyz(xyz.ExtractR(),xyz.ExtractG(),xyz.ExtractB(),expr=mt_polish("-0.66649797*x+1.615688601*y+0.016331712*z"),Y=3,U=0,V=0,use_expr=2)
    b=mt_lutxyz(xyz.ExtractR(),xyz.ExtractG(),xyz.ExtractB(),expr=mt_polish("0.017715392*x-0.04312261*y+0.942482633*z"),Y=3,U=0,V=0,use_expr=2)
    MergeRGB(r,g,b)
}

function xyz_to_lms(clip xyz){
    #values from oklab
    l=mt_lutxyz(xyz.ExtractR(),xyz.ExtractG(),xyz.ExtractB(),expr=mt_polish("0.8189330101*x+0.3618667424*y-0.1288597137*z"),Y=3,U=0,V=0,use_expr=2)
    m=mt_lutxyz(xyz.ExtractR(),xyz.ExtractG(),xyz.ExtractB(),expr=mt_polish("0.0329845436*x+0.9293118715*y+0.0361456387*z"),Y=3,U=0,V=0,use_expr=2)
    s=mt_lutxyz(xyz.ExtractR(),xyz.ExtractG(),xyz.ExtractB(),expr=mt_polish("0.0482003018*x+0.2643662691*y+0.6338517070*z"),Y=3,U=0,V=0,use_expr=2)
    MergeRGB(l,m,s)
}
function lms_to_xyz(clip lms){
    #values from inverse calculation by my excel sheet
    x=mt_lutxyz(lms.ExtractR(),lms.ExtractG(),lms.ExtractB(),expr=mt_polish("1.227013851*x-0.557799981*y+0.281256149*z"),use_expr=2)
    y=mt_lutxyz(lms.ExtractR(),lms.ExtractG(),lms.ExtractB(),expr=mt_polish("-0.040580178*x+1.11225687*y-0.071676679*z"),use_expr=2)
    z=mt_lutxyz(lms.ExtractR(),lms.ExtractG(),lms.ExtractB(),expr=mt_polish("-0.076381285*x-0.421481978*y+1.58616322*z"),use_expr=2)
    MergeRGB(x,y,z)
}

function lms_to_oklab(clip lms){
    #values from oklab
    lms=lms.mt_lut(expr=mt_polish("x^(1/3)"),Y=3,U=3,V=3,use_expr=2)
    L=mt_lutxyz(lms.ExtractR(),lms.ExtractG(),lms.ExtractB(),expr=mt_polish("0.2104542553*x+0.7936177850*y-0.0040720468*z"),use_expr=2)
    a=mt_lutxyz(lms.ExtractR(),lms.ExtractG(),lms.ExtractB(),expr=mt_polish("1.9779984951*x-2.4285922050*y+0.4505937099*z"),use_expr=2)
    b=mt_lutxyz(lms.ExtractR(),lms.ExtractG(),lms.ExtractB(),expr=mt_polish("0.0259040371*x+0.7827717662*y-0.8086757660*z"),use_expr=2)
    MergeRGB(L,a,b)
}
function oklab_to_lms(clip Lab){
    #values from inverse calculation by my excel sheet
    l=mt_lutxyz(Lab.ExtractR(),Lab.ExtractG(),Lab.ExtractB(),expr=mt_polish("1.000000000*x+0.396337792*y+0.215803758*z"),Y=3,U=0,V=0,use_expr=2)
    m=mt_lutxyz(Lab.ExtractR(),Lab.ExtractG(),Lab.ExtractB(),expr=mt_polish("1.000000000*x-0.105561342*y-0.063854175*z"),Y=3,U=0,V=0,use_expr=2)
    s=mt_lutxyz(Lab.ExtractR(),Lab.ExtractG(),Lab.ExtractB(),expr=mt_polish("1.000000000*x-0.089484182*y-1.291485538*z"),Y=3,U=0,V=0,use_expr=2)
    MergeRGB(l,m,s).mt_lut(expr=mt_polish("x^3"),Y=3,U=3,V=3,use_expr=2)
}

function xyzGrayscale(clip xyz){
    lms=xyz_to_lms(xyz)
    Lab=lms_to_oklab(lms)
    Lab=Lab.mt_lut(expr="x",Y=3,U=0,V=0,use_expr=2)
    lms=oklab_to_lms(Lab)
    lms_to_xyz(lms)
}

function betterGrayscale(clip org, string "primariesIn", string "transferIn", string "scaleIn", string "primariesOut", string "transferOut", string "scaleOut", string "pixelTypeOut", int "outWidth", int "outHeight"){
    primariesOut = default( primariesOut, primariesIn )
    transferOut = default( transferOut, transferIn )
    scaleOut = default( scaleOut, scaleIn )
    pixelTypeOut = default( pixelTypeOut, org.PixelType() )
    outWidth = default( outWidth, org.Width() )
    outHeight = default( outHeight, org.Height() )
    chromaLoc = (transferIn=="st2084" || transferIn=="std-b67")?"top_left":"left"
    
    rgb=org.IsRGB()?org:org.z_ConvertFormat(pixel_type="RGBPS",colorspace_op=primariesIn+":"+transferIn+":"+primariesIn+":"+scaleIn+"=>rgb:linear:"+primariesIn+":full",dither_type="none",resample_filter="bilinear",resample_filter_uv="bilinear", chromaloc_op=chromaLoc+"=>center", width=outWidth, height=outHeight) 
    xyz=primariesIn=="2020"?wide_rgb_to_xyz(rgb):narrow_rgb_to_xyz(rgb)
    xyz=xyz.xyzGrayscale()
    rgb=primariesIn=="2020"?xyz_to_wide_rgb(xyz):xyz_to_narrow_rgb(xyz)
    org.IsRGB()?rgb:rgb.z_ConvertFormat(pixel_type=pixelTypeOut,colorspace_op="rgb:linear:"+primariesIn+":full=>"+primariesOut+":"+transferOut+":"+primariesOut+":"+scaleOut,dither_type="none",resample_filter="bilinear",resample_filter_uv="bilinear", chromaloc_op="center=>"+chromaLoc)
}

function betterNtscGrayscale(clip org){
    betterGrayscale(org, "170m", "170m", "limited")
}
function betterPalGrayscale(clip org){
    betterGrayscale(org, "470bg", "470bg", "limited")
}
function better709Grayscale(clip org){
    betterGrayscale(org, "709", "709", "limited")
}
function better2020Grayscale(clip org){
    betterGrayscale(org, "2020", "2020", "limited")
}
function betterPqGrayscale(clip org){
    betterGrayscale(org, "2020", "st2084", "limited")
}
function betterHlgGrayscale(clip org){
    betterGrayscale(org, "2020", "std-b67","limited")
}
