# scene change detection, strongly inspired by: https://forum.doom9.org/showthread.php?t=170532
# needs MVTools2 and avsresize
# writes the scene files which can then be read by DoViSceneFileReader

LoadPlugin("DGDecodeNV.dll")
DGSource("bl.dgi")
crop(0, 280, 0, -280)

#LoadPlugin("ffms2.dll")
#FFVideoSource("bl.ts", cachefile="index.ffindex", threads=1, track=0) 
#crop(0, 280, 0, -280)

#LoadPlugin("LSMASHSource.dll")
#LWLibavVideoSource("bl.ts", format="YUV420P8", prefer_hw=0, stream_index=0)
#crop(0, 280, 0, -280)

global LumaDiffFac = 4.5 # 4.0 for SDR, 4.5 for HDR
global LumaDiffMin = 7.8 # 11.0 for SDR, 8.0 for HDR (the usefull range for this is somewhere between 7 and 9 for HDR)

scdResizeFac = (width()<1100)?1:(width()<2100)?2:4
cllResizeFac = (width()<1100)?1:(width()<2100)?1:2
global scd=z_ConvertFormat(width=width()/scdResizeFac,height=height()/scdResizeFac,pixel_type="Y8",dither_type="none",resample_filter="bilinear",resample_filter_uv="bilinear",chromaloc_op="top_left=>top_left",colorspace_op="2020ncl:st2084:2020:limited=>2020ncl:st2084:2020:full")
global rgb=z_ConvertFormat(width=width()/cllResizeFac,height=height()/cllResizeFac,pixel_type="RGBP16",dither_type="none",resample_filter="bilinear",resample_filter_uv="bilinear",chromaloc_op="top_left=>center",colorspace_op="2020ncl:st2084:2020:limited=>rgb:st2084:2020:full")
scd_super = scd.MSuper(pel=2, sharp=0, hpad=16, vpad=16)
bac2_vec = MAnalyse(scd_super, isb=true, delta=2, blksize=16, search=4, searchparam=2)
back_vec = MAnalyse(scd_super, isb=true, delta=1, blksize=16, search=4, searchparam=2)
fore_vec = MAnalyse(scd_super, isb=false, delta=1, blksize=16, search=4, searchparam=2)

global bac2_compensated=scd.MCompensate(scd_super, bac2_vec, thSAD=-1, thSCD1=350, thSCD2=256) 
global back_compensated=scd.MCompensate(scd_super, back_vec, thSAD=-1, thSCD1=350, thSCD2=256) 
global fore_compensated=scd.MCompensate(scd_super, fore_vec, thSAD=-1, thSCD1=350, thSCD2=256) 

global writeCheck = """(YDifferenceToNext(scd)>LumaDiffMin+0.5) && LumaDifference(back_compensated,scd)>LumaDiffMin && (Min(LumaDifference(back_compensated,scd),LumaDifference(bac2_compensated,scd))/LumaDifference(fore_compensated,scd)>LumaDiffFac) || (AverageLuma(scd)<9 && YPlaneMinMaxDifference(scd)<44 && YDifferenceFromPrevious(scd)<0.01 && YDifferenceToNext(scd)>0.01)"""
#                                                                                                                       ^ scene change when drastic changes in luma-difference, except it the change returns back one frame later.
#                       ^ no scene change when luma-difference is too low, while YDifferenceToNext processes much faster then LumaDifference
#                                                                                                                                                                                                                                                            ^ very dark images with low variations and very low differences to previous and next images are fade ins and fade outs, these must also be scene changes
global newSceneNextFrame = "(" + writeCheck + ")?1:0"
scd

#Uncomment following line for troubleshooting. It puts the scene change threshold numbers on the screen.
#StackHorizontal(fore_compensated,scd,back_compensated)
#ScriptClip("Subtitle(String("+writeCheck+"),align=8)")
#ScriptClip("""Subtitle(String(YDifferenceFromPrevious(scd))+" "+String(YDifferenceToNext(scd))+" "+String(YDifferenceToNext(scd)/YDifferenceFromPrevious(scd)),align=4)""")
#ScriptClip("""Subtitle(String(LumaDifference(fore_compensated,scd))+" "+String(LumaDifference(back_compensated,scd))+" "+String(LumaDifference(bac2_compensated,scd))+" "+String(LumaDifference(back_compensated,scd)/LumaDifference(fore_compensated,scd))+" "+String(LumaDifference(bac2_compensated,scd)/LumaDifference(fore_compensated,scd)),align=6)""")
#ScriptClip("""Subtitle(String(AverageLuma(scd))+" "+String(YPlaneMinMaxDifference(scd)),align=2)""")
#MShow(scd_super,fore_vec,thSCD1=350,thSCD2=256)

#This line sets the frame property "_SceneChangeNext"
#ScriptClip(writeCheck+"""?propSet("_SceneChangeNext",1):propSet("_SceneChangeNext",0)""")

#This line writes frame numbers of scene changes to file. These can be imported into Vegas or other editor
WriteFileIf("sceneCuts.txt", writeCheck, "current_frame+1", flush=false, append=false)

function BitShiftLeft(int value, int shift){
    # bi-directional left shift function with rounding to nearest when downshifting
    (shift < 0) ? BitRShiftU(value + BitLShift(1,-shift-1), -shift) : BitLShift(value, shift)
}
function BitShiftRight(int value, int shift){
    # bi-directional right shift function with rounding to nearest when downshifting
    (shift > 0) ? BitRShiftU(value + BitLShift(1,shift-1), shift) : BitLShift(value, -shift)
}

#This line writes max RGB pixel value for each frame
WriteFile("maxPq.txt", "current_frame", """ " " """, "BitShiftLeft(Max(RPlaneMax(rgb),GPlaneMax(rgb),BPlaneMax(rgb)),12-BitsPerComponent(rgb))", flush=false, append=false)

