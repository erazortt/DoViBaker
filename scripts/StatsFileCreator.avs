# scene change detection, strongly inspired by: https://forum.doom9.org/showthread.php?t=170532
# writes the scene files which can then be read by DoViSceneFileReader
# cropping must be set very precisely
# needs MVTools2 and avsresize

#LoadPlugin("avsresize.dll")
#LoadPlugin("mvtools2.dll")


#### Select one of the source libraries ####

#LoadPlugin("DGDecodeNV.dll")
#DGSource("bl.dgi")
#crop(0, 280, 0, -280)

#LoadPlugin("ffms2.dll")
#FFVideoSource("clip.ts", cachefile="index.ffindex", threads=1, track=0) 
#crop(0, 280, 0, -280)

#LoadPlugin("LSMASHSource.dll")
#LWLibavVideoSource("clip.ts", format="YUV420P10", prefer_hw=0, stream_index=0)
#crop(0, 280, 0, -280)

############################################


global LumaDiffFac = 4.5 # 4.0 for SDR, 4.5 for HDR
global LumaDiffMin = 7.8 # 11.0 for SDR, 8.0 for HDR (the usefull range for this is somewhere between 7 and 9 for HDR)

resizeFac = (width()<1100)?1:(width()<2100)?2:4
global scd=z_ConvertFormat(width=width()/resizeFac,height=height()/resizeFac,pixel_type="Y8",dither_type="none",resample_filter="bilinear",resample_filter_uv="bilinear",chromaloc_op="top_left=>top_left",colorspace_op="2020ncl:st2084:2020:limited=>2020ncl:st2084:2020:full")
global rgb=z_ConvertFormat(width=width()/resizeFac,height=height()/resizeFac,pixel_type="RGBP16",dither_type="none",resample_filter="bilinear",resample_filter_uv="bilinear",chromaloc_op="top_left=>center",colorspace_op="2020ncl:st2084:2020:limited=>rgb:st2084:2020:full")
scd_super = scd.MSuper(pel=2, sharp=0, hpad=16, vpad=16)
bac2_vec = MAnalyse(scd_super, isb=true, delta=2, blksize=16, search=4, searchparam=2)
back_vec = MAnalyse(scd_super, isb=true, delta=1, blksize=16, search=4, searchparam=2)
fore_vec = MAnalyse(scd_super, isb=false, delta=1, blksize=16, search=4, searchparam=2)

global bac2_compensated=scd.MCompensate(scd_super, bac2_vec, thSAD=-1, thSCD1=350, thSCD2=256) 
global back_compensated=scd.MCompensate(scd_super, back_vec, thSAD=-1, thSCD1=350, thSCD2=256) 
global fore_compensated=scd.MCompensate(scd_super, fore_vec, thSAD=-1, thSCD1=350, thSCD2=256) 

global newSceneNextFrame = """(YDifferenceToNext(scd)>LumaDiffMin+0.5) && LumaDifference(back_compensated,scd)>LumaDiffMin && (Min(LumaDifference(back_compensated,scd),LumaDifference(bac2_compensated,scd))/LumaDifference(fore_compensated,scd)>LumaDiffFac) || ((AverageLuma(scd)-YPlaneMin(scd))<1.5 && YPlaneMinMaxDifference(scd)<44 && YDifferenceFromPrevious(scd)<0.0075 && YDifferenceToNext(scd)>0.0075)"""
#                                                                                                                              ^ scene change when drastic changes in luma-difference, except when the change returns back one frame later.
#                              ^ no scene change when luma-difference is too low. YDifferenceToNext processes much faster then LumaDifference thus it is checked first as a filter
#                                                                                                                                                                                                                                                                   ^ very dark images with low variations and very low differences to previous and next images are fade ins and fade outs, these must also be scene changes
global newSceneNextFrameInt = "(" + newSceneNextFrame + ")?1:0"
scd

#Uncomment following line for troubleshooting. It puts the scene change threshold numbers on the screen.
#StackHorizontal(fore_compensated,scd,back_compensated)
#ScriptClip("Subtitle(String("+newSceneNextFrame+"),align=8)")
#ScriptClip("""Subtitle(String(YDifferenceFromPrevious(scd))+" "+String(YDifferenceToNext(scd))+" "+String(YDifferenceToNext(scd)/YDifferenceFromPrevious(scd)),align=4)""")
#ScriptClip("""Subtitle(String(LumaDifference(fore_compensated,scd))+" "+String(LumaDifference(back_compensated,scd))+" "+String(LumaDifference(bac2_compensated,scd))+" "+String(LumaDifference(back_compensated,scd)/LumaDifference(fore_compensated,scd))+" "+String(LumaDifference(bac2_compensated,scd)/LumaDifference(fore_compensated,scd)),align=6)""")
#ScriptClip("""Subtitle(String(AverageLuma(scd))+" "+String(YPlaneMinMaxDifference(scd)),align=2)""")
#MShow(scd_super,fore_vec,thSCD1=350,thSCD2=256)

function signal2pq(clip c, int value){
    Int(value/Float(BitLShift(1,BitsPerComponent(c))-1)*4095)
}

#This line writes the frame numbers where scene changes occur to file. These can be imported into Vegas or other editor
#WriteFileIf("sceneCuts.txt", newSceneNextFrame, "current_frame+1", flush=false, append=false)

#This line writes max RGB pixel value, scene change detection and max frame brightness for each frame
WriteFile("statsFile.txt", "current_frame", """ " " """, "Eval(newSceneNextFrameInt)", """ " " """, "rgb.signal2pq(Max(RPlaneMax(rgb),GPlaneMax(rgb),BPlaneMax(rgb)))", """ " " """, "rgb.signal2pq(Min(RPlaneMin(rgb),GPlaneMin(rgb),BPlaneMin(rgb)))", flush=false, append=false)

