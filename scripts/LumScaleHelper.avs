# compares the base brightness of a pq stream with a sdr stream.
# needs masktools2, avsresize and BetterGrayscale.avsi
# you might need to use the trim command so that the frames of the two streams match exactly

#LoadPlugin("avsresize.dll")
#LoadPlugin("masktools2.dll")
#Import("BetterGrayscale.avsi")


#### Select one of the source libraries ####

#LoadPlugin("DGDecodeNV.dll")
#pq=DGSource("pq.dgi")
#sdr=DGSource("sdr.dgi")

#LoadPlugin("ffms2.dll")
#pq=FFVideoSource("pq.ts", cachefile="pq.index.ffindex", threads=1, track=0)
#sdr=FFVideoSource("sdr.ts", cachefile="sdr.index.ffindex", threads=1, track=0)

#LoadPlugin("LSMASHSource.dll")
#pq=LWLibavVideoSource("pq.ts", format="YUV420P10", prefer_hw=0, stream_index=0)
#sdr=LWLibavVideoSource("pq.ts", format="YUV420P10", prefer_hw=0, stream_index=0)

############################################


resizeFac=(pq.width()<1100)?1:(pq.width()<2100)?2:4
pqg=pq.betterGrayscale(primariesIn="2020", transferIn="st2084", scaleIn="limited", scaleOut="full", pixelTypeOut="Y32", outWidth=pq.width()/resizeFac, outHeight=pq.height()/resizeFac)
sdg=sdr.betterGrayscale(primariesIn="709", transferIn="709", scaleIn="limited", scaleOut="full", pixelTypeOut="Y32", outWidth=pq.width()/resizeFac, outHeight=pq.height()/resizeFac)
pqLin=mt_lut(pqg,mt_polish("(x>0)?100*(max(x^0.012683314-0.8359375,0)/(18.8515625-18.6875*(x^0.012683314)))^6.277394636:0"),use_expr=2)
sdLin=mt_lut(sdg,mt_polish("(x>0)?2.676*(x)^2.4:0"),use_expr=2)
global pqc=mt_lutxy(pqLin,sdg,mt_polish("(y>0&y<0.5)?x:0"),use_expr=2)
global sdc=mt_lutxy(sdLin,sdg,mt_polish("(y>0&y<0.5)?x:0"),use_expr=2)

pqLinMultiplied=mt_lut(pqLin,use_expr=2,mt_polish("x*1.00")) #<---- change multiplication factor here to the shown value of lumScale to see the resulting brightness comparison

rePq=mt_lut(pqLinMultiplied,mt_polish("(x<0.5)?max(x/2.676,0)^(1/2.4):0.5"),use_expr=2).Subtitle("multiplied PQ")
reSdr=mt_lut(sdLin,mt_polish("(x<0.5)?max(x/2.676,0)^(1/2.4):0.5"),use_expr=2).Subtitle("sdr")
Interleave(reSdr,rePq).SelectEvery(4,0,1).ScriptClip("""Subtitle(align=8,"lumScale="+String(AverageLuma(pqc)>0?Float(AverageLuma(sdc))/AverageLuma(pqc):1.0))""")