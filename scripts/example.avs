bl=FFVideoSource("clip.ts", threads=1, track=0)
el=FFVideoSource("clip.ts", threads=1, track=1)

#if downsizing by factor of two (4K->FullHD) is wished, it should be done here
bl=bl.z_ConvertFormat(chromaloc_op="top_left=>top_left",width=bl.Width()/2,height=bl.Height()/2,resample_filter="spline36",resample_filter_uv="spline36")
#the el clip does not need downsizing since its usually already at that size

DoViBaker(bl,el)
DoViTonemap(masterMaxNits=-1, masterMinNits=0, targetMaxNits=1000, targetMinNits=0, lumScale=1.0).AVSCube("PQ2HLG.cube")
z_ConvertFormat(chromaloc_op="center=>top_left",pixel_type="YUV420P16",colorspace_op="rgb:std-b67:2020:full=>2020ncl:std-b67:2020:limited",dither_type="none",resample_filter="spline36",resample_filter_uv="spline36")

ScriptClip("""
mcll=propGetInt("_dovi_dynamic_max_content_light_level")
mpq=propGetInt("_dovi_dynamic_max_pq")
subtitle("maxpq = " + string(mpq) + "; maxcll = " + string(mcll))
""")
