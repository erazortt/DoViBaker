bl=FFVideoSource("clip.ts", threads=1, track=0)
el=FFVideoSource("clip.ts", threads=1, track=1)

#if downsizing by factor of two (4K->FullHD) is wished, it should be done here
bl=bl.z_ConvertFormat(chromaloc_op="top_left=>top_left",width=bl.Width()/2,height=bl.Height()/2,resample_filter="spline36",resample_filter_uv="spline36")

DoViBaker(bl,el,cubes_basepath="C:\lut_",cubes="1000.cube;1414.cube;2000.cube;2828.cube;4000.cube;5656.cube;8000.cube;10000.cube",mclls="1000;1414;2000;2828;4000;5656;8000")
z_ConvertFormat(chromaloc_op="center=>top_left",pixel_type="YUV420P16",colorspace_op="rgb:std-b67:2020:full=>2020ncl:std-b67:2020:limited",dither_type="none",resample_filter="spline36",resample_filter_uv="spline36")

ScriptClip("""
mcll=propGetInt("_dovi_max_content_light_level")
mpq=propGetInt("_dovi_max_pq")
subtitle("maxpq = " + string(mpq) + "; maxcll = " + string(mcll))
""")
