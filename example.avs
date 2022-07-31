bl=DGSource("C:\blclip.dgi")
el=DGSource("C:\elclip.dgi")

bl=bl.Spline36Resize(1920,1080) #if downsizing by factor of two (4K->FullHD) is wished, it should be done here
el=el                           #since the el clip is usually already only at FullHD size, no downsizing is needed here

DoViBaker(bl,el,rpu="C:\RPU.bin",cubes_basepath="C:\lut_1.75_",cubes="1000.cube;1414.cube;2000.cube;2828.cube;4000.cube;5656.cube;8000.cube;10000.cube;",mclls="1000;1414;2000;2828;4000;5656;8000")
z_ConvertFormat(pixel_type="YUV420P16",colorspace_op="rgb:std-b67:2020:full=>2020ncl:std-b67:2020:limited",dither_type="none",resample_filter="spline36",resample_filter_uv="spline36")

ScriptClip("""
mcll=propGetInt("_dovi_max_content_light_level")
mpq=propGetInt("_dovi_max_pq")
subtitle("maxpq = " + string(mpq) + "; maxcll = " + string(mcll))
""")
