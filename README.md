# DoViBaker
Bake the DoVi into your clip

This avisynth plugin reads the Base Layer, Enhancement Layer and RPU data from a profile 7 DolbyVision stream to create a clip with the DolbyVision data baked in.

1. Get dovi_tool and libdovi and save the dovi.dll such that it is reachable for avisynth (like in C:\Windows\System32)
  - https://github.com/quietvoid/dovi_tool/releases/tag/1.5.5
  - https://github.com/quietvoid/dovi_tool/releases/tag/libdovi-1.6.7
2. Extract the Base and Enhancement Layers separately from the initial profile 7 stream
3. Extract the RPU data from the Enhancement Layer using dovi_tool
4. Write a Avisynth script like the example below

example.avs:
```
bl=DGSource("blclip.dgi")
el=DGSource("elclip.dgi")
DoViBaker(bl,el,rpu="RPU.bin")
```

This plugin uses the metadata from the RPU file to compose the DolbyVision HDR picture out of the Base Layer (BL) and Enhancement Layer (EL). Display Management (DM) metadata will not be processed. It is however possible to use level 1 maximal pixel brightness data from DM by providing a collection of LUTs and limits of validity measured in nits of max-content-light-level. These will then be processed internally. (The LUT processing implentation is based on: https://github.com/sekrit-twc/timecube).
```
bl=DGSource("blclip.dgi")
el=DGSource("elclip.dgi")
DoViBaker(bl,el,rpu="RPU.bin",cubes="lut_1000.cube;lut_2000.cube;lut_3000.cube",mclls="1000;2000",cubes_basepath="C:\")
```
This will use the file lut_1000.cube for frames where the max-content-light-level is below or equal to 1000nits, the file lut_2000.cube for above 1000 but below or equal 2000 nits and lut_3000.cube for all frames above 2000nits. All cube files must be available in the path given to cubes_basepath, in this example it would be "C:\\".

You can get the current tonemapping value of max-content-light-level by reading the frame property "\_dovi_max_content_light_level":
```
ScriptClip("""
mcll=propGetInt("_dovi_max_content_light_level")
subtitle("maxcll = " + string(mcll))
""")
```

