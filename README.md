# DoViBaker
Bake-in the DoVi into your clip

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

DoViBaker(bl,el,"RPU.bin")
```

Tonemapping needs to be taken care of externally, which can be done using the exported frame property "_dovi_max_content_light_level":
```
ScriptClip("""
mcll=propGetInt("_dovi_max_content_light_level")
nits=mcll <= 1000 ? "1000" : mcll <= 1400 ? "1400" : mcll <= 2000 ? "2000" : mcll <= 2800 ? "2800" : "4000"
Cube("Z:\lut_"+nits+".cube",fullrange=true)
subtitle("maxcll = " + string(mcll))
""")
```
