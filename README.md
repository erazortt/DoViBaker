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

Currently the speed is pretty bad. Even though the plugin does support the typical quarter size el clip, it is suggested to equalize the sizes of the el clip and the bl clip prior to handing them to the plugin to improve the performance. Either by downsizeing the bl clip to the el clip size or upsizing the el clip size to the bl clip size.
