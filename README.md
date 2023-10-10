# DoViBaker
Bake the DoVi into your clip

This avisynth plugin reads the Base Layer, Enhancement Layer and RPU data from a profile 7 DolbyVision stream to create a clip with the DolbyVision data baked in.

1. Get dovi_tool
  - https://github.com/quietvoid/dovi_tool/releases/tag/2.1.0
2. Extract the Base and Enhancement Layers separately from the initial profile 7 stream
3. Extract the RPU data from the Enhancement Layer using dovi_tool
4. Write a Avisynth script like the example below

example.avs:
```
bl=DGSource("blclip.dgi")
el=DGSource("elclip.dgi")
DoViBaker(bl,el,rpu="RPU.bin")
```

This plugin uses the metadata from the RPU file to compose the DolbyVision HDR picture out of the Base Layer (BL) and Enhancement Layer (EL). Display Management (DM) metadata will not be processed per default. It is however possible to use level 1 maximal pixel brightness data from DM by providing a collection of LUTs and limits of validity measured in nits of max-content-light-level. These will then be processed internally. (The LUT processing implentation is based on: https://github.com/sekrit-twc/timecube).
```
bl=DGSource("blclip.dgi")
el=DGSource("elclip.dgi")
DoViBaker(bl,el,rpu="RPU.bin",cubes="lut_1000.cube;lut_2000.cube;lut_4000.cube",mclls="1010;2020",cubes_basepath="C:\")
```
This example will use the file lut_1000.cube for frames where the max-content-light-level is below or equal to 1010 nits, the file lut_2000.cube for above 1010 but below or equal 2020 nits and lut_4000.cube for all frames above 2020 nits. All cube files must be available in the path given to cubes_basepath, in this example it would be "C:\\".

You can get the current tonemapping value of max-content-light-level by reading the frame property "\_dovi_max_content_light_level":
```
ScriptClip("""
mcll=propGetInt("_dovi_max_content_light_level")
subtitle("maxcll = " + string(mcll))
""")
```

## Trims
Also it is possible to apply the trims available in the stream. Select which trim to apply using the `trimPq` argument and set `targetMaxNits` and `targetMinNits` as necessary. Be warned however, only the typical CM v2.9 processing is implemented thus far, and most streams have not very optimized parameters, producing suboptimal results. Thus this feature is experimental only!

## Frame Properties
The following frame properties will be set:
- `_Matrix` set to 0, representing that the output is RGB
- `_ColorRange` set to 0 in case of full range and 1 in case of limited range
- `_SceneChangePrev` set to 1 at frames where the stream indicates a scene change
- `_dovi_max_pq` the max_pq of the current scene read from the stream
- `_dovi_max_content_light_level` the equivalend value of maximal nits of the current scene

# DoViAnalyzer
This application analyzes the RPU.bin file in order to show information relevant to deciding whether it is worth to use DoViBaker or if this can be skipped completly and the Base Layer can be used directly.

```
usage: DoViAnalyzer.exe <path_to_rpu.bin_file> <optional_scenecutfile.txt>
```
The output will show the following attributes:
1. clip length
2. overall max-content-light-level
3. unusual color matrices defined by the RPU file
4. mapping non-identity introduced by the RPU file
5. enabled processing of the Enhancement Layer
6. available trims, with the trimPq of each trim being shown in brackets
7. optional warning in case that the RGB output is of limited range

Pay attention to 3-5 since these will indicate if the look of the clip will be different when DolbyVision is taken into account compared to just playing the Base Layer clip. If you are using LUTs, pay attention to 7) since in this case the LUTs provided will need to be different.

Additionally it is possible to generate a scenecutfile based on the information from the RPU file. This might be given to the encoder to improve the scene detection (using the parameter --qpfile for x265).

# Remarks concerning compilation
I had some issues linking against Timecube. I was constantly getting the following error:
```
fatal error C1083: Cannot open compiler generated file: 'x64\Release\timecube.asm': No such file or directory
```

It turnes out that this is related to the following setting: "Properies" (of the timecube project) -> "C++" -> "Output Files" -> "Assembler Output".
Setting this to "No Listing" resolves the issue.

