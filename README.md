# DoViBaker
Bake the DoVi into your clip

This avisynth plugin reads the Base Layer, Enhancement Layer and RPU data from a profile 7 DolbyVision stream to create a clip with the DolbyVision data baked in.

## General information
This plugin uses the metadata from and RPU file or from the inside stream itself to compose the DolbyVision HDR picture out of the Base Layer (BL) and Enhancement Layer (EL). Display Management (DM) metadata will not be processed per default. It is however possible to further process the clip using DM data by explicitly enabling [Trims](#trims) or by the means of [DoViTonemap](#dovitonemap) or [DoViCubes](#dovicubes). 

## Feeding the plugin 
To my knowledge there are currently three source libraries that can be used. It is advisable to choose one of them in a speed test on your machine.

### [LSMASHSource](https://github.com/HomeOfAviSynthPlusEvolution/L-SMASH-Works)

example.avs:
```
bl=LWLibavVideoSource("clip.ts", format="YUV420P10", stream_index=0)
el=LWLibavVideoSource("clip.ts", format="YUV420P10", stream_index=1)
DoViBaker(bl,el)
```
### [FFmpegSource](https://codeberg.org/StvG/ffms2)

example.avs:
```
bl=FFVideoSource("clip.ts", threads=1, track=0)
el=FFVideoSource("clip.ts", threads=1, track=1)
DoViBaker(bl,el)
```

### [DGDecNV](https://www.rationalqm.us/dgdecnv/binaries/)

1. Get dovi_tool: https://github.com/quietvoid/dovi_tool/releases/tag/2.1.0
2. Extract the Base and Enhancement Layers separately from the initial profile 7 stream
3. Extract the RPU data from the Enhancement Layer using dovi_tool
4. Write a Avisynth script like the example below

example.avs:
```
bl=DGSource("blclip.dgi")
el=DGSource("elclip.dgi")
DoViBaker(bl,el,rpu="RPU.bin")
```

## Trims
Also it is possible to apply the trims available in the DolbyVision substream. Select which trim to apply using the `trimPq` argument and set `targetMaxNits` and `targetMinNits` as necessary. Be warned however, only the typical CM v2.9 processing is implemented thus far, and most streams don't have very optimized parameters, producing suboptimal results. Thus this feature is experimental only!

Typical trim targets usually available are:
* 100 nits, with a `trimPq` of 2081
* 600 nits, with a `trimPq` of 2851
* 1000 nits, with a `trimPq` of 3079

In comparison to trims and especially for higher brightness targets like 600 nits and above, results might be better using [DoViTonemap](#dovitonemap) with both `masterMaxNits` and `masterMinNits` set to `-1`.

## Frame Properties
The following frame properties will be set:
- `_Matrix` set to 0, representing that the output is RGB
- `_ColorRange` set to 0 in case of full range and 1 in case of limited range
- `_SceneChangePrev` set to 1 for the first frame in a scene
- `_dovi_dynamic_min_pq` the min_pq value of the current scene
- `_dovi_dynamic_max_pq` the max_pq value of the current scene
- `_dovi_dynamic_max_content_light_level` the equivalend value of maximal nits of the current scene
- `_dovi_static_max_pq` the max_pq value of the whole stream
- `_dovi_static_max_content_light_level` the value of maximal nits of the whole stream

You can get the current tonemapping value of max-content-light-level by reading the frame property `_dovi_dynamic_max_content_light_level`:
```
ScriptClip("""
mcll=propGetInt("_dovi_dynamic_max_content_light_level")
subtitle("maxcll = " + string(mcll))
""")
```

# DoViTonemap
This plugin processes the tonemapping of any HDR PQ streams to lower dynamic range targets. The implementation is based on ITU-R BT.2408-7 Annex 5 (was in ITU-R BT.2390 until revision 7), with the addition of an optional luminosity factor which scales the brightness linearily. 

It is not expected to give good results for low brightness targets of 400 nits and below. Also color space is preserved and not converted to narrower gamut. For conversions to HLG or SDR additional processing is required, see [DoViLutGen](#dovilutgen).

The following arguments control the tonemapping function: 
- `masterMaxNits` and `masterMinNits` set the white and black brightness value of the source. The values for the master brightness can be either given explicitly or `masterMaxNits` and `masterMinNits` can both be set to `-1` which will indicate that the actual values are read from the related frame properties `_dovi_dynamic_max_pq` and `_dovi_dynamic_min_pq` which are set by [DoViBaker](#dovibaker) or [DoViStatsFileLoader](#dovistatsfileloader), leading to a dynamic tonemapping. If not given these will default to `-1`.
- `targetMaxNits` and `targetMinNits` set the desired target capabilities. These must be given explicitly.
- `lumScale` changes the total brightness, this can be usefull since many HDR PQ and DV streams are actually too dark, darker then the respective SDR streams. To find the proper `lumScale` factor you might use the script [LumScaleHelper.avs](#lumscalehelperavs). It is also possible to read the luminosity factor from the frame property `_dovi_dynamic_luminosity_scale` by setting `lumSacle` to `-1`. When not given explicitly the default of `1.0` is used.
- `kneeOffset` is a parameter of the tonemapping curve, which governs the size of the region where the tonemapping function is flattened (see figure below). The mathematical validity range is [0.5, 2.0]. In report BT.2408 this value is fixed at 0.5, which leads to low highlight details while favoring max brightness. Here the default value used is `0.75` which should be a better compromise overall, especially when using dynamic tonemapping.
- `normalizeOutput` normalizes the output from the range `[targetMinNits, targetMaxNits]` to the full range. This can be usefull when the output is just an intermediate result which is further processed, since the usage of the full value range decreases rounding errors down the line. Default is `false`.

The following example applies a dynamic tonemapping to a 1000nits target while reading the current max and min brightness values off the frame properties which are set by [DoViBaker](#dovibaker). The luminosity scale is not given thus the default of 1.0 is used.
```
DoViBaker(bl,el)
DoViTonemap(targetMaxNits=1000, targetMinNits=0)
```

If your source is just PQ and doesn't have a DolbyVision substream, there are two options:
- use static tonemapping by explicitly defining `masterMaxNits` and `masterMinNits` to [DoViTonemap](#dovitonemap)
- analyse the source using [StatsFileCreator.avs](#statsfilecreatoravs) and provide the create stats file to [DoViStatsFileLoader](#dovistatsfileloader) for a dynamic tonemapping with [DoViTonemap](#dovitonemap)

Shown below is the functional form of the tonemapping curve with the following parameters: masterMaxNits=10000, targetMaxNits=1000, masterMinNits=0, targetMinNits=0.1, lumscale=1.
![Tonemapping function](EETF.png "Tonemapping function")

## Frame Properties
The following frame properties will be consumed (if the related arguments `masterMaxNits`, `masterMinNits` and `lumScale` are set to `-1`):
- `_dovi_dynamic_max_pq` the max_pq value of the current scene
- `_dovi_dynamic_min_pq` the min_pq value of the current scene
- `_dovi_dynamic_luminosity_scale` the luminosity scaling factor of the current scene
- `_ColorRange` both limited and full range RGB inputs are supported

The following frame properties will be set:
- `_ColorRange` set to 0, since the output is always full range RGB independently of the input

# DoViLutGen
This application generates LUTs for conversions from BT.2100 PQ to BT.2100 HLG or to BT.2020 SDR. The PQ to HLG conversion is based on BT.2408-7 in conjunction with BT.2100-2. The LUTs will only process input values up to 1000 nits and will clip anything above that. If the PQ source has brightness levels above that, use [DoViTonemap](#dovitonemap) to tonemap the PQ stream to 1000 nits.

The generated SDR LUTs provide no colorspace conversion, and create a BT.2020 output. For conversions to BT.709 an additional color conversion is necessary. This can be done using http://avisynth.nl/index.php/Avsresize. 

```
usage: DoViLutGen.exe <output_file> <lut_size> (<normalized_input>) (<sdr>) (<sdr_gain>) (<sdr_compression>)
```

The meaning of the expected arguments:
- `output file` the name of the to-be-generated LUT file
- `lut size` generally a bigger LUT, is a better LUT. A good size is `65`.
- `normalized input` if this optional argument is set to `1`, the generated LUT will expect that the input PQ was re-normalized to 1000 nits max brightness. LUTs for re-normalized inputs can be of smaller size than normal LUTs while still providing better quality. A good size for such a LUT is `50`. When not given, this will default to `0`.
- `sdr` if this optional argument is set to `1`, the generatewd LUT will convert to BT.2020 SDR. Default is `0`, with the generated LUT converting to BT.2100 HLG.
- `sdr_gain` this optional argument adjusts the SDR mapping function, by setting the amount of gain of the bright midtones. Value range is [0.0, 1.0], default is `0.0`.
- `sdr_compression` this optional argument adjusts the SDR mapping function, by setting the amount of compression of the very bright highlights. Value range is [0.0, 1.0], default is `1.0`.

## Workflow for conversion to HLG
Generate the LUT by the following command:
```
DoViLutGen.exe pq2hlg_normalizedInput.cube 50 1
```

Create the following avisyth script:
```
DoViBaker(bl,el)
DoViTonemap(targetMaxNits=1000, targetMinNits=0, normalizeOutput=true)
AVSCube("pq2sdr_normalizedInput.cube")
z_ConvertFormat(pixel_type="YUV420P16",colorspace_op="rgb:std-b67:2020:full=>2020ncl:std-b67:2020:limited",chromaloc_op="center=>left")
```
## Workflow for conversion to SDR
Generate the LUT by the following command:
```
DoViLutGen.exe pq2sdr_normalizedInput.cube 50 1 1
```

Create the following avisyth script:
```
DoViBaker(bl,el)
DoViTonemap(targetMaxNits=1000, targetMinNits=0, normalizeOutput=true)
AVSCube("pq2sdr_normalizedInput.cube")
z_ConvertFormat(pixel_type="YUV420P8",colorspace_op="rgb:709:2020:full=>709:709:709:limited",chromaloc_op="center=>left")
```

# DoViCubes
This plugin provides LUT processing capabilites based on the frame property `_dovi_dynamic_max_content_light_level` set by either [DoViBaker](#dovibaker) or [DoViStatsFileReader](#dovistatsfilereader). Different LUTs are applied based adjustable thresholds. This is done by providing a collection of LUTs and limits of validity measured in nits of max-content-light-level. (The LUT processing implentation is based on: https://github.com/sekrit-twc/timecube).
```
DoViBaker(bl,el)
DoViCubes(cubes="lut_1000.cube;lut_2000.cube;lut_4000.cube",mclls="1010;2020",cubes_basepath="C:\")
```
This example will use the file lut_1000.cube for frames where the max-content-light-level is below or equal to 1010 nits, the file lut_2000.cube for above 1010 but below or equal 2020 nits and lut_4000.cube for all frames above 2020 nits. All cube files must be available in the path given to cubes_basepath, in this example it would be "C:\\".

## Frame Properties
The following frame properties will be consumed:
- `_dovi_dynamic_max_content_light_level` the maximal nits value of the current scene

# DoViStatsFileLoader
This plugin reads the stats file generated by the avisynth script [StatsFileCreator.avs](statsfilecreatoravs). It can be used for sources which do not have any DolbyVision substream, but where a processing by [DoViCubes](#dovicubes) or [DoViTonemap](#dovitonemap) is still desired.
The format of each line of the stats file needs to be, last entry is optional:  
`<frame_number> <decision_if_frame_is_last_in_scene> <frame_max_pq> <frame_min_pq> <frame_lum_scale>`

Additionally it is also possible to provide another scene cut file, created by other means than through [StatsFileCreator.avs](statsfilecreatoravs). In this case the scene cuts are going to be taken for that file and the stream statistics from the stats file. The format of each line of the optional alternative scene cut file needs to be:  
`<frame_number_of_first_frame_after_scene_cut>`

In this example the input stats file is read feeding [DoViTonemap](#dovitonemap):
```
DoViStatsFileReader("statsFile.txt")
DoViTonemap(targetMaxNits=1000, targetMinNits=0)
```

## Frame Properties
The following frame properties will be set:
- `_SceneChangePrev` set to 1 for the first frame in a scene
- `_SceneChangeNext` set to 1 for the last frame in a scene
- `_dovi_dynamic_min_pq` the min_pq value of the current scene
- `_dovi_dynamic_max_pq` the max_pq value of the current scene
- `_dovi_dynamic_max_content_light_level` the equivalend value of maximal nits of the current scene
- `_dovi_dynamic_luminosity_scale` the optional luminosity scaling factor of the current scene
- `_dovi_static_max_pq` the max_pq value of the whole stream
- `_dovi_static_max_content_light_level` the value of maximal nits of the whole stream

# StatsFileCreator.avs
This avisynth script scans through the clip and writes the stats file needed for [DoViStatsFileLoader](#dovistatsfileloader). The stats file includes scene cuts and per-frame max brightness values. The scene cut detection algorithm is a rather simple implementation, which is however good enough for most cases.

The format of each line of the stats file created is:  
`<frame_number> <decision_if_frame_is_last_in_scene> <frame_max_pq> <frame_min_pq>`

The format of each line of the optional alternative scene cut file created is:  
`<frame_number_of_first_frame_after_scene_cut>`

# LumScaleHelper.avs
Used to find `lumScale` for [DoViTonemap](#dovitonemap) manually. This is the factor by which to mutiply the brightness of the PQ stream such that its base brightness matches that of the SDR stream. Typical factors can be 1.0 all the way to up 5.0 in very extreme cases. Also this factor might fluctuate from scene to scene. In this case it is advisable to use one best fitting factor thoughout the whole stream in order to maintain the creator's intent. For low brightness targets it might however be necessary have the factor adjusted from scene to scene.

# BetterGrayscale.avsi
Needed by [LumScaleHelper.avs](#lumscalehelperavs) for showing a more correct and better comparable grayscale of PQ and SDR sources.

# DoViAnalyzer
This application analyzes the RPU.bin file in order to show information relevant to deciding whether it is worth to use [DoViBaker](#dovibaker) or if this can be skipped completely and the Base Layer can be used directly.

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

Additionally it is possible to generate a scenecutfile based on the information from the RPU file. This might be used as the optional scene cut file by [DoViStatsFileLoader](#dovistatsfileloader). Or it might be given to the encoder to improve the scene detection (using the parameter --qpfile for x265). In this case add " K" to the end of each line of the file.

# AVSCube
This is a simplified version of http://avisynth.nl/index.php/AVSCube with exactly no adjustability whatsoever. It provides the same quality of image processing as the LUT processing done by [DoViCubes](#dovicubes), while supporting only a single LUT. The desire for this implementation stems from the fact that the original version has non-optimal default settings.

# Remarks concerning compilation
I had some issues linking against Timecube. I was constantly getting the following error:
```
fatal error C1083: Cannot open compiler generated file: 'x64\Release\timecube.asm': No such file or directory
```

It turnes out that this is related to the following setting: "Properties" (of the timecube project) -> "C/C++" -> "Output Files" -> "Assembler Output".
Setting this to "No Listing" resolves the issue.
