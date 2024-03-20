DoViLutGen.exe pq2hlg.cube -s 33 -i 0 -o 0
DoViLutGen.exe pq2sdr2020.cube -s 33 -i 0 -o 1
DoViLutGen.exe pq2sdr709.cube -s 33 -i 0 -o 3
DoViLutGen.exe hlg2sdr2020.cube -s 26 -i 2 -o 1
DoViLutGen.exe hlg2sdr709.cube -s 26 -i 2 -o 3
DoViLutGen.exe bt2020to709.cube -s 26 -i 3 -o 3
