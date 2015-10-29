# itkN3
Using ITK's N3 implementation for reading DICOM images, applying non-uniform intensity normalization and writing out a DICOM dataset again

This program requires the ITK toolkit (http://www.itk.org/) to build. The binary files are statically linked against itk and therefore should not require a separate ITK install on your machine.

   https://github.com/HaukeBartsch/itkN3/raw/master/binary/Linux/itkN3

   https://github.com/HaukeBartsch/itkN3/raw/master/binary/MacOS/itkN3
  

```
 Command tags:
    [ -s [ shrinkFactor ] ]
         = The shrink factor will make the problem easier to handle (sub-sample data). The larger the value the faster.
         With: shrinkFactor (Default = 3)
    [ -i [ iterations ] ]
         = Number of iterations "100x50x50".
         With: iterations (Default = 100x50x50)
    [ -n [ seriesname ] ]
         = Select series by series name (if more than one series is present).
    [ -b < biasfieldfilename > ]
         = Save the biasfield as a nifty file in the current directory
    [ -n < niftyfilename > ]
         = Save the corrected dataset as a nifty image to the current directory
    [ -v ]
         = Print more verbose output
 Command fields:
    < indir >
    = Directory with input DICOM image series.
    < outdir >
    = Directory for output DICOM image series.											  
```

If you have DICOM images in a directory called 'data' you can start the program by running:
```
  ./itkN3 data output
```
The program will parse the data folder, read the first series and run the N3 algorithm. Output DICOM images will appear in the output folder.

You can speed up processing by specifying less iterations and a larger shinkFactor:
```
  ./itkN3 -s 4 -i 50x25x25 data output
```
How to compile
--------------

Compile a static version of the ITK library (on Linux):
```
mkdir ~/InsightToolkit-4.8.1/bin
cd ~/InsightToolkit-4.8.1/bin
cmake -DBUILD_SHARED_LIBS=OFF -DCMAKE_EXE_LINKER_FLAGS="-static" -DCMAKE_FIND_LIBRARY_SUFFIXES=".a" ..
cd ~/
git clone https://github.com/HaukeBartsch/itkN3.git itkN3
cd ~/itkN3
cmake -DITK_DIR=/root/InsightToolkit-4.8.1/bin -DCMAKE_EXE_LINKER_FLAGS="-static" -DCMAKE_FIND_LIBRARY_SUFFIXES=".a" .
file ./itkN3
./itkN3: ELF 64-bit LSB  executable, x86-64, version 1 (GNU/Linux), statically linked, for GNU/Linux 2.6.24, BuildID[sha1]=c59064b228f72b7395060eb8b372b6dbade54980, not stripped
```

On MacOS compile ITK without dynamically linked libraries (edit the CMakeFile after the cmake step). No options are required in that case to produce a static binary which includes the itk libraries.
