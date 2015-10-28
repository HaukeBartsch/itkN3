# itkN3
Using ITK's N3 implementation for reading DICOM images, applying non-uniform intensity normalization and writing out a DICOM dataset again

This program requires the ITK toolkit (http://www.itk.org/).

How to compile
--------------

```
cmake .
./itkN3
Option indir is required but not defined
Option outdir is required but not defined
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
