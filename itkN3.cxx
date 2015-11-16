#include "itkGDCMImageIO.h"
#include "itkGDCMSeriesFileNames.h"
#include "itkSmoothingRecursiveGaussianImageFilter.h"
#include "itkImageSeriesReader.h"
#include "itkImageFileWriter.h"
#include "itkMetaDataObject.h"

#include "itkN4BiasFieldCorrectionImageFilter.h"
#include "itkOtsuThresholdImageFilter.h"
#include "itkShrinkImageFilter.h"
#include "itkConstantPadImageFilter.h"
#include "gdcmUIDGenerator.h"
#include "itkMetaDataDictionary.h"

#include "metaCommand.h"

// forward declaration
void CopyDictionary (itk::MetaDataDictionary &fromDict, itk::MetaDataDictionary &toDict);

template<typename TFilter>
class CommandIterationUpdate : public itk::Command
{
public:
  typedef CommandIterationUpdate  Self;
  typedef itk::Command            Superclass;
  typedef itk::SmartPointer<Self> Pointer;
  itkNewMacro( Self );

protected:
  CommandIterationUpdate() {}

public:

  virtual void Execute(itk::Object *caller, const itk::EventObject & event) ITK_OVERRIDE
    {
    Execute( (const itk::Object *) caller, event);
    }

  virtual void Execute(const itk::Object * object, const itk::EventObject & event) ITK_OVERRIDE
    {
    const TFilter * filter =
      dynamic_cast< const TFilter * >( object );

    if ( typeid( event ) != typeid( itk::IterationEvent ) )
                                        { return; }
    if ( filter->GetElapsedIterations() == 1 ) {
      std::cout << "Current level = " << filter->GetCurrentLevel() + 1
                << std::endl;
    }
    std::cout << "  Iteration " << filter->GetElapsedIterations()
              << " (of "
              << filter->GetMaximumNumberOfIterations()[ filter->GetCurrentLevel() ]
              << ").  ";
    std::cout << " Current convergence value = "
              << filter->GetCurrentConvergenceMeasurement()
              << " (threshold = " << filter->GetConvergenceThreshold()
              << ")" << std::endl;
    }

};


template<typename TValue>
TValue Convert( std::string optionString )
{
  TValue             value;
  std::istringstream iss( optionString );

  iss >> value;
  return value;
}

template<typename TValue>
std::vector<TValue> ConvertVector( std::string optionString )
{
  std::vector<TValue>    values;
  std::string::size_type crosspos = optionString.find( 'x', 0 );

  if ( crosspos == std::string::npos )
    {
    values.push_back( Convert<TValue>( optionString ) );
    }
  else
    {
    std::string        element = optionString.substr( 0, crosspos );
    TValue             value;
    std::istringstream iss( element );
    iss >> value;
    values.push_back( value );
    while ( crosspos != std::string::npos )
      {
      std::string::size_type crossposfrom = crosspos;
      crosspos = optionString.find( 'x', crossposfrom + 1 );
      if ( crosspos == std::string::npos )
        {
        element = optionString.substr( crossposfrom + 1, optionString.length() );
        }
      else
        {
        element = optionString.substr( crossposfrom + 1, crosspos );
        }
      std::istringstream iss2( element );
      iss2 >> value;
      values.push_back( value );
      }
    }
  return values;
}


int main( int argc, char* argv[] ) {
  
  MetaCommand command;
  command.SetAuthor("Hauke Bartsch");
  command.SetDescription("N3 intensity non-uniformity correction on DICOM images. Read DICOM image series and apply N3 (intensity non-uniformity correction). Exports a new DICOM series.");
  command.AddField("indir", "Directory with input DICOM image series.", MetaCommand::STRING, true);
  command.AddField("outdir", "Directory for output DICOM image series.", MetaCommand::STRING, true);  

  command.SetOption("Mask", "m", false, "Provide a mask image as an additional input. If not supplied Otzu will be used to create a mask.");
  command.AddOptionField("Mask", "mask", MetaCommand::STRING, false);

  command.SetOption("Write", "s", false, "The shrink factor will make the problem easier to handle (sub-sample data). The larger the value the faster.");
  command.AddOptionField("Write", "shrinkFactor", MetaCommand::INT, false, "3");

  command.SetOption("Iterations", "i", false, "Number of iterations \"100x50x50\".");
  command.AddOptionField("Iterations", "iterations", MetaCommand::STRING, false, "100x50x50"); 

  command.SetOption("SeriesName", "n", false, "Select series by series name (if more than one series is present).");
  command.AddOptionField("SeriesName", "seriesname", MetaCommand::STRING, false);
  
  command.SetOption("SaveBiasfield", "b", false, "Save the biasfield as a nifty file in the current directory");
  command.AddOptionField("SaveBiasfield", "biasfieldfilename", MetaCommand::STRING,true);

  command.SetOption("SaveNifty", "n", false, "Save the corrected dataset as a nifty image to the current directory");
  command.AddOptionField("SaveNifty", "niftyfilename", MetaCommand::STRING, true);

  command.SetOption("Verbose", "v", false, "Print more verbose output");


  if ( !command.Parse(argc,argv) ) {
      return 1;
  }

  std::string input  = command.GetValueAsString("indir");
  std::string output = command.GetValueAsString("outdir");
  int shrinkFactor = 3;
  if (command.GetOptionWasSet("Write")) {
     shrinkFactor = command.GetValueAsInt("Write", "shrinkFactor");
  }

  bool verbose = false;
  if (command.GetOptionWasSet("Verbose"))
     verbose = true;

  std::string iterations = "100x50x50";
  if (command.GetOptionWasSet("Iterations")) {
      iterations = command.GetValueAsString("Iterations", "iterations");
  }
  bool saveBiasField = false;  
  bool saveNifty = false;
  bool seriesIdentifierFlag = false;

  if ( command.GetOptionWasSet("SaveBiasfield") )
    saveBiasField = true;
  if ( command.GetOptionWasSet("SaveNifty") )
    saveNifty = true;
  if ( command.GetOptionWasSet("SeriesName") )
    seriesIdentifierFlag = true;
  std::string biasfieldfilename = command.GetValueAsString("SaveBiasfield", "biasfieldfilename");
  std::string niftyfilename = command.GetValueAsString("SaveNifty", "niftyfilename");
  std::string seriesName    = command.GetValueAsString("SeriesName", "seriesname");
  
  
  typedef signed short    PixelType;
  const unsigned int      Dimension = 3;
  
  typedef itk::Image< PixelType, Dimension >         ImageType;
  
  typedef itk::ImageSeriesReader< ImageType >        ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
    
  typedef itk::GDCMImageIO       ImageIOType;
  ImageIOType::Pointer dicomIO = ImageIOType::New();
  
  reader->SetImageIO( dicomIO );
  
  
  typedef itk::GDCMSeriesFileNames NamesGeneratorType;
  NamesGeneratorType::Pointer nameGenerator = NamesGeneratorType::New();
  
  nameGenerator->SetUseSeriesDetails( true );
  nameGenerator->AddSeriesRestriction("0008|0021" );
  
  nameGenerator->SetDirectory( input );
    
  try {
      std::cout << std::endl << "The directory: " << std::endl;
      std::cout << std::endl << input << std::endl << std::endl;
      std::cout << "Contains the following DICOM Series: ";
      std::cout << std::endl << std::endl;
      
      typedef std::vector< std::string >    SeriesIdContainer;
      
      const SeriesIdContainer & seriesUID = nameGenerator->GetSeriesUIDs();
      
      SeriesIdContainer::const_iterator seriesItr = seriesUID.begin();
      SeriesIdContainer::const_iterator seriesEnd = seriesUID.end();
      while( seriesItr != seriesEnd ) {
          std::cout << seriesItr->c_str() << std::endl;
          ++seriesItr;
      }     
      
      std::string seriesIdentifier;
      
      SeriesIdContainer runThese;
      if( seriesIdentifierFlag ) { // If no optional series identifier
          // seriesIdentifier = seriesName;
          runThese.push_back( seriesName );
      } else {
          // Todo: here we select only the first series. We should run
          // N3 on all series.
        
          seriesItr = seriesUID.begin();
          seriesEnd = seriesUID.end();
          // seriesIdentifier = seriesUID.begin()->c_str();
          while (seriesItr != seriesEnd) {
            runThese.push_back( seriesItr->c_str() );
            ++seriesItr;
          }
          // Todo: if we have multiple phases they will all be in the same series. 
          // It does not make sense to handle them here as one big volume, we should
          // look for the slice locations (consecutive) identify the first volume
          // and run N3 on that one. The resulting bias field should be applied to
          // all phases of the series.
      }

      seriesItr = runThese.begin();
      seriesEnd = runThese.end();
      while( seriesItr != seriesEnd) {
        seriesIdentifier = seriesItr->c_str();
        ++seriesItr;
      
        std::cout << std::endl << std::endl;
        std::cout << "Processing series: " << std::endl;
        std::cout << seriesIdentifier << std::endl;
        std::cout << std::endl << std::endl;
      
        typedef std::vector< std::string >   FileNamesContainer;
        FileNamesContainer fileNames;
      
        fileNames = nameGenerator->GetFileNames( seriesIdentifier );
      
        // here we read in all the slices as a single volume for processing
        // if we want to write them back out we have to read them slice by
        // slice and get a copy of the meta data for each slice
        reader->SetFileNames( fileNames );
      
        try {
          reader->Update();
        } catch (itk::ExceptionObject &ex) {
          std::cout << ex << std::endl;
          return EXIT_FAILURE;
        }

        // read the data dictionary      
        ImageType::Pointer inputImage = reader->GetOutput();
        typedef itk::MetaDataDictionary   DictionaryType;
        DictionaryType & dictionary = inputImage->GetMetaDataDictionary();

        // overwrite a value in the dicom dictionary
        //std::string entryId( "0010|0010" );
        //std::string value( "MYNAME" );
        //itk::EncapsulateMetaData<std::string>( dictionary, entryId, value );

        //      
        // replace this with the N4 algorithm instead of the gaussian filter
        //       
        typedef itk::Image<unsigned char, Dimension> MaskImageType;
        MaskImageType::Pointer maskImage = ITK_NULLPTR;

        // if we have a mask on the command line      
        if (command.GetOptionWasSet("Mask")) {
          std::string maskName    = command.GetValueAsString("Mask", "mask");
          typedef itk::ImageFileReader<MaskImageType> MaskReaderType;
          MaskReaderType::Pointer maskreader = MaskReaderType::New();
          maskreader->SetFileName( maskName );
          try {
            maskreader->Update();
            maskImage = maskreader->GetOutput();
            maskImage->DisconnectPipeline();
          } catch( ... ) {
            maskImage = ITK_NULLPTR;
          }
        }

        if ( !maskImage ) {
           //std::cout << "Mask not read.  Creating Otsu mask." << std::endl;
           typedef itk::OtsuThresholdImageFilter<ImageType, MaskImageType> ThresholderType;
           ThresholderType::Pointer otsu = ThresholderType::New();
           otsu->SetInput( inputImage );
           // otsu->SetNumberOfHistogramBins( 200 );
           otsu->SetInsideValue( 0 );
           otsu->SetOutsideValue( 1 );

           otsu->Update();
           maskImage = otsu->GetOutput();
           maskImage->DisconnectPipeline();
        }

        typedef itk::N4BiasFieldCorrectionImageFilter<ImageType, MaskImageType,
                                                    ImageType> CorrecterType;
        CorrecterType::Pointer correcter = CorrecterType::New();
        correcter->SetMaskLabel( 1 );
        correcter->SetSplineOrder( 3 );
        correcter->SetWienerFilterNoise( 0.01 );
        correcter->SetBiasFieldFullWidthAtHalfMaximum( 0.15 );
        correcter->SetConvergenceThreshold( 0.0000001 );

        std::vector<unsigned int> numIters = ConvertVector<unsigned int>( iterations );
        //if( argc > 5 ) {
        //   numIters = ConvertVector<unsigned int>( argv[5] );
        //}
      
        CorrecterType::VariableSizeArrayType  maximumNumberOfIterations( numIters.size() );
        for ( unsigned int d = 0; d < numIters.size(); d++ ) {
           maximumNumberOfIterations[d] = numIters[d];
        }
        correcter->SetMaximumNumberOfIterations( maximumNumberOfIterations );

        CorrecterType::ArrayType numberOfFittingLevels;
        numberOfFittingLevels.Fill( numIters.size() );
        correcter->SetNumberOfFittingLevels( numberOfFittingLevels );


        ImageType::PointType newOrigin = inputImage->GetOrigin();
        CorrecterType::ArrayType numberOfControlPoints;

        float splineDistance = 200;
        //if( argc > 7 )
        //{
        //   splineDistance = atof( argv[7] );
        //}

        itk::SizeValueType lowerBound[Dimension];
        itk::SizeValueType upperBound[Dimension];
        typedef float RealType;

        for( unsigned int d = 0; d < Dimension; d++ ) {
           float domain = static_cast<RealType>( inputImage->GetLargestPossibleRegion().GetSize()[d] - 1 ) * inputImage->GetSpacing()[d];
           unsigned int numberOfSpans = static_cast<unsigned int>(  std::ceil( domain / splineDistance ) );
           unsigned long extraPadding = static_cast<unsigned long>( ( numberOfSpans * splineDistance - domain ) / inputImage->GetSpacing()[d] + 0.5 );
           lowerBound[d] = static_cast<unsigned long>( 0.5 * extraPadding );
           upperBound[d] = extraPadding - lowerBound[d];
           newOrigin[d] -= ( static_cast<RealType>( lowerBound[d] ) * inputImage->GetSpacing()[d] );
           numberOfControlPoints[d] = numberOfSpans + correcter->GetSplineOrder();
        }

        typedef itk::ConstantPadImageFilter<ImageType, ImageType> PadderType;
        PadderType::Pointer padder = PadderType::New();
        padder->SetInput( inputImage );
        padder->SetPadLowerBound( lowerBound );
        padder->SetPadUpperBound( upperBound );
        padder->SetConstant( 0 );
        padder->Update(); 

        inputImage = padder->GetOutput();
        inputImage->DisconnectPipeline();

        typedef itk::ConstantPadImageFilter<MaskImageType, MaskImageType>
        MaskPadderType;
        MaskPadderType::Pointer maskPadder = MaskPadderType::New();
        maskPadder->SetInput( maskImage );
        maskPadder->SetPadLowerBound( lowerBound );
        maskPadder->SetPadUpperBound( upperBound );
        maskPadder->SetConstant( 0 );
        maskPadder->Update(); 

        maskImage = maskPadder->GetOutput();
        maskImage->DisconnectPipeline();
 
        correcter->SetNumberOfControlPoints( numberOfControlPoints );

        // handle the shrink factor
        typedef itk::ShrinkImageFilter<ImageType, ImageType> ShrinkerType;
        ShrinkerType::Pointer shrinker = ShrinkerType::New();
        shrinker->SetInput( inputImage );
        shrinker->SetShrinkFactors( shrinkFactor );

        typedef itk::ShrinkImageFilter<MaskImageType, MaskImageType> MaskShrinkerType;
        MaskShrinkerType::Pointer maskshrinker = MaskShrinkerType::New();
        maskshrinker->SetInput( maskImage );
        maskshrinker->SetShrinkFactors( shrinkFactor );

        //if( argc > 4 ) {
        //  shrinker->SetShrinkFactors( atoi( argv[4] ) );
        //  maskshrinker->SetShrinkFactors( atoi( argv[4] ) );
        //}
        shrinker->Update();
        inputImage = shrinker->GetOutput();
        inputImage->DisconnectPipeline(); 

        maskshrinker->Update();
        maskImage = maskshrinker->GetOutput();
        maskImage->DisconnectPipeline();

        // set the input image and mask image
        correcter->SetInput( inputImage );
        correcter->SetMaskImage( maskImage );

        typedef CommandIterationUpdate<CorrecterType> CommandType;
        CommandType::Pointer observer = CommandType::New();
        correcter->AddObserver( itk::IterationEvent(), observer );

        try {
           correcter->Update();
        } catch( itk::ExceptionObject &excep ) {
           std::cerr << "Exception caught !" << std::endl;
           std::cerr << excep << std::endl;
           return EXIT_FAILURE;
        }
        std::cout << "Final convergence measure: " << correcter->GetCurrentConvergenceMeasurement() << std::endl;

        if (verbose)
          correcter->Print( std::cout, 3 );

        // Test the reconstruction of the log bias field
        typedef itk::Image< float, Dimension >    BiasFieldImageType;
        typedef BiasFieldImageType::Pointer BiasFieldImagePointer;
        typedef ImageType::Pointer  ImagePointer;
        ImagePointer originalInputImage = reader->GetOutput();
        reader->UpdateOutputInformation();
        typedef itk::BSplineControlPointImageFilter <CorrecterType::BiasFieldControlPointLatticeType, CorrecterType::ScalarImageType> BSplinerType;
        BSplinerType::Pointer bspliner = BSplinerType::New();
        bspliner->SetInput( correcter->GetLogBiasFieldControlPointLattice() );
        bspliner->SetSplineOrder( correcter->GetSplineOrder() );
        bspliner->SetSize( originalInputImage->GetLargestPossibleRegion().GetSize() );
        bspliner->SetOrigin( originalInputImage->GetOrigin() );
        bspliner->SetDirection( originalInputImage->GetDirection() );
        bspliner->SetSpacing( originalInputImage->GetSpacing() );
        bspliner->Update();

        BiasFieldImageType::Pointer logField = BiasFieldImageType::New();
        logField->SetOrigin( bspliner->GetOutput()->GetOrigin() );
        logField->SetSpacing( bspliner->GetOutput()->GetSpacing() );
        logField->SetRegions( bspliner->GetOutput()->GetLargestPossibleRegion().GetSize() );
        logField->SetDirection( bspliner->GetOutput()->GetDirection() );
        logField->Allocate();

        itk::ImageRegionIterator<CorrecterType::ScalarImageType> ItB(
                bspliner->GetOutput(),
                bspliner->GetOutput()->GetLargestPossibleRegion() );
        itk::ImageRegionIterator<BiasFieldImageType> ItF( logField, logField->GetLargestPossibleRegion() );
        for ( ItB.GoToBegin(), ItF.GoToBegin(); !ItB.IsAtEnd(); ++ItB, ++ItF )  {
            ItF.Set( ItB.Get()[0] );
        }

        typedef itk::ExpImageFilter<BiasFieldImageType, BiasFieldImageType> ExpFilterType;
        ExpFilterType::Pointer expFilter = ExpFilterType::New();
        expFilter->SetInput( logField );
        expFilter->Update();

        typedef itk::DivideImageFilter<ImageType, BiasFieldImageType, BiasFieldImageType> DividerType;
        DividerType::Pointer divider = DividerType::New();
        divider->SetInput1( reader->GetOutput() );
        divider->SetInput2( expFilter->GetOutput() );
        divider->Update();
      
        // cast the output of divider back to ImageType to get 16bit int again
        typedef itk::Image< ImageType, 3 > OutputImageType;
        typedef itk::CastImageFilter< 
                        BiasFieldImageType,
                        ImageType > CastFilterType;
        CastFilterType::Pointer  caster =  CastFilterType::New();
        caster->SetInput( divider->GetOutput() );
        caster->Update();

        //typedef itk::ImageFileWriter<CorrecterType::BiasFieldControlPointLatticeType> WriterType2;
        // write out the bias field (should be optional)
        if (saveBiasField) {
           typedef itk::ImageFileWriter<BiasFieldImageType> WriterType2;
           WriterType2::Pointer writer2 = WriterType2::New();
           writer2->SetFileName( biasfieldfilename );
           // writer2->SetInput( correcter->GetLogBiasFieldControlPointLattice() );
           writer2->SetInput( expFilter->GetOutput() );
           writer2->Update();
        }

        if (saveNifty) {
           typedef itk::ImageFileWriter< ImageType > WriterType;
           WriterType::Pointer writer = WriterType::New();      
           writer->SetFileName( niftyfilename );
           writer->SetInput( caster->GetOutput() );
     
           //////////////////////////////////////////////  
 
           /* typedef itk::SmoothingRecursiveGaussianImageFilter< ImageType, ImageType >  FilterType;
      
           FilterType::Pointer filter = FilterType::New();
      
           filter->SetInput( reader->GetOutput() );
      
           const double sigma = atof( argv[3] );
           filter->SetSigma( sigma );

      
           typedef itk::ImageFileWriter< ImageType > WriterType;
           WriterType::Pointer writer = WriterType::New();
      
           writer->SetFileName( argv[2] );
           writer->SetInput( filter->GetOutput() );
           //writer->SetImageIO( dicomIO ); */
      
           std::cout  << "Writing the image as " << std::endl << std::endl;
           std::cout  << argv[2] << std::endl << std::endl;      
      
           try  {
             writer->Update();
           } catch (itk::ExceptionObject &ex) {
             std::cout << ex << std::endl;
             return EXIT_FAILURE;
           }
        }
      
        // now save as DICOM
        gdcm::UIDGenerator suid;
        std::string newSeriesUID = suid.Generate();
        gdcm::UIDGenerator fuid;
        std::string frameOfReferenceUID = fuid.Generate();
 
        // create the output directory for the DICOM data     
        itksys::SystemTools::MakeDirectory( output );
      
        // now read in each input file in a loop, copy the result data over and write out as DICOM
        for (int i = 0; i < fileNames.size(); i++) {
           //std::cout << "use slice: " << fileNames[i] << " as template for output" << std::endl;
        
           // this is 2D work
           typedef signed short InputPixelType;
           const unsigned int   Dimension = 2;
           typedef itk::Image< InputPixelType, Dimension > InputImageType;

           typedef itk::ImageFileReader< InputImageType > ReaderType;
           ReaderType::Pointer reader = ReaderType::New();
           reader->SetFileName( fileNames[i] );

           typedef itk::GDCMImageIO           ImageIOType;
           ImageIOType::Pointer gdcmImageIO = ImageIOType::New();
           reader->SetImageIO( gdcmImageIO );
       
           try {
              reader->Update();
           } catch (itk::ExceptionObject & e) {
              std::cerr << "exception in file reader " << std::endl;
              std::cerr << e.GetDescription() << std::endl;
              std::cerr << e.GetLocation() << std::endl;
              return EXIT_FAILURE;
           }
           // ReaderType::DictionaryRawPointer inputDict = (*(reader->GetMetaDataDictionaryArray()))[0];        
        
           InputImageType::Pointer inputImage = reader->GetOutput();
           InputImageType::RegionType region;
           region = inputImage->GetBufferedRegion();
           InputImageType::SizeType size  = region.GetSize();
           // std::cout << "size is: " << size[0] << " " << size[1] << std::endl;

           InputImageType::PixelContainer* container;
           container = inputImage->GetPixelContainer();
           container->SetContainerManageMemory( false );
           unsigned int bla = sizeof(InputImageType::PixelType);
           InputImageType::PixelType* buffer2 = container->GetBufferPointer();
         
           ImageType::Pointer nImage = caster->GetOutput();
           InputImageType::PixelContainer* container2;
           container2 = nImage->GetPixelContainer();
           InputImageType::PixelType* buffer3 = container2->GetBufferPointer();
         
           memcpy(buffer2, &(buffer3[i*size[0]*size[1]]), size[0]*size[1]*bla);
         
           typedef itk::MetaDataDictionary DictionaryType;
           DictionaryType & dictionary = inputImage->GetMetaDataDictionary();
        
           std::string studyUID;
           std::string sopClassUID;
           itk::ExposeMetaData<std::string>(dictionary, "0020|000d", studyUID);
           itk::ExposeMetaData<std::string>(dictionary, "0008|0016", sopClassUID);
           gdcmImageIO->KeepOriginalUIDOn();


           gdcm::UIDGenerator sopuid;
           std::string sopInstanceUID = sopuid.Generate();
                
           //std::string entryId( "0008|103e" );
           //std::string value( "Intensity Corrected" );
           //itk::EncapsulateMetaData<std::string>( dictionary, entryId, value );
 
           //DictionaryType *dict = new DictionaryType();
 
           // Copy the dictionary from the first slice
           //CopyDictionary (dictionary, *dict);
  
           // Set the UID's for the study, series, SOP  and frame of reference
  
           itk::EncapsulateMetaData<std::string>(dictionary,"0020|000d", studyUID);
           itk::EncapsulateMetaData<std::string>(dictionary,"0020|000e", newSeriesUID);
           itk::EncapsulateMetaData<std::string>(dictionary,"0020|0052", frameOfReferenceUID);

           // these keys don't exist - results in error
           //itk::EncapsulateMetaData<std::string>(dictionary,"0020|0052", "0"); // Intercept
           //itk::EncapsulateMetaData<std::string>(dictionary,"0020|0053", "1"); // Slope

           std::string oldSeriesDesc;
           itk::ExposeMetaData<std::string>(dictionary, "0008|103e", oldSeriesDesc);
 
           std::ostringstream value;
           value.str("");
           value << oldSeriesDesc  << " (intensity corrected)";
           // This is a long string and there is a 64 character limit in the
           // standard
           unsigned lengthDesc = value.str().length();
  
           std::string seriesDesc( value.str(), 0,
                              lengthDesc > 64 ? 64
                              : lengthDesc);
           itk::EncapsulateMetaData<std::string>(dictionary,"0008|103e", seriesDesc);
     
           // copy the values for this slice over
           //CopyDictionary (*dict, dictionary);
        
           // write out the result as a DICOM again
           typedef itk::ImageFileWriter< InputImageType >  Writer1Type;
           Writer1Type::Pointer writer1 = Writer1Type::New();
         
           writer1->SetInput( inputImage );
           std::ostringstream o;
           o << output << "/dicom" << i << ".dcm";
           writer1->SetFileName( o.str() );
           writer1->SetImageIO( gdcmImageIO );
           writer1->Update();
          //std::cout << "done with writing the image...";
        
        }
      } // loop over series 
      
  } catch (itk::ExceptionObject &ex) {
      std::cout << ex << std::endl;
      return EXIT_FAILURE;
  }
    
  return EXIT_SUCCESS;
}


 
void CopyDictionary (itk::MetaDataDictionary &fromDict, itk::MetaDataDictionary &toDict)
{
  typedef itk::MetaDataDictionary DictionaryType;
 
  DictionaryType::ConstIterator itr = fromDict.Begin();
  DictionaryType::ConstIterator end = fromDict.End();
  typedef itk::MetaDataObject< std::string > MetaDataStringType;
 
  while( itr != end )
    {
    itk::MetaDataObjectBase::Pointer  entry = itr->second;
 
    MetaDataStringType::Pointer entryvalue =
      dynamic_cast<MetaDataStringType *>( entry.GetPointer() ) ;
    if( entryvalue )
      {
      std::string tagkey   = itr->first;
      std::string tagvalue = entryvalue->GetMetaDataObjectValue();
      itk::EncapsulateMetaData<std::string>(toDict, tagkey, tagvalue);
      }
    ++itr;
    }
}