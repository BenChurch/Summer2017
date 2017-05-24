#include "PlusConfigure.h"
#include "plusTrackedFrame.h"
#include "vtkImageCast.h"
#include "vtkImageData.h"
#include "vtkPlusMetaImageSequenceIO.h"
#include "vtkSmartPointer.h"
#include "vtkPlusTrackedFrameList.h"
#include "vtkUsBoneFilter.h"

#include "vtksys/CommandLineArguments.hxx"

int main(int argc, char **argv)
{
  bool printHelp = false;
  vtksys::CommandLineArguments args;
  
  std::string inputFileName;
  std::string outputFileName;
  std::string configFileName;
  std::string linesImageFileName;
  std::string intermediateImageFileName;
  std::string processedLinesImageFileName;
  int verboseLevel=vtkPlusLogger::LOG_LEVEL_UNDEFINED;

  args.Initialize(argc, argv);
  args.AddArgument("--help", vtksys::CommandLineArguments::NO_ARGUMENT, &printHelp, "Print this help");
  args.AddArgument("--input-seq-file", vtksys::CommandLineArguments::EQUAL_ARGUMENT, &inputFileName, "The filename for the input ultrasound sequence to process.");
  args.AddArgument("--config-file", vtksys::CommandLineArguments::EQUAL_ARGUMENT, &configFileName, "The filename for input config file.");
  args.AddArgument("--output-seq-file", vtksys::CommandLineArguments::EQUAL_ARGUMENT, &outputFileName, "The filename to write the processed sequence to.");
  args.AddArgument("--lines-image-file", vtksys::CommandLineArguments::EQUAL_ARGUMENT, &linesImageFileName, "Optional output files for subsampled input image");
  args.AddArgument("--intermediate-image-file", vtksys::CommandLineArguments::EQUAL_ARGUMENT, &intermediateImageFileName, "Optional output file");
  args.AddArgument("--processedlines-image-file", vtksys::CommandLineArguments::EQUAL_ARGUMENT, &processedLinesImageFileName, "Optional output files for processed subsampled image");
  args.AddArgument("--verbose", vtksys::CommandLineArguments::EQUAL_ARGUMENT, &verboseLevel, "Verbose level (1=error only, 2=warning, 3=info, 4=debug, 5=trace)");

  if (!args.Parse())
  {
    std::cerr << "Error parsing arguments" << std::endl;
    std::cout << "Help: " << args.GetHelp() << std::endl;
    exit(EXIT_FAILURE);
  }

  if (printHelp)
  {
    std::cout << args.GetHelp() << std::endl;
    return EXIT_SUCCESS;
  }

  vtkPlusLogger::Instance()->SetLogLevel(verboseLevel);

  if (inputFileName.empty())
  {
    std::cerr << "--input-seq-file not found!" << std::endl;
    return EXIT_FAILURE;
  }
  
  if (configFileName.empty())
  {
    std::cerr << "--config-file not found!" << std::endl;
    return EXIT_FAILURE;
  }
  
  if (outputFileName.empty())
  {
    std::cerr << "--output-seq-file not found!" << std::endl;
    return EXIT_FAILURE;
  }

  // Read config file

  LOG_DEBUG("Reading config file...")
  vtkSmartPointer<vtkXMLDataElement> configRootElement = vtkSmartPointer<vtkXMLDataElement>::New();
  if (PlusXmlUtils::ReadDeviceSetConfigurationFromFile(configRootElement, configFileName.c_str())==PLUS_FAIL)
  {  
    LOG_ERROR("Unable to read configuration from file " << configFileName.c_str()); 
    return EXIT_FAILURE;
  }
  LOG_DEBUG("Reading config file finished.");

  vtkXMLDataElement* boneFilterConfig = configRootElement->FindNestedElementWithName("UsBoneFilter");
  if (boneFilterConfig == NULL)
  {
    LOG_ERROR("Cannot find UsBoneFilter element in XML tree!");
    return PLUS_FAIL;
  }

  // Read the input sequence.

  vtkSmartPointer<vtkPlusTrackedFrameList> trackedFrameList = vtkSmartPointer<vtkPlusTrackedFrameList>::New();
  trackedFrameList->ReadFromSequenceMetafile(inputFileName.c_str());

  vtkSmartPointer<vtkImageCast> castToUchar = vtkSmartPointer<vtkImageCast>::New();
  castToUchar->SetOutputScalarTypeToUnsignedChar();


  int numberOfFrames = trackedFrameList->GetNumberOfTrackedFrames();
  std::cout << "Number of frames in input: " << numberOfFrames << std::endl;

  // Bone filter.
  
  vtkSmartPointer<vtkUsBoneFilter> boneFilter = vtkSmartPointer<vtkUsBoneFilter>::New();
  
  if (! linesImageFileName.empty())
  {
    boneFilter->SetLinesImageFileName( linesImageFileName );
  }

  if (!intermediateImageFileName.empty())
  {
    boneFilter->SetIntermediateImageFileName(intermediateImageFileName);
  }

  if (! processedLinesImageFileName.empty())
  {
    boneFilter->SetProcessedLinesImageFileName( processedLinesImageFileName );
  }

  boneFilter->SetInputFrames( trackedFrameList );
  boneFilter->ReadConfiguration(boneFilterConfig);
  
  PlusStatus filterStatus = boneFilter->Update();
  if ( filterStatus != PlusStatus::PLUS_SUCCESS )
  {
	  std::cout << "Failed processing frames" << std::endl;
    return EXIT_FAILURE;
  }

  std::cout << "Writing output to file." << std::endl;
  vtkPlusLogger::Instance()->SetLogLevel(3);
  boneFilter->GetOutputFrames()->SaveToSequenceMetafile(outputFileName.c_str());
  
  return EXIT_SUCCESS;
}