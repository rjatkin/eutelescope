// -*- mode: c++; mode: auto-fill; mode: flyspell-prog; -*-
// Author Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
// Version $Id: EUTelPedestalNoiseProcessor.cc,v 1.23 2007-09-26 15:15:52 bulgheroni Exp $
/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

// eutelescope includes ".h" 
#include "EUTelExceptions.h"
#include "EUTelRunHeaderImpl.h"
#include "EUTelEventImpl.h"
#include "EUTelPedestalNoiseProcessor.h"
#include "EUTELESCOPE.h"

// marlin includes ".h"
#include "marlin/Processor.h"
#include "marlin/Exceptions.h"
#include "marlin/Global.h"

#ifdef MARLIN_USE_AIDA
// aida includes <.h>
#include <marlin/AIDAProcessor.h>
#include <AIDA/IHistogramFactory.h>
#include <AIDA/IHistogram1D.h>
#include <AIDA/IHistogram2D.h>
#include <AIDA/IProfile2D.h>
#include <AIDA/ITree.h>
#endif

// lcio includes <.h> 
#include <lcio.h>
#include <UTIL/LCTOOLS.h>
#include <UTIL/CellIDEncoder.h>
#include <UTIL/LCTime.h>
#include <IMPL/LCEventImpl.h>
#include <IMPL/TrackerRawDataImpl.h>
#include <IMPL/TrackerDataImpl.h> 

// system includes <>
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <memory>

using namespace std;
using namespace lcio;
using namespace marlin;
using namespace eutelescope;


// definition of static members mainly used to name histograms
#ifdef MARLIN_USE_AIDA
std::string EUTelPedestalNoiseProcessor::_pedeDistHistoName   = "PedeDist";
std::string EUTelPedestalNoiseProcessor::_noiseDistHistoName  = "NoiseDist";
std::string EUTelPedestalNoiseProcessor::_commonModeHistoName = "CommonMode";
std::string EUTelPedestalNoiseProcessor::_pedeMapHistoName    = "PedeMap";
std::string EUTelPedestalNoiseProcessor::_noiseMapHistoName   = "NoiseMap";
std::string EUTelPedestalNoiseProcessor::_statusMapHistoName  = "StatusMap";
std::string EUTelPedestalNoiseProcessor::_tempProfile2DName   = "TempProfile2D";
std::string EUTelPedestalNoiseProcessor::_fireFreqHistoName   = "Firing frequency";
#endif

EUTelPedestalNoiseProcessor::EUTelPedestalNoiseProcessor () :Processor("EUTelPedestalNoiseProcessor") {

  // modify processor description
  _description =
    "EUTelPedestalNoiseProcessor computes the pedestal and noise values of a pixel detector";

  // first of all we need to register the input collection
  registerInputCollection (LCIO::TRACKERRAWDATA, "RawDataCollectionName",
			   "Input raw data collection",
			   _rawDataCollectionName, string ("rawdata"));

  // register compulsory parameters
  registerProcessorParameter ("CalculationAlgorithm",
			      "Select the algorithm for pede/noise calculation",
			      _pedestalAlgo,
			      string (EUTELESCOPE::MEANRMS));
  registerProcessorParameter("BadPixelMaskingAlgorithm",
			     "Select the algorithm for bad pixel masking",
			     _badPixelAlgo,
			     string (EUTELESCOPE::NOISEDISTRIBUTION));
  registerProcessorParameter ("NoOfCMIteration",
			      "Number of common mode suppression iterations",
			      _noOfCMIterations, static_cast < int >(1));
  registerProcessorParameter ("HitRejectionCut",
			      "Threshold for rejection of hit pixel (SNR units)",
			      _hitRejectionCut, static_cast < float >(4));
  registerProcessorParameter ("MaxNoOfRejectedPixels",
			      "Maximum allowed number of rejected pixels per event",
			      _maxNoOfRejectedPixels,
			      static_cast < int >(1000));
  registerProcessorParameter ("BadPixelMaskCut",
			      "Threshold for bad pixel identification",
			      _badPixelMaskCut, static_cast < float >(3.5));
  registerProcessorParameter ("FirstEvent", 
			      "First event for pedestal calculation",
			      _firstEvent, static_cast < int > (0));
  registerProcessorParameter ("LastEvent",
			      "Last event for pedestal calculation",
			      _lastEvent, static_cast < int > (-1));
  registerProcessorParameter ("OutputPedeFile","The filename (w/o .slcio) to store the pedestal file",
			      _outputPedeFileName , string("outputpede")); 

  registerProcessorParameter ("ASCIIOutputSwitch","Set to true if the pedestal should also be saved as ASCII files",
			      _asciiOutputSwitch, static_cast< bool > ( true ) );
    
  // now the optional parameters
  registerOptionalParameter ("PedestalCollectionName",
			     "Pedestal collection name",
			     _pedestalCollectionName, string ("pedestalDB"));
  registerOptionalParameter ("NoiseCollectionName",
			     "Noise collection name", _noiseCollectionName,
			     string ("noiseDB"));
  registerOptionalParameter ("StatusCollectionName",
			     "Status collection name",
			     _statusCollectionName, string ("statusDB"));
  
  registerOptionalParameter ("AdditionalMaskingLoop",
			     "Perform an additional loop for bad pixel masking",
			     _additionalMaskingLoop, static_cast<bool> ( true ) );
  
  _histogramSwitch = true;

}


void EUTelPedestalNoiseProcessor::init () {
  // this method is called only once even when the rewind is active
  // usually a good idea to
  printParameters ();

  // set to zero the run and event counters
  _iRun = 0;
  _iEvt = 0;

  // set the pedestal flag to true and the loop counter to zero
  _doPedestal = true;
  _iLoop = 0;

  if ( _pedestalAlgo == EUTELESCOPE::MEANRMS ) {
    // reset the temporary arrays
    _tempPede.clear ();
    _tempNoise.clear ();
    _tempEntries.clear ();
  }

#ifndef MARLIN_USE_AIDA
  _histogramSwitch = false;
  if ( _pedestalAlgo == EUTELESCOPE::AIDAPROFILE ) {
    streamlog_out ( ERROR0 )  << "The " << EUTELESCOPE::AIDAPROFILE
			      << " algorithm cannot be applied since Marlin is not using AIDA" << endl
			      << " Algorithm changed to " << EUTELESCOPE::MEANRMS << endl;
    _pedestalAlgo = EUTELESCOPE::MEANRMS;
  }
#endif
  

  // reset all the final arrays
  _pedestal.clear();
  _noise.clear();
  _status.clear();

  // clean up things related to CM 
  _noOfSkippedEvent = 0;

}

void EUTelPedestalNoiseProcessor::processRunHeader (LCRunHeader * rdr) {

  _detectorName = rdr->getDetectorName();

  // to make things easier re-cast the input header to the EUTelRunHeaderImpl
  auto_ptr<EUTelRunHeaderImpl> runHeader ( new EUTelRunHeaderImpl(rdr)) ;
  runHeader->addProcessor( type() );

  // increment the run counter
  ++_iRun;

  // let me get from the run header all the available parameter
  _noOfDetector = runHeader->getNoOfDetector();

  // now the four vectors containing the first and the last pixel
  // along both the directions
  _minX = runHeader->getMinX();
  _maxX = runHeader->getMaxX();
  _minY = runHeader->getMinY();
  _maxY = runHeader->getMaxY();

  // make some test on parameters
  if ( ( _pedestalAlgo != EUTELESCOPE::MEANRMS ) &&
       ( _pedestalAlgo != EUTELESCOPE::AIDAPROFILE) 
       ) {
    throw InvalidParameterException(string("_pedestalAlgo cannot be " + _pedestalAlgo));
  }

  if ( ( _badPixelAlgo != EUTELESCOPE::NOISEDISTRIBUTION ) &&
       ( _badPixelAlgo != EUTELESCOPE::ABSOLUTENOISEVALUE )
       ) {
    throw InvalidParameterException(string("_badPixelAlgo cannot be " + _badPixelAlgo));
  }

  
  // the user can decide to limit the pedestal calculation on a
  // sub range of events for many reasons. The most common one is that
  // a run is started with the beam off and some hundreds of events
  // are taken on purpose before switching the beam on. In this way
  // the same file contains both pedestal and data, with pedestal
  // events within a specific event window. 
  //
  // From the Marlin steering file the best way the user has to select
  // this range is using the _firstEvent and _lastEvent parameter of
  // the processor it self. 
  // There is another variable that can influence this behavior and it
  // is the global MaxRecordNumber. Being global, of course it is
  // dominant with respect to the local _lastEvent setting. Once more,
  // if the EORE is found before the _lastEvent than finalize method
  // is called anyway.
  //


  int maxRecordNumber = Global::parameters->getIntVal("MaxRecordNumber");

  streamlog_out ( DEBUG4 )  << "Event range for pedestal calculation is from " << _firstEvent << " to " << _lastEvent 
			    << "\nMaxRecordNumber from the global section is   " << maxRecordNumber << endl;

  // check if the user wants an additional loop to remove the bad pixels
  int additionalLoop = 0; 
  if ( _additionalMaskingLoop ) additionalLoop = 1;

  if ( _lastEvent == -1 ) {
    // the user didn't select an upper limit for the event range, so
    // we don't know on how many events the calculation should be done
    //
    // if the global MaxRecordNumber has been set, so warn the user
    // that the procedure could be wrong due to a too low number of
    // records.    
    if ( maxRecordNumber != 0 ) {
      streamlog_out ( WARNING2 )  << "The MaxRecordNumber in the Global section of the steering file has been set to " 
				  << maxRecordNumber << ".\n" 
				  << "This means that in order to properly perform the pedestal calculation the maximum allowed number of events is "   
				  << maxRecordNumber / ( _noOfCMIterations + 1 + additionalLoop ) << ".\n"
				  << "Let's hope it is correct and try to continue." << endl;
    }
  } else {
    // ok we know on how many events the calculation should be done. 
    // we can compare this number with the maxRecordNumber if
    // different from 0
    if ( maxRecordNumber != 0 ) {
      if ( (_lastEvent - _firstEvent) * ( _noOfCMIterations + 1 + additionalLoop ) > maxRecordNumber ) {
	streamlog_out ( ERROR4 ) << "The pedestal calculation should be done on " << _lastEvent - _firstEvent 
				 << " times " <<  _noOfCMIterations + 1 << " iterations = " 
				 << (_lastEvent - _firstEvent) * ( _noOfCMIterations + 1 + additionalLoop ) << " records.\n" 
				 << "The global variable MarRecordNumber is limited to " << maxRecordNumber << endl;
      throw InvalidParameterException("MaxRecordNumber");
      }
    }
  }


  if ( _iLoop == 0 ) {
    // write the current header to the output condition file
    LCWriter * lcWriter = LCFactory::getInstance()->createLCWriter();

    try {
      lcWriter->open(_outputPedeFileName, LCIO::WRITE_NEW);
    } catch (IOException& e) {
      cerr << e.what() << endl;
      return;
    }
    
    LCRunHeaderImpl    * lcHeader  = new LCRunHeaderImpl;
    EUTelRunHeaderImpl * newHeader = new EUTelRunHeaderImpl (lcHeader);
    
    newHeader->lcRunHeader()->setRunNumber(runHeader->lcRunHeader()->getRunNumber());
    newHeader->lcRunHeader()->setDetectorName(runHeader->lcRunHeader()->getDetectorName());
    newHeader->setHeaderVersion(runHeader->getHeaderVersion());
    newHeader->setDataType(runHeader->getDataType());
    newHeader->setDateTime();
    newHeader->setDAQHWName(runHeader->getDAQHWName());
    newHeader->setDAQHWVersion(runHeader->getDAQHWVersion());
    newHeader->setDAQSWName(runHeader->getDAQSWName());
    newHeader->setDAQSWVersion(runHeader->getDAQSWVersion());  
    newHeader->setNoOfEvent(runHeader->getNoOfEvent());
    newHeader->setNoOfDetector(runHeader->getNoOfDetector());
    newHeader->setMinX(runHeader->getMinY());
    newHeader->setMaxX(runHeader->getMaxX());
    newHeader->setMinY(runHeader->getMinY());
    newHeader->setMaxY(runHeader->getMaxY());
    newHeader->addProcessor(name());
    
    lcWriter->writeRunHeader(lcHeader);
    delete newHeader;
    delete lcHeader;

    lcWriter->close();

    // also book histos
    bookHistos();

  }
}


void EUTelPedestalNoiseProcessor::processEvent (LCEvent * evt) {

  EUTelEventImpl * eutelEvent = static_cast<EUTelEventImpl*> (evt);
  EventType type              = eutelEvent->getEventType();
  

  if ( type == kUNKNOWN ) {
    streamlog_out ( WARNING2 ) << "Event number " << evt->getEventNumber() << " in run " << evt->getRunNumber()
			       << " is of unknown type. Continue considering it as a normal Data Event." << endl;
  }

  if ( _iLoop == 0 ) firstLoop(evt);
  else if ( _additionalMaskingLoop ) {
    if ( _iLoop == _noOfCMIterations + 1 ) {
      additionalMaskingLoop(evt);
    } else {
      otherLoop(evt);
    }
  } else { 
    otherLoop(evt);
  }

}
  


void EUTelPedestalNoiseProcessor::check (LCEvent * evt) {
  // nothing to check here - could be used to fill check plots in reconstruction processor
}


void EUTelPedestalNoiseProcessor::end() {

  if ( _iLoop == _noOfCMIterations + 1 )  streamlog_out ( MESSAGE4 ) << "Successfully finished" << endl;
  else {
    streamlog_out ( ERROR4 ) << "Not all the iterations have been done because of a too MaxRecordNumber.\n"
			     << "Try to increase it in the global section of the steering file." << endl;
    exit(-1);
  }

}

void EUTelPedestalNoiseProcessor::fillHistos() {
  
#ifdef MARLIN_USE_AIDA
  streamlog_out ( MESSAGE2 ) << "Filling final histograms " << endl;

  string tempHistoName;
  for (int iDetector = 0; iDetector < _noOfDetector; iDetector++) {
    int iPixel = 0;
    for (int yPixel = _minY[iDetector]; yPixel <= _maxY[iDetector]; yPixel++) {
      for (int xPixel = _minX[iDetector]; xPixel <= _maxX[iDetector]; xPixel++) {
	if ( _histogramSwitch ) {
	  {
	    stringstream ss;
	    ss << _statusMapHistoName << "-d" << iDetector << "-l" << _iLoop;
	    tempHistoName = ss.str();
	  }
	  if ( AIDA::IHistogram2D * histo = dynamic_cast<AIDA::IHistogram2D*>(_aidaHistoMap[tempHistoName]) ) 
	    histo->fill(static_cast<double>(xPixel), static_cast<double>(yPixel), static_cast<double> (_status[iDetector][iPixel]));
	  else {
	    streamlog_out ( ERROR1 )  << "Not able to retrieve histogram pointer for " << tempHistoName 
				      << ".\nDisabling histogramming from now on " << endl;
	    _histogramSwitch = false;
	  }
	}
	

	if ( _status[iDetector][iPixel] == EUTELESCOPE::GOODPIXEL) {
	  if ( _histogramSwitch ) {
	    {
	      stringstream ss;
	      ss << _pedeDistHistoName << "-d" << iDetector << "-l" << _iLoop;
	      tempHistoName = ss.str();
	    }
	    if ( AIDA::IHistogram1D * histo = dynamic_cast<AIDA::IHistogram1D*> (_aidaHistoMap[tempHistoName]) )
	      histo->fill(_pedestal[iDetector][iPixel]);
	    else {
	      streamlog_out ( ERROR1 )  << "Not able to retrieve histogram pointer for " << tempHistoName 
				      << ".\nDisabling histogramming from now on " << endl;
	      _histogramSwitch = false;
	    }
	  }  
	  
	  if ( streamlog::out.write< streamlog::DEBUG0 > () ) {	 
	    if ( (xPixel == 10) && (yPixel == 10 )) {
	      streamlog::out()  << "Detector " << iDetector << " pedestal " << (_pedestal[iDetector][iPixel]) << endl;
	    } 
	  }

	  if ( _histogramSwitch ) {
	    {
	      stringstream ss;
	      ss << _noiseDistHistoName << "-d" << iDetector << "-l" << _iLoop;
	      tempHistoName = ss.str();
	    }
	    if ( AIDA::IHistogram1D * histo = dynamic_cast<AIDA::IHistogram1D*> (_aidaHistoMap[tempHistoName]) )
	      histo->fill(_noise[iDetector][iPixel]);
	    else {
	      streamlog_out ( ERROR1 )  << "Not able to retrieve histogram pointer for " << tempHistoName 
					<< ".\nDisabling histogramming from now on " << endl;
	      _histogramSwitch = false;
	    }
	  }
	  
	  if ( _histogramSwitch ) {
	    {
	      stringstream ss;
	      ss << _pedeMapHistoName << "-d" << iDetector << "-l" << _iLoop;
	      tempHistoName = ss.str();
	    }
	    if ( AIDA::IHistogram2D * histo = dynamic_cast<AIDA::IHistogram2D*>(_aidaHistoMap[tempHistoName]) )
	      histo-> fill(static_cast<double>(xPixel), static_cast<double>(yPixel), _pedestal[iDetector][iPixel]);
	    else {
	      streamlog_out ( ERROR1 )  << "Not able to retrieve histogram pointer for " << tempHistoName 
					<< ".\nDisabling histogramming from now on " << endl;
	      _histogramSwitch = false;
	    }
	  }
	  
	  if ( _histogramSwitch ) {
	    {
	      stringstream ss;
	      ss << _noiseMapHistoName << "-d" << iDetector << "-l" << _iLoop;
	      tempHistoName = ss.str();
	    } 
	    if ( AIDA::IHistogram2D* histo = dynamic_cast<AIDA::IHistogram2D*>(_aidaHistoMap[tempHistoName]) ) 
	      histo->fill(static_cast<double>(xPixel), static_cast<double>(yPixel), _noise[iDetector][iPixel]);
	    else {
	      streamlog_out ( ERROR1 )  << "Not able to retrieve histogram pointer for " << tempHistoName 
					<< ".\nDisabling histogramming from now on " << endl;
	      _histogramSwitch = false;
	    }
	  }
	}
	++iPixel;
      }
    }
  }
  

#else
  if ( _iEvt == 0 ) streamlog_out ( MESSAGE4 )  << "No histogram produced because Marlin doesn't use AIDA" << endl;
#endif
  
}

void EUTelPedestalNoiseProcessor::maskBadPixel() {


  if ( ( !_additionalMaskingLoop ) ||
       ( _iLoop < _noOfCMIterations + 1 )) {

    vector<double > thresholdVec;
    int    badPixelCounter = 0;
    
    
    if ( _badPixelAlgo == EUTELESCOPE::NOISEDISTRIBUTION ) {
      // to do it we need to know the mean value and the RMS of the noise
      // vector.
      
      for (int iDetector = 0; iDetector < _noOfDetector; iDetector++) {
	
	double sumw  = 0;
	double sumw2 = 0;
	double num   = 0;
	
	// begin a first loop on all pixel to calculate the masking threshold
	for (unsigned int iPixel = 0; iPixel < _status[iDetector].size(); iPixel++) {
	  if ( _status[iDetector][iPixel] == EUTELESCOPE::GOODPIXEL ) {
	    sumw  += _noise[iDetector][iPixel];
	    sumw2 += pow(_noise[iDetector][iPixel],2);
	    ++num;
	  }
	} 
	double meanw      = sumw  / num;
	double meanw2     = sumw2 / num;
	double rms        = sqrt( meanw2 - pow(meanw,2));
	thresholdVec.push_back(meanw + (rms * _badPixelMaskCut) );
	streamlog_out ( DEBUG4 ) << "Mean noise value is " << meanw << " ADC\n"
	  "RMS of noise is " << rms << " ADC\n"
	  "Masking threshold is set to " << thresholdVec[iDetector] << endl;
	
	
      }
    } else if ( _badPixelAlgo == EUTELESCOPE::ABSOLUTENOISEVALUE ) {
      thresholdVec.push_back(_badPixelMaskCut);
    }
    
    const float lowerThreshold = 0.2;
    streamlog_out ( MESSAGE2 ) << "Marking as bad also dead pixels (noise < " << lowerThreshold << " ADC)" << endl;
    
    
    // scan the noise vector again and apply the cut
    for (int iDetector = 0; iDetector < _noOfDetector; iDetector++) {
      for (unsigned int iPixel = 0; iPixel < _status[iDetector].size(); iPixel++) {
	if ( ( 
	      ( _noise[iDetector][iPixel] > thresholdVec[iDetector] ) || 
	      ( _noise[iDetector][iPixel] < lowerThreshold ) 
	       ) && 
	     ( _status[iDetector][iPixel] == EUTELESCOPE::GOODPIXEL ) ) {
	  _status[iDetector][iPixel] = EUTELESCOPE::BADPIXEL;
	  streamlog_out ( DEBUG0 ) <<  "Masking pixel number " << iPixel 
				   << " on detector " << iDetector 
				   << " (" << _noise[iDetector][iPixel] 
				   << " > " << thresholdVec[iDetector] << ")" << endl;
	  ++badPixelCounter;
	}
      }
    } // end loop on detector;
    streamlog_out ( MESSAGE4 )  << "Masked " << badPixelCounter << " bad pixels " << endl;
  } 




  if ( (  _additionalMaskingLoop ) &&
       ( _iLoop == _noOfCMIterations + 1 )) { 
    // now masking relying on the additional loop
    // for the time being this is hardcoded
    float _maxFreq = 0.25;
    int badPixelCounter = 0;
    for ( int iDetector = 0 ; iDetector < _noOfDetector; iDetector++ ) {
      
      for (unsigned int iPixel = 0; iPixel < _status[iDetector].size(); iPixel++) {
#ifdef MARLIN_USE_AIDA
	  if ( _histogramSwitch ) {
	    string tempHistoName;
	    {
	      stringstream ss; 
	      ss << _fireFreqHistoName << "-d" << iDetector << "-l" << _iLoop;
	      tempHistoName = ss.str();
	    }
	    if ( AIDA::IHistogram1D * histo = dynamic_cast<AIDA::IHistogram1D*> ( _aidaHistoMap[ tempHistoName ] ))
	      histo->fill( (static_cast<double> ( _hitCounter[ iDetector ][ iPixel ] )) / _iEvt );
	    
	  }
#endif
	  if ( _hitCounter[ iDetector ][ iPixel ] > _maxFreq * _iEvt ) {
	    _status[ iDetector ][ iPixel ] = EUTELESCOPE::BADPIXEL;
	    ++badPixelCounter;
	  }
      }
      
    } // end loop on detector
    streamlog_out ( MESSAGE4 ) << "Masked " << badPixelCounter << " bad pixels " << endl;
  }
}


void EUTelPedestalNoiseProcessor::firstLoop(LCEvent * event) {

  EUTelEventImpl * evt = static_cast<EUTelEventImpl*> (event);


  // do some checks in order to see if we have to continue or to stop
  // with the processor.
  //
  // 1. we have to go immediately to the finalize if this is a EORE event
  // 2. we have to go to the finalize if the user select an event
  // range for pedestal calculation and the current event number is
  // out of range
  // 3. we have to skip this event if _iEvt is < than the first event
  // selected for pedestal calculation

  
  if ( evt->getEventType() == kEORE ) {
    streamlog_out ( DEBUG4 ) << "EORE found: calling finalizeProcessor()." << endl;
    finalizeProcessor();
  }

  if ( ( _lastEvent != -1 ) && ( _iEvt >= _lastEvent ) ) {
    streamlog_out ( DEBUG4 ) << "Looping limited by _lastEvent: calling finalizeProcessor()." << endl;
    finalizeProcessor();
  }

  if ( _iEvt < _firstEvent ) {
    ++_iEvt;
    throw SkipEventException(this);
  }

  if (_iEvt % 10 == 0) 
    streamlog_out( MESSAGE4 ) << "Processing event " 
			      << setw(6) << setiosflags(ios::right) << event->getEventNumber() << " in run "
			      << setw(6) << setiosflags(ios::right) << setfill('0')  << event->getRunNumber() << setfill(' ')
			      << " (Total = " << setw(10) << _iEvt << ")" << resetiosflags(ios::left) 
			      << " - loop " << _iLoop <<  endl;

  // let me get the rawDataCollection. This is should contain a TrackerRawDataObject
  // for each detector plane in the telescope.
  try { 
    LCCollectionVec *collectionVec = dynamic_cast < LCCollectionVec * >(evt->getCollection (_rawDataCollectionName));

  
    if ( isFirstEvent() ) {
      // the collection contains several TrackerRawData

      for (int iDetector = 0; iDetector < _noOfDetector; iDetector++) {
	// _tempPedestal, _tempNoise, _tempEntries are vector of vector.
	// they have been already cleared in the init() method we are
	// already looping on detectors, so we just need to push back a
	// vector empty for each cycle
	// 
	// _tempPedestal should be initialized with the adcValues, while
	// _tempNoise and _tempEntries must be initialized to zero. Since
	// adcValues is a vector of shorts, we need to copy each
	// elements into _tempPedestal with a suitable re-casting
	
	// get the TrackerRawData object from the collection for this detector
	
	TrackerRawData *trackerRawData = dynamic_cast < TrackerRawData * >(collectionVec->getElementAt (iDetector));
	ShortVec adcValues = trackerRawData->getADCValues ();
	
	if ( _pedestalAlgo == EUTELESCOPE::MEANRMS ) {
	  // in the case of MEANRMS we have to deal with the standard
	  // vectors
	  ShortVec::iterator iter = adcValues.begin();
	  FloatVec tempDoubleVec;
	  while ( iter != adcValues.end() ) {
	    tempDoubleVec.push_back( static_cast< double > (*iter));
	    ++iter;
	  }
	  _tempPede.push_back(tempDoubleVec);
	  
	  // initialize _tempNoise and _tempEntries with all zeros and
	  // ones
	  _tempNoise.push_back(FloatVec(adcValues.size(), 0.));
	  _tempEntries.push_back(IntVec(adcValues.size(), 1));
	  

	} else if ( _pedestalAlgo == EUTELESCOPE::AIDAPROFILE ) {
#ifdef MARLIN_USE_AIDA
	  // in the case of AIDAPROFILE we don't need any vectors since
	  // everything is done by the IProfile2D automatically
	  int iPixel = 0;
	  stringstream ss;
	  ss << _tempProfile2DName << "-d" << iDetector;
	  for (int yPixel = _minY[iDetector]; yPixel <= _maxY[iDetector]; yPixel++) {
	    for (int xPixel = _minX[iDetector]; xPixel <= _maxX[iDetector]; xPixel++) {
	      double temp = static_cast<double> (adcValues[iPixel]);
	      if ( AIDA::IProfile2D * profile = dynamic_cast<AIDA::IProfile2D*> (_aidaHistoMap[ss.str()]) )
		profile ->fill(static_cast<double> (xPixel), static_cast<double> (yPixel), temp);
	      else {
		streamlog_out ( ERROR4 )  << "Irreversible error: " << ss.str() << " is not available. Sorry for quitting." << endl;
		exit(-1);
	      }
	      ++iPixel;
	    }
	  }
#endif
	}
	
	// the status vector can be initialize as well with all
	// GOODPIXEL
	_status.push_back(ShortVec(adcValues.size(), EUTELESCOPE::GOODPIXEL));

	// if the user wants to add an additional loop on events to
	// mask pixels singing too loud, so the corresponding counter
	// vector should be reset
	if ( _additionalMaskingLoop ) _hitCounter.push_back( ShortVec( adcValues.size(), 0) );

      } // end of detector loop
      
      // nothing else to do in the first event
      _isFirstEvent = false;

    } else {     // end of first event
      
      // after the firstEvent all temp vectors and the status one have
      // the correct number of entries for both indeces
      // loop on the detectors
      for (int iDetector = 0; iDetector < _noOfDetector; iDetector++) {
	
	// get the TrackerRawData object from the collection for this plane
	TrackerRawData *trackerRawData = dynamic_cast < TrackerRawData * >(collectionVec->getElementAt (iDetector));
	ShortVec adcValues = trackerRawData->getADCValues ();
	
	
	if ( _pedestalAlgo == EUTELESCOPE::MEANRMS ) {
	  
	  // start looping on all pixels
	  int iPixel = 0;
	  for (int yPixel = _minY[iDetector]; yPixel <= _maxY[iDetector]; yPixel++) {
	    for (int xPixel = _minX[iDetector]; xPixel <= _maxX[iDetector]; xPixel++) {
	      _tempEntries[iDetector][iPixel] =  _tempEntries[iDetector][iPixel] + 1;
	      _tempPede[iDetector][iPixel]    = ((_tempEntries[iDetector][iPixel] - 1) * _tempPede[iDetector][iPixel]
						 + adcValues[iPixel]) / _tempEntries[iDetector][iPixel];
	      _tempNoise[iDetector][iPixel]   = sqrt(((_tempEntries[iDetector][iPixel] - 1) * pow(_tempNoise[iDetector][iPixel],2) 
						      + pow(adcValues[iPixel] - _tempPede[iDetector][iPixel], 2)) / 
						     _tempEntries[iDetector][iPixel]);
	      ++iPixel;
	    } // end loop on xPixel
	  }   // end loop on yPixel	  
	  
	} else if ( _pedestalAlgo == EUTELESCOPE::AIDAPROFILE ) {
	  
#ifdef MARLIN_USE_AIDA	
	  stringstream ss;
	  ss << _tempProfile2DName << "-d" << iDetector;
	  
	  int iPixel = 0;
	  for (int yPixel = _minY[iDetector]; yPixel <= _maxY[iDetector]; yPixel++) {
	    for (int xPixel = _minX[iDetector]; xPixel <= _maxX[iDetector]; xPixel++) {
	      if ( AIDA::IProfile2D* profile = dynamic_cast<AIDA::IProfile2D*> (_aidaHistoMap[ss.str()]) )
		profile->fill(static_cast<double> (xPixel), static_cast<double> (yPixel), static_cast<double> (adcValues[iPixel]));
	      else {
		streamlog_out ( ERROR ) << "Irreversible error: " << ss.str() << " is not available. Sorry for quitting." << endl;
		exit(-1);
	      }
	      ++iPixel;
	    }
	  }
#endif
	}
	
      }     // end loop on detectors
    }
    
    // increment the event counter
    ++_iEvt;
  } catch (DataNotAvailableException& e) {
    streamlog_out ( WARNING2 ) << "No input collection " << _rawDataCollectionName << " is not available in the current event" << endl;
  }

}

void EUTelPedestalNoiseProcessor::otherLoop(LCEvent * event) {

  EUTelEventImpl * evt = static_cast<EUTelEventImpl*> (event);

  // do some checks in order to see if we have to continue or to stop
  // with the processor.
  //
  // 1. we have to go immediately to the finalize if this is a EORE event
  // 2. we have to go to the finalize if the user select an event
  // range for pedestal calculation and the current event number is
  // out of range
  // 3. we have to skip this event if _iEvt is < than the first event
  // selected for pedestal calculation

  
  if ( evt->getEventType() == kEORE ) finalizeProcessor();
  if ( ( _lastEvent != -1 ) && ( _iEvt >= _lastEvent ) ) finalizeProcessor();
  if ( _iEvt < _firstEvent ) {
    ++_iEvt;
    throw SkipEventException(this);
  }
  
  // keep the user updated
  if ( _iEvt % 10 == 0 ) 
    streamlog_out( MESSAGE4 ) << "Processing event " 
			      << setw(6) << setiosflags(ios::right) << event->getEventNumber() << " in run "
			      << setw(6) << setiosflags(ios::right) << setfill('0')  << event->getRunNumber() << setfill(' ')
			      << " (Total = " << setw(10) << _iEvt << ")" << resetiosflags(ios::left) 
			      << " - loop " << _iLoop <<  endl;
  
  // let me get the rawDataCollection. This is should contain a TrackerRawDataObject
  // for each detector plane in the telescope.
  try {
    LCCollectionVec *collectionVec = dynamic_cast < LCCollectionVec * >(evt->getCollection (_rawDataCollectionName));
    
    for (int iDetector = 0; iDetector < _noOfDetector; iDetector++) {
      // get the TrackerRawData object from the collection for this detector
      TrackerRawData *trackerRawData = dynamic_cast < TrackerRawData * >(collectionVec->getElementAt (iDetector));
      ShortVec adcValues = trackerRawData->getADCValues ();
      
      // prepare stuff for common mode calculation: pixelSum is the
      // sum of all good pixel signals. Pixels are identify as good if
      // their status is good and it they are not recognized as hit
      // pixel by the hit rejection threshold. goodPixel is the number
      // of good pixel in this detector used for common mode
      // calculation. commonMode is, indeed, the mean value of the
      // pixel signals pedestal sub'ed. iPixel is a pixel counter
      double pixelSum     = 0.;
      double commonMode   = 0.;
      int    goodPixel    = 0;
      int    skippedPixel = 0;
      int    iPixel       = 0;
      
      // start looping on all pixels for hit rejection
      for (int yPixel = _minY[iDetector]; yPixel <= _maxY[iDetector]; yPixel++) {
	for (int xPixel = _minX[iDetector]; xPixel <= _maxX[iDetector]; xPixel++) {
	  bool isHit  = ( ( adcValues[iPixel] - _pedestal[iDetector][iPixel] ) > _hitRejectionCut * _noise[iDetector][iPixel] );
	  bool isGood = ( _status[iDetector][iPixel] == EUTELESCOPE::GOODPIXEL );
	  if ( !isHit && isGood ) {
	    pixelSum += adcValues[iPixel] - _pedestal[iDetector][iPixel];
	    ++goodPixel;
	  } else if ( isHit ) {
	    ++skippedPixel;
	  }
	  ++iPixel;
	} 
      }
      
      if ( ( skippedPixel < _maxNoOfRejectedPixels ) &&
	   ( goodPixel != 0 ) ) {
	
	commonMode = pixelSum / goodPixel;
#ifdef MARLIN_USE_AIDA      
	stringstream ss;
	ss << _commonModeHistoName << "-d" << iDetector << "-l" << _iLoop;
	(dynamic_cast<AIDA::IHistogram1D*>(_aidaHistoMap[ss.str()]))->fill(commonMode);
#endif
	iPixel = 0;
	for (int yPixel = _minY[iDetector]; yPixel <= _maxY[iDetector]; yPixel++) {
	  for (int xPixel = _minX[iDetector]; xPixel <= _maxX[iDetector]; xPixel++) {
	    if ( _status[iDetector][iPixel] == EUTELESCOPE::GOODPIXEL ) {
	      double pedeCorrected = adcValues[iPixel] - commonMode;
	      if ( _pedestalAlgo == EUTELESCOPE::MEANRMS ) {
		
		_tempEntries[iDetector][iPixel] = _tempEntries[iDetector][iPixel] + 1;
		_tempPede[iDetector][iPixel]    = ((_tempEntries[iDetector][iPixel] - 1) * _tempPede[iDetector][iPixel]
						   + pedeCorrected) / _tempEntries[iDetector][iPixel];
		_tempNoise[iDetector][iPixel]   = sqrt(((_tempEntries[iDetector][iPixel] - 1) * pow(_tempNoise[iDetector][iPixel],2) 
							+ pow(pedeCorrected - _tempPede[iDetector][iPixel], 2)) / 
						       _tempEntries[iDetector][iPixel]);
		
	      } else if ( _pedestalAlgo == EUTELESCOPE::AIDAPROFILE) {
#ifdef MARLIN_USE_AIDA	      
		stringstream ss;
		ss << _tempProfile2DName << "-d" << iDetector;
		(dynamic_cast<AIDA::IProfile2D*> (_aidaHistoMap[ss.str()]))
		  ->fill(static_cast<double> (xPixel), static_cast<double> (yPixel), pedeCorrected);
#endif
	      }
	    }
	    ++iPixel;
	  }
	}
      } else {
	streamlog_out ( WARNING2 ) <<  "Skipping event " << _iEvt << " because of max number of rejected pixels exceeded. (" 
				   << skippedPixel << ")" << endl;
	++_noOfSkippedEvent;
      }
    }	
    ++_iEvt;
  } catch (DataNotAvailableException& e) {
    streamlog_out ( WARNING2 ) << "No input collection " << _rawDataCollectionName << " is not available in the current event" << endl;
  }
}

void EUTelPedestalNoiseProcessor::bookHistos() {

#ifdef MARLIN_USE_AIDA
  // histograms are grouped in loops and detectors
  streamlog_out ( MESSAGE2 ) << "Booking histograms " << endl;


  string tempHistoName;

  // start looping on the number of loops. Remember that we have one
  // loop more than the number of common mode iterations
  for (int iLoop = 0; iLoop < _noOfCMIterations + 1; iLoop++) {
      
    // prepare the name of the current loop directory and add it the
    // the current ITree
    string loopDirName;
    {
      stringstream ss;
      ss << "loop-" << iLoop;
      loopDirName = ss.str();
    }
    AIDAProcessor::tree(this)->mkdir(loopDirName.c_str());
      
    // start looping on detectors
    for (int iDetector = 0; iDetector < _noOfDetector;  iDetector++) {
	
      // prepare the name of the current detector and add it to the
      // current ITree inside the current loop folder
      string detectorDirName;
      {
	stringstream ss;
	ss << "detector-" << iDetector;
	detectorDirName = ss.str();
      }
      string basePath = loopDirName + "/" + detectorDirName + "/";
      AIDAProcessor::tree(this)->mkdir(basePath.c_str());
	
      // book an histogram for the pedestal distribution
      const int    pedeDistHistoNBin   = 100; 
      const double pedeDistHistoMin    = -20.;
      const double pedeDistHistoMax    =  29.;
      {
	stringstream ss;
	ss << _pedeDistHistoName << "-d" << iDetector << "-l" << iLoop;
	tempHistoName = ss.str();
      } 
      AIDA::IHistogram1D * pedeDistHisto = 
	AIDAProcessor::histogramFactory(this)->createHistogram1D( (basePath + tempHistoName).c_str(), 
								  pedeDistHistoNBin, pedeDistHistoMin, pedeDistHistoMax);
      if ( pedeDistHisto ) {
	_aidaHistoMap.insert(make_pair(tempHistoName, pedeDistHisto));
	pedeDistHisto->setTitle("Pedestal distribution");
      } else {
	streamlog_out ( ERROR1 ) << "Problem booking the " << (basePath + tempHistoName) << ".\n"
				 << "Very likely a problem with path name. Switching off histogramming and continue w/o" << endl;
	_histogramSwitch = false;
      }
	
      // book an histogram for the noise distribution
      const int    noiseDistHistoNBin  =  100;
      const double noiseDistHistoMin   =  -5.;
      const double noiseDistHistoMax   =  15.;
      {
	stringstream ss;
	ss << _noiseDistHistoName << "-d" << iDetector << "-l" << iLoop;
	tempHistoName = ss.str();
      }
      AIDA::IHistogram1D * noiseDistHisto =
	AIDAProcessor::histogramFactory(this)->createHistogram1D( (basePath + tempHistoName).c_str(),
								  noiseDistHistoNBin, noiseDistHistoMin, noiseDistHistoMax);
      if ( noiseDistHisto ) {
	_aidaHistoMap.insert(make_pair(tempHistoName, noiseDistHisto));
	noiseDistHisto->setTitle("Noise distribution");
      }	else {
	streamlog_out ( ERROR1 ) << "Problem booking the " << (basePath + tempHistoName) << ".\n"
				 << "Very likely a problem with path name. Switching off histogramming and continue w/o" << endl;
	_histogramSwitch = false;
      }
	
      // book a 1d histo for common mode only if loop >= 1
      if (iLoop >= 1) {
	const int    commonModeHistoNBin = 100;
	const double commonModeHistoMin  =  -2;
	const double commonModeHistoMax  =   2;
	{
	  stringstream ss;
	  ss << _commonModeHistoName << "-d" << iDetector << "-l" << iLoop;
	  tempHistoName = ss.str();
	}
	AIDA::IHistogram1D * commonModeHisto =
	  AIDAProcessor::histogramFactory(this)->createHistogram1D( (basePath + tempHistoName).c_str(),
								    commonModeHistoNBin, commonModeHistoMin, commonModeHistoMax);
	if ( commonModeHisto ) {
	  _aidaHistoMap.insert(make_pair(tempHistoName, commonModeHisto));
	  commonModeHisto->setTitle("Common mode distribution");
	} else {
	  streamlog_out ( ERROR1 ) << "Problem booking the " << (basePath + tempHistoName) << ".\n"
				   << "Very likely a problem with path name. Switching off histogramming and continue w/o" << endl;
	  _histogramSwitch = false;
	}
      }
	
      // book a 2d histogram for pedestal map
      const int    xNoOfPixel = abs( _maxX[iDetector] - _minX[iDetector] + 1);
      const int    yNoOfPixel = abs( _maxY[iDetector] - _minY[iDetector] + 1);
      const double xMin       = _minX[iDetector] - 0.5;
      const double xMax       =  xMin + xNoOfPixel;
      const double yMin       = _minY[iDetector] - 0.5;
      const double yMax       =  yMin + yNoOfPixel;
      {
	stringstream ss;
	ss << _pedeMapHistoName << "-d" << iDetector << "-l" << iLoop;
	tempHistoName = ss.str();
      }
      AIDA::IHistogram2D * pedeMapHisto =
	AIDAProcessor::histogramFactory(this)->createHistogram2D( (basePath + tempHistoName).c_str(),
								  xNoOfPixel, xMin, xMax, yNoOfPixel, yMin, yMax);
      if ( pedeMapHisto ) {
	_aidaHistoMap.insert(make_pair(tempHistoName, pedeMapHisto));
	pedeMapHisto->setTitle("Pedestal map");
      }	else {
	streamlog_out ( ERROR1 ) << "Problem booking the " << (basePath + tempHistoName) << ".\n"
				 << "Very likely a problem with path name. Switching off histogramming and continue w/o" << endl;
	_histogramSwitch = false;
      }
	
      // book a 2d histogram for noise map
      {
	stringstream ss;
	ss << _noiseMapHistoName << "-d" << iDetector << "-l" << iLoop;
	tempHistoName = ss.str();
      }
      AIDA::IHistogram2D * noiseMapHisto =
	AIDAProcessor::histogramFactory(this)->createHistogram2D( (basePath + tempHistoName).c_str(),
								  xNoOfPixel, xMin, xMax, yNoOfPixel, yMin, yMax);	
      if ( noiseMapHisto ) {
	_aidaHistoMap.insert(make_pair(tempHistoName, noiseMapHisto));      
	noiseMapHisto->setTitle("Noise map");
      }	else {
	streamlog_out ( ERROR1 ) << "Problem booking the " << (basePath + tempHistoName) << ".\n"
				 << "Very likely a problem with path name. Switching off histogramming and continue w/o" << endl;
	_histogramSwitch = false;
      }

      // book a 2d histogram for status map
      {
	stringstream ss;
	ss << _statusMapHistoName << "-d" << iDetector << "-l" << iLoop;
	tempHistoName = ss.str();
      }
      AIDA::IHistogram2D * statusMapHisto =
	AIDAProcessor::histogramFactory(this)->createHistogram2D( (basePath + tempHistoName).c_str(),
								  xNoOfPixel, xMin, xMax, yNoOfPixel, yMin, yMax);	
      if ( statusMapHisto ) {
	_aidaHistoMap.insert(make_pair(tempHistoName, statusMapHisto));
	statusMapHisto->setTitle("Status map");
      }	else {
	streamlog_out ( ERROR1 ) << "Problem booking the " << (basePath + tempHistoName) << ".\n"
				 << "Very likely a problem with path name. Switching off histogramming and continue w/o" << endl;
	_histogramSwitch = false;
      }	

      if ( ( _pedestalAlgo == EUTELESCOPE::AIDAPROFILE ) &&
	   ( iLoop == 0 ) ) {
	// we just need to prepare such a 2d profile only in the case
	// we are using the AIDAPROFILE calculation algorithm and we
	// need just one copy of it for each detector.
	{
	  stringstream ss;
	  ss << _tempProfile2DName << "-d" << iDetector;
	  tempHistoName = ss.str();
	}
	AIDA::IProfile2D * tempProfile2D =
	  AIDAProcessor::histogramFactory(this)->createProfile2D( (basePath + tempHistoName).c_str(),
								  xNoOfPixel, xMin, xMax, yNoOfPixel, yMin, yMax,-1000,1000);
	if ( tempProfile2D ) {
	  _aidaHistoMap.insert(make_pair(tempHistoName, tempProfile2D));
	  tempProfile2D->setTitle("Temp profile for pedestal calculation");
	} else {
	  streamlog_out ( ERROR1 ) << "Problem booking the " << (basePath + tempHistoName) << ".\n"
				   << "Very likely a problem with path name. Switching off histogramming and continue w/o" << endl;
	  _histogramSwitch = false;
	  exit(-1);
	}
      }
	
    }
  } // end on iLoop

  if ( _additionalMaskingLoop ) {
    
    int iLoop = _noOfCMIterations + 1 ;
    string loopDirName;
    {
      stringstream ss;
      ss << "loop-" << iLoop;
      loopDirName = ss.str();
    }
    AIDAProcessor::tree(this)->mkdir(loopDirName.c_str());
    
    for ( int iDetector = 0; iDetector < _noOfDetector; iDetector++ ) {

 // prepare the name of the current detector and add it to the
      // current ITree inside the current loop folder
      string detectorDirName;
      {
	stringstream ss;
	ss << "detector-" << iDetector;
	detectorDirName = ss.str();
      }
      string basePath = loopDirName + "/" + detectorDirName + "/";
      AIDAProcessor::tree(this)->mkdir(basePath.c_str());

      // book an histogram for the firing frequency
      const int    fireFreqHistoNBin   = 100; 
      const double fireFreqHistoMin    =  0.0;
      const double fireFreqHistoMax    =  1.0;
      {
	stringstream ss;
	ss << _fireFreqHistoName << "-d" << iDetector << "-l" << iLoop;
	tempHistoName = ss.str();
      } 
      AIDA::IHistogram1D * fireFreqHisto = 
	AIDAProcessor::histogramFactory(this)->createHistogram1D( (basePath + tempHistoName).c_str(), 
								  fireFreqHistoNBin, fireFreqHistoMin, fireFreqHistoMax);
      if ( fireFreqHisto ) {
	_aidaHistoMap.insert(make_pair(tempHistoName, fireFreqHisto));
	fireFreqHisto->setTitle("Firing frequency distribution");
      } else {
	streamlog_out ( ERROR1 ) << "Problem booking the " << (basePath + tempHistoName) << ".\n"
				 << "Very likely a problem with path name. Switching off histogramming and continue w/o" << endl;
	_histogramSwitch = false;
      }
      
    }
  }


#endif // MARLIN_USE_AIDA

}

void EUTelPedestalNoiseProcessor::finalizeProcessor(bool fromMaskingLoop) {
  

  if ( ! fromMaskingLoop ) {

    if ( _pedestalAlgo == EUTELESCOPE::MEANRMS ) {
      
      // the loop on events is over so we need to move temporary vectors
      // to final vectors
      _pedestal = _tempPede;
      _noise    = _tempNoise;
      
      // clear the temporary vectors
      _tempPede.clear();
      _tempNoise.clear();
      _tempEntries.clear();
      
    } else if ( _pedestalAlgo == EUTELESCOPE::AIDAPROFILE ) {
#ifdef MARLIN_USE_AIDA
      _pedestal.clear();
      _noise.clear();
      for (int iDetector = 0; iDetector < _noOfDetector; iDetector++) {
	stringstream ss;
	ss << _tempProfile2DName << "-d" << iDetector;
	FloatVec tempPede;
	FloatVec tempNoise;
	for (int yPixel = _minY[iDetector]; yPixel <= _maxY[iDetector]; yPixel++) {
	  for (int xPixel = _minX[iDetector]; xPixel <= _maxX[iDetector]; xPixel++) {
	    if ( AIDA::IProfile2D * profile = dynamic_cast<AIDA::IProfile2D*> (_aidaHistoMap[ss.str()]) ) {
	      tempPede.push_back((float) profile->binHeight(xPixel,yPixel));
	      // WARNING: the noise part of this algorithm is still not
	      // working probably because of a bug in RAIDA implementation
	      tempNoise.push_back((float) profile->binRms(xPixel,yPixel));
	      //cout << xPixel << " " << yPixel << " " << tempPede.back() << " " << tempNoise.back() << endl;
	    } else {
	      streamlog_out ( ERROR4 )  << "Problem with the AIDA temporary profile.\n"
					<< "Sorry for quitting... " << endl;
	      exit(-1);
	    }
	  }
	}
	_pedestal.push_back(tempPede);
	_noise.push_back(tempNoise);
      }
#endif
    }
    
    // mask the bad pixels here
    maskBadPixel();
    
    // fill in the histograms
    fillHistos();
  
  } else {

    // here refill the status histoMap
    maskBadPixel();
  }


  // increment the loop counter
  ++_iLoop;
  
  // check if we need another loop or we can finish. Remember that we
  // have a total number of loop of _noOfCMIteration + 1 + eventually
  // the additional loop on bad pixel masking
  int additionalLoop = 0;
  if ( _additionalMaskingLoop ) additionalLoop = 1;
  if ( _iLoop == _noOfCMIterations + 1 + additionalLoop ) {
    // ok this was last loop  
    

    // what we need to do now is to check if the user wants an
    // additional loop for more accurate bad pixel masking.
    

    streamlog_out ( MESSAGE4 ) << "Writing the output condition file" << endl;

    LCWriter * lcWriter = LCFactory::getInstance()->createLCWriter();
    
    try {
      lcWriter->open(_outputPedeFileName,LCIO::WRITE_APPEND);
    } catch (IOException& e) {
      cerr << e.what() << endl;
      return;
    }
    
    LCEventImpl * event = new LCEventImpl();
    event->setDetectorName(_detectorName);
    event->setRunNumber(_iRun);
    
    LCTime * now = new LCTime;
    event->setTimeStamp(now->timeStamp());
    delete now;


    LCCollectionVec * pedestalCollection = new LCCollectionVec(LCIO::TRACKERDATA);
    LCCollectionVec * noiseCollection    = new LCCollectionVec(LCIO::TRACKERDATA);
    LCCollectionVec * statusCollection   = new LCCollectionVec(LCIO::TRACKERRAWDATA);
    
    for (int iDetector = 0; iDetector < _noOfDetector; iDetector++) {
      
      TrackerDataImpl    * pedestalMatrix = new TrackerDataImpl;
      TrackerDataImpl    * noiseMatrix    = new TrackerDataImpl;
      TrackerRawDataImpl * statusMatrix   = new TrackerRawDataImpl;
      
      CellIDEncoder<TrackerDataImpl>    idPedestalEncoder(EUTELESCOPE::MATRIXDEFAULTENCODING, pedestalCollection);
      CellIDEncoder<TrackerDataImpl>    idNoiseEncoder(EUTELESCOPE::MATRIXDEFAULTENCODING, noiseCollection);
      CellIDEncoder<TrackerRawDataImpl> idStatusEncoder(EUTELESCOPE::MATRIXDEFAULTENCODING, statusCollection);
      
      idPedestalEncoder["sensorID"] = iDetector;
      idNoiseEncoder["sensorID"]    = iDetector;
      idStatusEncoder["sensorID"]   = iDetector;
      idPedestalEncoder["xMin"]     = _minX[iDetector];
      idNoiseEncoder["xMin"]        = _minX[iDetector];
      idStatusEncoder["xMin"]       = _minX[iDetector];
      idPedestalEncoder["xMax"]     = _maxX[iDetector];
      idNoiseEncoder["xMax"]        = _maxX[iDetector];
      idStatusEncoder["xMax"]       = _maxX[iDetector];
      idPedestalEncoder["yMin"]     = _minY[iDetector];
      idNoiseEncoder["yMin"]        = _minY[iDetector];
      idStatusEncoder["yMin"]       = _minY[iDetector];
      idPedestalEncoder["yMax"]     = _maxY[iDetector];
      idNoiseEncoder["yMax"]        = _maxY[iDetector];
      idStatusEncoder["yMax"]       = _maxY[iDetector];
      idPedestalEncoder.setCellID(pedestalMatrix);
      idNoiseEncoder.setCellID(noiseMatrix);
      idStatusEncoder.setCellID(statusMatrix);
      
      pedestalMatrix->setChargeValues(_pedestal[iDetector]);
      noiseMatrix->setChargeValues(_noise[iDetector]);
      statusMatrix->setADCValues(_status[iDetector]);
      
      pedestalCollection->push_back(pedestalMatrix);
      noiseCollection->push_back(noiseMatrix);
      statusCollection->push_back(statusMatrix);

      if ( _asciiOutputSwitch ) {
	if ( iDetector == 0 ) streamlog_out ( MESSAGE4 ) << "Writing the ASCII pedestal files" << endl;
	stringstream ss;
	ss << _outputPedeFileName << "-b" << iDetector << ".dat";
	ofstream asciiPedeFile(ss.str().c_str());
	asciiPedeFile << "# Pedestal and noise for board number " << iDetector << endl
		      << "# calculated from run " << _outputPedeFileName << endl;

	const int subMatrixWidth = 3;
	const int xPixelWidth    = 4;
	const int yPixelWidth    = 4;
	const int pedeWidth      = 15;
	const int noiseWidth     = 15;
	const int statusWidth    = 3;
	const int precision      = 8;

	int iPixel = 0;
	for (int yPixel = _minY[iDetector]; yPixel <= _maxY[iDetector]; yPixel++) {
	  for (int xPixel = _minX[iDetector]; xPixel <= _maxX[iDetector]; xPixel++) {
	    asciiPedeFile << setiosflags(ios::left) 
			  << setw(subMatrixWidth) << iDetector 
			  << setw(xPixelWidth)    << xPixel
			  << setw(yPixelWidth)    << yPixel
			  << resetiosflags(ios::left) << setiosflags(ios::fixed) << setprecision(precision) 
			  << setw(pedeWidth)      << _pedestal[iDetector][iPixel]
			  << setw(noiseWidth)     << _noise[iDetector][iPixel]
			  << resetiosflags(ios::fixed) 
			  << setw(statusWidth)    << _status[iDetector][iPixel]
			  << endl;
	    ++iPixel;
	  }
	}
	asciiPedeFile.close();
      }
    }

    event->addCollection(pedestalCollection, _pedestalCollectionName);
    event->addCollection(noiseCollection, _noiseCollectionName);
    event->addCollection(statusCollection, _statusCollectionName);
    
    lcWriter->writeEvent(event);
    delete event;

    lcWriter->close();

    throw StopProcessingException(this);
    setReturnValue("IsPedestalFinished", true);
  } else if ( _iLoop < _noOfCMIterations + 1 ) {
    // now we need to loop again
    // so reset the event counter
    _iEvt = 0;

    // prepare everything for the next loop
    if ( _pedestalAlgo == EUTELESCOPE::MEANRMS ) {

      // the collection contains several TrackerRawData
      // move back the _pedestal and _noise to the _temp vector
      _tempPede  = _pedestal;
      _tempNoise = _noise;
      
      for (int iDetector = 0; iDetector < _noOfDetector; iDetector++) {
	_tempEntries.push_back(IntVec( _noise[iDetector].size(), 1));
      }
      
    } else if ( _pedestalAlgo == EUTELESCOPE::AIDAPROFILE ) {
      // in case the AIDAPROFILE algorithm is used, the only thing we
      // need to do is to clean up the previous loop histograms
      // remember to loop over all detectors
#ifdef MARLIN_USE_AIDA
      for (int iDetector = 0; iDetector < _noOfDetector; iDetector++) {
	stringstream ss;
	ss << _tempProfile2DName << "-d" << iDetector;
	if ( AIDA::IProfile2D* profile = dynamic_cast<AIDA::IProfile2D*> (_aidaHistoMap[ss.str()]) ) 
	  profile->reset();
	else {
	  streamlog_out ( ERROR4 ) << "Unable to reset the AIDA temporary profile.\n"
				   << "Sorry for quitting..." << endl;
	  exit(-1);
	}
      }
#endif
    }
    setReturnValue("IsPedestalFinished", false);
    throw RewindDataFilesException(this);
  } else if ( ( _additionalMaskingLoop ) &&
	      ( _iLoop ==  _noOfCMIterations + 1 ) ) {
    // additional loop! 
    // now we need to loop again
    // so reset the event counter
    _iEvt = 0;
    setReturnValue("IsPedestalFinished", false);
    throw RewindDataFilesException(this);

  }
}

void EUTelPedestalNoiseProcessor::additionalMaskingLoop(LCEvent * event) {

  EUTelEventImpl * evt = static_cast<EUTelEventImpl*> (event);

  if ( evt->getEventType() == kEORE ) finalizeProcessor(true);
  if ( ( _lastEvent != -1 ) && ( _iEvt >= _lastEvent ) ) finalizeProcessor(true);
  if ( _iEvt < _firstEvent ) {
    ++_iEvt;
    throw SkipEventException(this);
  }
  
  // keep the user updated
  if ( _iEvt % 10 == 0 ) 
    streamlog_out( MESSAGE4 ) << "Processing event " 
			      << setw(6) << setiosflags(ios::right) << event->getEventNumber() << " in run "
			      << setw(6) << setiosflags(ios::right) << setfill('0')  << event->getRunNumber() << setfill(' ')
			      << " (Total = " << setw(10) << _iEvt << ")" << resetiosflags(ios::left) 
			      << " - loop " << _iLoop <<  endl;
  
  // let me get the rawDataCollection. This is should contain a TrackerRawDataObject
  // for each detector plane in the telescope.
  try {
    LCCollectionVec *collectionVec = dynamic_cast < LCCollectionVec * >(evt->getCollection (_rawDataCollectionName));
    
    for (int iDetector = 0; iDetector < _noOfDetector; iDetector++) {
      // get the TrackerRawData object from the collection for this detector
      TrackerRawData *trackerRawData = dynamic_cast < TrackerRawData * >(collectionVec->getElementAt (iDetector));
      ShortVec adcValues = trackerRawData->getADCValues ();
      for ( unsigned int iPixel = 0 ; iPixel < adcValues.size(); iPixel++ ) {
	if ( _status[iDetector][iPixel] == EUTELESCOPE::GOODPIXEL ) {
	  float correctedValue = adcValues[iPixel] - _pedestal[iDetector][iPixel];
	  float threshold      = _noise[iDetector][iPixel] *  (0.5 * _hitRejectionCut );
	  if ( correctedValue > threshold ) {
	    _hitCounter[iDetector][iPixel]++;
	  }
	}
      }
    }
    ++_iEvt;
  } catch (DataNotAvailableException& e) {
    streamlog_out ( WARNING2 ) << "No input collection " << _rawDataCollectionName << " is not available in the current event" << endl;
  }
}
