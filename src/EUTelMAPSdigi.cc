// -*- mode: c++; mode: auto-fill; mode: flyspell-prog; -*-
// Author Aleksander Zarnecki, University of Warsaw <mailto:zarnecki@fuw.edu.pl>
// Version $Id: EUTelMAPSdigi.cc,v 1.2 2008-11-11 19:12:32 zarnecki Exp $
/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

// built only if GEAR is used
#ifdef USE_GEAR

// eutelescope includes ".h"
#include "EUTelMAPSdigi.h"
#include "EUTelRunHeaderImpl.h"
#include "EUTelEventImpl.h"
#include "EUTELESCOPE.h"
#include "EUTelExceptions.h"

// Include for Track Detailed Simulation (temporarily in Eutelescope)

// #include "TDSPixelsChargeMap.h"

// marlin includes ".h"
#include "marlin/Processor.h"
#include "marlin/Global.h"

// gear includes <.h>
#include <gear/GearMgr.h>
#include <gear/SiPlanesParameters.h>

// aida includes <.h>
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
#include <marlin/AIDAProcessor.h>
#include <AIDA/IHistogramFactory.h>
#include <AIDA/IHistogram1D.h>
#include <AIDA/IHistogram2D.h>
#include <AIDA/IHistogram3D.h>
#include <AIDA/ITree.h>
#endif

// lcio includes <.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/SimTrackerHitImpl.h>
#include <UTIL/CellIDDecoder.h>

// system includes <>
#include <string>
#include <vector>
#include <algorithm>
#include <memory>
#include <iostream>
#include <iomanip>

using namespace std;
using namespace marlin;
using namespace gear;
using namespace eutelescope;

// definition of static members mainly used to name histograms
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
std::string EUTelMAPSdigi::_hitHistoLocalName          = "HitHistoLocal";
std::string EUTelMAPSdigi::_hitHistoTelescopeName      = "HitHistoTelescope";
std::string EUTelMAPSdigi::_pixelHistoName             = "PixelHisto";
#endif

EUTelMAPSdigi::EUTelMAPSdigi () : Processor("EUTelMAPSdigi") {

  // modify processor description
  _description =
    "EUTelMAPSdigi is used to calculate expected charge sharing between pixels"
    "for simulated Mokka hits";

  registerInputCollection(LCIO::SIMTRACKERHIT,"SimHitCollectionName",
                          "Simulated (Mokka) hit collection name",
                          _simhitCollectionName, string ( "telEUTelescopeCollection" ));
 
  registerOutputCollection(LCIO::TRACKERHIT,"PixelCollectionName",
                           "Collection name for simulated raw data",
                           _pixelCollectionName, string ("simEUTelescopeData"));

  // Track Detailes Simulation parameters

  registerProcessorParameter ("MaxStepInCharge",
                              "Maximum charge per simulation step [e]",
                              _maxStepInCharge,  static_cast < double > (300.));

  registerProcessorParameter ("MaxStepInLength",
                              "Maximum step length for single hit",
                              _maxStepInLength,  static_cast < double > (0.001));
  
  registerProcessorParameter ("IntegMaxNumberPixelsAlongL",
                              "Maximum diffusion range in pixels along sensor length ",
                              _integMaxNumberPixelsAlongL,  static_cast < int > (9));
  
  registerProcessorParameter ("IntegMaxNumberPixelsAlongW",
                              "Maximum diffusion range in pixels along sensor width ",
                              _integMaxNumberPixelsAlongW,  static_cast < int > (9));
  
  registerProcessorParameter ("ChargeAttenuationLength",
                              "Charge attenuation length in diffusion",
                              _chargeAttenuationLength,  static_cast < double > (0.038));
  

  registerProcessorParameter ("ChargeReflectedContribution",
                              "Charge reflection coefficient",
                              _chargeReflectedContribution,  static_cast < double > (1.0));
  
  registerProcessorParameter ("GSL function calls",
                              "Number of function calls in one GSL integration",
                              _gslFunctionCalls,  static_cast < int > (5000));
  
  registerProcessorParameter ("IntegPixelSegmentsAlongL",
                              "Number of bins along sensor lenght for storing integration results",
                              _integPixelSegmentsAlongL,  static_cast < int > (16));
  

  registerProcessorParameter ("IntegPixelSegmentsAlongW",
                              "Number of bins along sensor width for storing integration results",
                              _integPixelSegmentsAlongW,  static_cast < int > (16));
  

  registerProcessorParameter ("IntegPixelSegmentsAlongH",
                              "Number of bins along sensor depth for storing integration results",
                              _integPixelSegmentsAlongH,  static_cast < int > (16));
  

  registerProcessorParameter ("Use common integration storage",
         "Common integration storage can be used for all sensors if they have same geometry",
                             _useCommonIntegrationStorage,  static_cast < bool > (true));


  registerProcessorParameter ("IonizationEnergy",
                              "Ionization energy in silicon [eV]",
                              _ionizationEnergy,  static_cast < double > (3.6));
  


  registerProcessorParameter ("DepositedChargeScaling",
                              "Scaling of the deposited charge",
                             _depositedChargeScaling ,  static_cast < double > (1.0));
  


  registerProcessorParameter ("ApplyPoissonSmearing",
                              "Flag for applying Poisson smearing of the collected charge",
                              _applyPoissonSmearing,  static_cast < bool > (true));


  registerProcessorParameter ("AdcGain",
                              "ADC gain  in ADC counts per unit charge",
                              _adcGain,  static_cast < double > (0.01));
  


  registerProcessorParameter ("AdcGainVariation",
                              "ADC gain variation ",
                             _adcGainVariation ,  static_cast < double > (0.));
  


  registerProcessorParameter ("AdcNoise",
                              "ADC noise in ADC counts",
                             _adcNoise ,  static_cast < double > (1.0));
  


  registerProcessorParameter ("AdcOffset",
                              "Constant pedestal value in ADC couts",
                             _adcOffset ,  static_cast < double > (0.));
  

  registerProcessorParameter ("ZeroSuppressionThreshold",
                              "Threshold (in ADC counts) for removing empty pixels",
                             _zeroSuppressionThreshold ,  static_cast < double > (3.));
  


// Other processor parameters


  registerProcessorParameter ("SingleHitOutput",
                              "Flag controling pixel list output after each Mokka hit",
                              _singleHitOutput,  static_cast < bool > (true));


  registerProcessorParameter ("DebugEventCount",
                              "Print out every DebugEnevtCount event",
                              _debugCount,  static_cast < int > (10));

}


void EUTelMAPSdigi::init() {
  // this method is called only once even when the rewind is active
  // usually a good idea to
  printParameters ();

  // set to zero the run and event counters
  _iRun = 0;
  _iEvt = 0;


  // check if Marlin was built with GEAR support or not
#ifndef USE_GEAR

  streamlog_out ( ERROR4 ) <<  "Marlin was not built with GEAR support." << endl
                           <<  "You need to install GEAR and recompile Marlin with -DUSE_GEAR before continue." << endl;

  // I'm thinking if this is the case of throwing an exception or
  // not. This is a really error and not something that can
  // exceptionally happens. Still not sure what to do
  exit(-1);

#else

  // check if the GEAR manager pointer is not null!
  if ( Global::GEAR == 0x0 ) {
    streamlog_out ( ERROR4 ) <<  "The GearMgr is not available, for an unknown reason." << endl;
    exit(-1);
  }

  _siPlanesParameters  = const_cast<SiPlanesParameters* > (&(Global::GEAR->getSiPlanesParameters()));
  _siPlanesLayerLayout = const_cast<SiPlanesLayerLayout*> ( &(_siPlanesParameters->getSiPlanesLayerLayout() ));

  _histogramSwitch = true;

#endif

  // Book common integration storage for all sensors, if this was
  // requested 

  if (_useCommonIntegrationStorage)
     {
     _integrationStorage = new TDSIntegrationStorage(_integPixelSegmentsAlongL,
                          _integPixelSegmentsAlongW, _integPixelSegmentsAlongH);
      streamlog_out( MESSAGE4 )  << " Common TDS integration storage initialized  " << endl;
     }

}

void EUTelMAPSdigi::processRunHeader (LCRunHeader * rdr) {

// This processor should only be used for simulated data
// Check if the input was generated with Mokka


streamlog_out( MESSAGE4 )  << "  Run : " << rdr->getRunNumber() 
                           << " - "      << rdr->getDetectorName() 
                           << ":  "      << rdr->getDescription()  << endl ;

string simulator = rdr->parameters().getStringVal("SimulatorName");

if(simulator != "")
  {
   streamlog_out( MESSAGE4 )  << simulator << " input file recognized, version  " 
                              << rdr->parameters().getStringVal("SimulatorVersion")  << endl;
  }
else
   {
   streamlog_out ( ERROR4 ) << "Error during consistency check: " << endl
                            << "EUTelMAPSdigi processor can only run on simulated data"  << endl
                            << "but SimulatorName not set !?" << endl;
   exit(-1);
   }


  auto_ptr<EUTelRunHeaderImpl> header ( new EUTelRunHeaderImpl (rdr) );
  header->addProcessor( type() );

  // the run header contains the number of detectors. This number
  // should be in principle be the same as the number of layers in the
  // geometry description. But it is not always set for simulation output
  if ( header->getNoOfDetector() != _siPlanesParameters->getSiPlanesNumber() ) 
    streamlog_out ( WARNING0 ) << "Warning during the geometry consistency check: " << endl
          << "The run header says there are " << header->getNoOfDetector() << " silicon detectors " << endl
          << "The GEAR description says     " << _siPlanesParameters->getSiPlanesNumber() << " silicon planes" << endl;
 

  // this is the right place also to check the geometry ID. This is a
  // unique number identifying each different geometry used at the
  // beam test. The same number should be saved in the run header and
  // in the xml file.  But it is not always set for simulation output


  if ( header->getGeoID() != _siPlanesParameters->getSiPlanesID() ) 
    streamlog_out ( WARNING0 ) <<  "Warning during the geometry consistency check: " << endl
                             << "The run header says the GeoID is " << header->getGeoID() << endl
                             << "The GEAR description says is     " << _siPlanesParameters->getSiPlanesID() << endl;

  

  // now book histograms plz...
  if ( isFirstEvent() )  bookHistos();

  // increment the run counter
  ++_iRun;
}


void EUTelMAPSdigi::processEvent (LCEvent * event) {

  bool debug = ( _debugCount>0 && _iEvt%_debugCount == 0);

  if (_iEvt % 10 == 0  || debug)
    streamlog_out( MESSAGE4 ) << "Processing event "
                              << setw(6) << setiosflags(ios::right) << event->getEventNumber() << " in run "
                              << setw(6) << setiosflags(ios::right) << setfill('0')  << event->getRunNumber() << setfill(' ')
                              << " (Total = " << setw(10) << _iEvt << ")" << resetiosflags(ios::left) << endl;

   ++_iEvt;

  EUTelEventImpl * evt = static_cast<EUTelEventImpl*> (event) ;

// Event type not set in the Mokka output
  if ( evt->getEventType() == kEORE ) {
    streamlog_out ( DEBUG4 ) << "EORE found: nothing else to do." << endl;
    return;
  } 


  try {

    LCCollectionVec * simhitCollection   = static_cast<LCCollectionVec*> (event->getCollection( _simhitCollectionName ));

    int detectorID    = -99; // it's a non sense
    int oldDetectorID = -100;
    int mapDetectorID = -100;

    int    layerIndex = -99;
    double xZero = 0., yZero = 0., zZero = 0. ;
    double xSize = 0., ySize = 0.;
    double zThickness = 0.;
    double xPitch = 0., yPitch = 0.;
    double xPointing[2] = { 1., 0. }, yPointing[2] = { 1., 0. };

    for ( int iHit = 0; iHit < simhitCollection->getNumberOfElements(); iHit++ ) {

      SimTrackerHitImpl * simhit   = static_cast<SimTrackerHitImpl*> ( simhitCollection->getElementAt(iHit) );

      // there could be several hits belonging to the same
      // detector. So update the geometry information only if this new
      // cluster belongs to a different detector.
      detectorID = simhit->getCellID();

      if ( detectorID != oldDetectorID ) {

	oldDetectorID = detectorID;

        if ( _conversionIdMap.size() != (unsigned) _siPlanesParameters->getSiPlanesNumber() ) {
          // first of all try to see if this detectorID already belong to
          if ( _conversionIdMap.find( detectorID ) == _conversionIdMap.end() ) {
            // this means that this detector ID was not already inserted,
            // so this is the right place to do that
            for ( int iLayer = 0; iLayer < _siPlanesLayerLayout->getNLayers(); iLayer++ ) {
              if ( _siPlanesLayerLayout->getSensitiveID(iLayer) == detectorID ) {
                _conversionIdMap.insert( make_pair( detectorID, iLayer ) );
                streamlog_out ( DEBUG4 ) << "Sensitive layer ID = " << detectorID << " found in GEAR" << endl;
               break;
              }
            }
          }
        }

        // perfect! The full geometry description is now coming from the
        // GEAR interface. Let's keep the finger xed!
        layerIndex   = _conversionIdMap[detectorID];

        xZero        = _siPlanesLayerLayout->getSensitivePositionX(layerIndex); // mm
        yZero        = _siPlanesLayerLayout->getSensitivePositionY(layerIndex); // mm
        zZero        = _siPlanesLayerLayout->getSensitivePositionZ(layerIndex); // mm
        zThickness   = _siPlanesLayerLayout->getSensitiveThickness(layerIndex); // mm
        xPitch       = _siPlanesLayerLayout->getSensitivePitchX(layerIndex);    // mm
        yPitch       = _siPlanesLayerLayout->getSensitivePitchY(layerIndex);    // mm
        xSize        = _siPlanesLayerLayout->getSensitiveSizeX(layerIndex);     // mm
        ySize        = _siPlanesLayerLayout->getSensitiveSizeY(layerIndex);     // mm
        xPointing[0] = _siPlanesLayerLayout->getSensitiveRotation1(layerIndex); // was -1 ;
        xPointing[1] = _siPlanesLayerLayout->getSensitiveRotation2(layerIndex); // was  0 ;
        yPointing[0] = _siPlanesLayerLayout->getSensitiveRotation3(layerIndex); // was  0 ;
        yPointing[1] = _siPlanesLayerLayout->getSensitiveRotation4(layerIndex); // was -1 ;

        if (  ( xPointing[0] == xPointing[1] ) && ( xPointing[0] == 0 ) ) {
          streamlog_out ( ERROR4 ) << "Detector " << detectorID << " has a singular rotation matrix. Sorry for quitting" << endl;
        }

        if (  ( yPointing[0] == yPointing[1] ) && ( yPointing[0] == 0 ) ) {
          streamlog_out ( ERROR4 ) << "Detector " << detectorID << " has a singular rotation matrix. Sorry for quitting" << endl;

        }

     }  // End of if ( detectorID != oldDetectorID ) 


      // get the position and momentum of particle generating energy deposit

     const double* hitpos = simhit->getPosition();
     const float*  hitmom = simhit->getMomentum();

  if (debug)
    streamlog_out( DEBUG4 ) << "SimHit in global frame  at X = " <<  hitpos[0] 
                                                     <<  " Y = " <<  hitpos[1]
                                                     <<  " Z = " <<  hitpos[2] << endl;

     // 2D histograms of simulated hits in global (telescope) reference frame

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
      if ( _histogramSwitch ) {
        string tempHistoName;
        {
          stringstream ss;
          ss << _hitHistoTelescopeName << "-" << layerIndex ;
          tempHistoName = ss.str();
        }
        AIDA::IHistogram2D * histo2D = dynamic_cast<AIDA::IHistogram2D*> (_aidaHistoMap[ tempHistoName ] );
        if ( histo2D ) histo2D->fill( hitpos[0], hitpos[1] );
        else {
          streamlog_out ( ERROR1 )  << "Not able to retrieve histogram pointer for " << tempHistoName
                                    << ".\nDisabling histogramming from now on " << endl;
          _histogramSwitch = false;
        }

      }
#endif

     // transform position and momentum to local (sensor) reference frame
     // sensor rotation is described by two vectors (from implementation in HitMaker):
     // (xPointing[0],yPointing[0],0) gives direction of sensor X axis in telescope reference frame
     // (xPointing[1],yPointing[1],0) gives direction of sensor Y axis in telescope reference frame
     // direction of Z axis is unchanged, only pointing can flip


      // first the translation to the frame with origin at the sensor center

      double senspos[3];

      senspos[0] = hitpos[0] -  xZero;
      senspos[1] = hitpos[1] -  yZero;
      senspos[2] = hitpos[2] -  zZero;

     // now move to the sensor corner. Sign determines which way to move
     
     // X axis
      double sign = 1;

      if      ( xPointing[0] < 0 )       sign = -1 ;
      else if ( xPointing[0] > 0 )       sign =  1 ;
      else {
        if       ( xPointing[1] < 0 )    sign = -1 ;
        else if  ( xPointing[1] > 0 )    sign =  1 ;
      }

      senspos[0] +=  sign * xSize/2;

      // Y axis

      if      ( yPointing[0] < 0 )       sign = -1 ;
      else if ( yPointing[0] > 0 )       sign =  1 ;
      else {
        if       ( yPointing[1] < 0 )    sign = -1 ;
        else if  ( yPointing[1] > 0 )    sign =  1 ;
      }

      senspos[1] +=  sign * ySize/2;

      // Z axis direction from determinant of the XY rotation matrix

      sign = xPointing[0] * yPointing[1] - xPointing[1] * yPointing[0] ;

      senspos[2] += sign * zThickness/2;

     // Last step: rotation taking into account sensor orientation.
     // Reverse rotation has to be aplied!
     // "sign" is the determinant of the rotation matrix


     _localPosition[0] =  yPointing[1]/sign * senspos[0] - xPointing[1]/sign * senspos[1] ;
     _localPosition[1] = -yPointing[0]/sign * senspos[0] + xPointing[0]/sign * senspos[1] ;
     _localPosition[2] = sign * senspos[2] ;

  if (debug)
    streamlog_out( DEBUG4 ) << "Position in local frame at X = " << _localPosition[0] 
                                                     <<  " Y = " << _localPosition[1]
                                                     <<  " Z = " << _localPosition[2] << endl;

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
      string tempHistoName;
      if ( _histogramSwitch ) {
        {
          stringstream ss;
          ss << _hitHistoLocalName << "-" << layerIndex ;
          tempHistoName = ss.str();
        }
        if ( AIDA::IHistogram2D* histo = dynamic_cast<AIDA::IHistogram2D*>(_aidaHistoMap[ tempHistoName ]) )
          histo->fill(_localPosition[0], _localPosition[1]);
        else {
          streamlog_out ( ERROR1 )  << "Not able to retrieve histogram pointer for " << tempHistoName
                                    << ".\nDisabling histogramming from now on " << endl;
          _histogramSwitch = false;
        }

      }
#endif

      // momentum transformation (only rotation)

 
     _localMomentum[0] =  yPointing[1]/sign * hitmom[0] - xPointing[1]/sign * hitmom[1] ;
     _localMomentum[1] = -yPointing[0]/sign * hitmom[0] + xPointing[0]/sign * hitmom[1] ;
     _localMomentum[2] = sign * hitmom[2] ;

    // Normalized direction vector 

    double ptot;
     ptot=sqrt(_localMomentum[0]*_localMomentum[0]+_localMomentum[1]*_localMomentum[1]+_localMomentum[2]*_localMomentum[2]);


    for(int idir=0; idir<3; idir++)
              _localDirection[idir]=_localMomentum[idir]/ptot;


     // Sensor size in the local (sensor) frame (X can swap with Y due to rotation!)

 
     _localSize[0] =  yPointing[1]/sign * xSize - xPointing[1]/sign * ySize ;
     if(_localSize[0] < 0)_localSize[0] = -_localSize[0];

     _localSize[1] = -yPointing[0]/sign * xSize + xPointing[0]/sign * ySize ;
     if(_localSize[1] < 0)_localSize[1] = -_localSize[1];

     

     _localSize[2] = zThickness;

    // From EUTelHitMaker code it looks like sensor pitches given in
    // GEAR are already in the local frame of reference:

     _localPitch[0] = xPitch;
     _localPitch[1] = yPitch;


    // Remaining information on the simulated hit:

     _mokkaDeposit =simhit->getdEdx();
     _mokkaPath=simhit->getPathLength();

    ///////////////////////////////////////////////////////////////////////////
    //
    // Now we are in the local coordinate frame, and we should have 
    // all information needed for digitization 
    //
    //  Hit position in the local frame of reference:  _localPosition[3]
    //  Particle momentum [MeV] in local frame:  _localMomentum[3]
    //  Particle direction in local frame:  _localDirection[3]
    //  Total energy deposit:  _mokkaDeposit
    //  Total particle path lenght:  _mokkaPath
    //
    //  Sensor boundaries in the local frame are:
    //    (0,0,0) to _localSize[3]
    //
    //  Sensor pitch is given by _localPitch[2] (X and Y only)
    //
    //////////////////////////////////////////////////////////////////////////



    // We have to initialize TDS pixel charge map if not done yet

      if ( detectorID != mapDetectorID ) {
          mapDetectorID = detectorID;

        if ( _pixelChargeMapCollection.size() != (unsigned) _siPlanesParameters->getSiPlanesNumber() ) {
          // first of all try to see if this detectorID already belong to
          if ( _pixelChargeMapCollection.find( detectorID ) == _pixelChargeMapCollection.end() ) {
            // this means that the map for was not defined yet,
            // so this is the right place to do that

	    // current TDS requires that negative sensor thickness is given!
          _pixelChargeMap = new TDSPixelsChargeMap(_localSize[0],_localSize[1],-_localSize[2]);

          _pixelChargeMapCollection.insert( make_pair( detectorID,_pixelChargeMap) );

	  // Set additional geometry parameters
          _pixelChargeMap->setPixelLength(_localPitch[0]);
          _pixelChargeMap->setPixelWidth(_localPitch[1]);

          // Set parameters of the charge diffusion model (input parameters)

          _pixelChargeMap->setLambda(_chargeAttenuationLength);
          _pixelChargeMap->setReflectedContribution(_chargeReflectedContribution);

          // Initialize integration parameters (all taken from input
          // processor parameters)

          _pixelChargeMap->initializeIntegration(_maxStepInLength,_maxStepInCharge,
                _integMaxNumberPixelsAlongL, _integMaxNumberPixelsAlongW, _gslFunctionCalls);


          streamlog_out( MESSAGE4 )  << " Pixel map initialized  for detectorID " 
                                       << detectorID << endl;


          // Initialize integration storage

          if (!_useCommonIntegrationStorage)
	    {
	      _integrationStorage = new TDSIntegrationStorage(_integPixelSegmentsAlongL,
                                   _integPixelSegmentsAlongW, _integPixelSegmentsAlongH);

            streamlog_out( MESSAGE4 )  << " TDS integration storage initialized  for detectorID " 
                                       << detectorID << endl;

            }

	  _pixelChargeMap->setPointerToIntegrationStorage(_integrationStorage);
	  }
	}

       // Get pointer to the TDS pixel charge map
        _pixelChargeMap = _pixelChargeMapCollection[detectorID];

      }


    // First step is to calculate track segment starting point and end point
    // One has to check if not going outside the sensor when moving by +/- path/2

     double pathFrac = CheckPathLimits();

    if (debug && pathFrac!=1.)
       streamlog_out( DEBUG4 ) << "Track path inside sensor reduced by factor " << pathFrac << endl;

    //
    //  Here comes the core of the algorithm
    //

    // Convert Mokka deposit [GeV] to charge in units of elementary
    // charge 

    _mokkaDeposit *= 1000000000./_ionizationEnergy;

    // Set the Geant/Mokka step to be considered
    // TDS requires position in depth to be negative ! Shift it by
    // thickness to keep direction unchanged

    TDSStep step(_localPosition[0], _localPosition[1], _localPosition[2]-_localSize[2],
		 _localDirection[0],_localDirection[1],_localDirection[2],
                 _mokkaPath, _mokkaDeposit);

    // distribute charge among pixels (here all the work is done!)
 
    _pixelChargeMap->update(step);

     if (debug)
	streamlog_out( DEBUG4 ) << "Total charge deposited: " << _pixelChargeMap->getTotalCharge()
                                << " (Mokka deposit: " << _mokkaDeposit << " )" << endl;

     // Charge processing

     // Scaling of charge deposited at each pixel

     if(_depositedChargeScaling!=1.)
       _pixelChargeMap->scaleCharge(_depositedChargeScaling);

 
     // Poisson smearing (if requested)

     if(_applyPoissonSmearing)
       {
        _pixelChargeMap->applyPoissonFluctuations(false);

       if (debug)
	 streamlog_out( DEBUG4 ) << "Charge after Poisson fluctuations: " 
                                 << _pixelChargeMap->getTotalCharge() << endl;

       }

     // ADC gain, noise, pedestal

     _pixelChargeMap->applyGain(_adcGain, _adcGainVariation, _adcNoise, _adcOffset);

     // TDS allow for negative charges and negative thresholds.
     // However we assume here that threshold has to be
     // positive. If zero or negative threshold value is set, not
     // threshold correction is applied. 

     if(_zeroSuppressionThreshold>0)
        _pixelChargeMap->applyThresholdCut(_zeroSuppressionThreshold);

     if (debug)
	 streamlog_out( DEBUG4 ) << "Total signal after gain and ZS: " 
                                 << _pixelChargeMap->getTotalCharge() << endl;

    // The last task is to put pixels with collected charge into
    // output data collection. This is done either for each Mokka hit (here)
    // or after the whole event is completed.

   if( _singleHitOutput)
     {
     // Get vector of fired pixels
     
     _vectorOfPixels = _pixelChargeMap->getVectorOfPixels();

     // Output to raw data structure
     //      To be added
     //      ===========
     // List all fired pixels
 
     if (debug)
        {
        int nPixel = _vectorOfPixels.size() ;

        streamlog_out( MESSAGE4 ) <<  "Detector ID = " << detectorID 
                                  << ", " << nPixel << " pixels fired " << endl;

       for(_pixelIterator = _vectorOfPixels.begin(); _pixelIterator != _vectorOfPixels.end(); _pixelIterator++)
	   streamlog_out( DEBUG4 ) <<  " Pixel at  (" << _pixelIterator->getIndexAlongL() 
                                   <<  "," <<  _pixelIterator->getIndexAlongW() 
				   <<  ") with Q = " << _pixelIterator->getCharge()  << endl;

// Store pixel charges in the histogram
// Only for events with debug output!

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
      string tempHistoName;
      if ( _histogramSwitch ) {

        {
          stringstream ss;
          ss << _pixelHistoName << "-" << layerIndex ;
          tempHistoName = ss.str();
        }
        if ( AIDA::IHistogram2D* histo = dynamic_cast<AIDA::IHistogram2D*>(_aidaHistoMap[ tempHistoName ]) )
           {
           for(_pixelIterator = _vectorOfPixels.begin(); _pixelIterator != _vectorOfPixels.end(); _pixelIterator++)
              {
              double xpixel = _pixelIterator->getIndexAlongL() + 0.5;
              double ypixel = _pixelIterator->getIndexAlongW() + 0.5;
              histo->fill(xpixel,ypixel,_pixelIterator->getCharge());
              }
           }
        else {
          streamlog_out ( ERROR1 )  << "Not able to retrieve histogram pointer for " << tempHistoName
                                    << ".\nDisabling histogramming from now on " << endl;
          _histogramSwitch = false;
        }

      }
#endif



	}  // end of if (debug)

       // clear pixel map after output 
       _pixelChargeMap->clear();
       _vectorOfPixels.clear();
      }


    }
       // end of loop over SimTrackerHit collection

// Output all collected pixels 
// If single hit flag is not set, then hits are written out after the whole event

  if( ! _singleHitOutput)
     {

     // Loop over defined detectors
 
     std::map< int,  TDSPixelsChargeMap *>::iterator mapIterator;

      for(mapIterator = _pixelChargeMapCollection.begin(); 
          mapIterator != _pixelChargeMapCollection.end(); mapIterator++)
        {
        detectorID = mapIterator->first;
        layerIndex   = _conversionIdMap[detectorID];

        _pixelChargeMap = mapIterator->second;

        _vectorOfPixels = _pixelChargeMap->getVectorOfPixels();

        // Output to raw data structure
        //      To be added
        //      ===========


     // List all fired pixels
 
     if (debug)
        {
        int nPixel = _vectorOfPixels.size() ;

        streamlog_out( MESSAGE4 ) <<  "Detector ID = " << detectorID 
                                  << ", " << nPixel << " pixels fired, " 
				  << "total charge deposited: " << _pixelChargeMap->getTotalCharge() << endl;

       for(_pixelIterator = _vectorOfPixels.begin(); _pixelIterator != _vectorOfPixels.end(); _pixelIterator++)
	   streamlog_out( DEBUG4 ) <<  " Pixel at  (" << _pixelIterator->getIndexAlongL() 
                                   <<  "," <<  _pixelIterator->getIndexAlongW() 
				   <<  ") with Q = " << _pixelIterator->getCharge()  << endl;

// Store pixel charges in the histogram
// Only for events with debug output!

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
      string tempHistoName;
      if ( _histogramSwitch ) {

        {
          stringstream ss;
          ss << _pixelHistoName << "-" << layerIndex ;
          tempHistoName = ss.str();
        }
        if ( AIDA::IHistogram2D* histo = dynamic_cast<AIDA::IHistogram2D*>(_aidaHistoMap[ tempHistoName ]) )
           {
           for(_pixelIterator = _vectorOfPixels.begin(); _pixelIterator != _vectorOfPixels.end(); _pixelIterator++)
              {
              double xpixel = _pixelIterator->getIndexAlongL() + 0.5;
              double ypixel = _pixelIterator->getIndexAlongW() + 0.5;
              histo->fill(xpixel,ypixel,_pixelIterator->getCharge());
              }
           }
        else {
          streamlog_out ( ERROR1 )  << "Not able to retrieve histogram pointer for " << tempHistoName
                                    << ".\nDisabling histogramming from now on " << endl;
          _histogramSwitch = false;
        }

      }
#endif


       }// end of if (debug)

      // clear stored pixel map
       _pixelChargeMap->clear();
       _vectorOfPixels.clear();

      } // end of loop over defined pixel maps


     } // end of  if( ! _singleHitOutput)



    if ( isFirstEvent() ) _isFirstEvent = false;
  } catch (DataNotAvailableException& e  ) {
    streamlog_out  ( WARNING2 ) <<  "No input collection found on event " << event->getEventNumber()
                                << " in run " << event->getRunNumber() << endl;
  }

}

void EUTelMAPSdigi::end() {

  streamlog_out ( MESSAGE4 )  << "Successfully finished" << endl;
}

void EUTelMAPSdigi::bookHistos() {

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)

  try {
    streamlog_out ( MESSAGE4 ) <<  "Booking histograms" << endl;

    string tempHistoName;

    // histograms are grouped into folders named after the
    // detector. This requires to loop on detector now.
    for (int iDet = 0 ; iDet < _siPlanesParameters->getSiPlanesNumber(); iDet++) {

      string basePath;
      {
        stringstream ss ;
        ss << "plane-" << iDet;
        basePath = ss.str();
      }
      AIDAProcessor::tree(this)->mkdir(basePath.c_str());
      basePath = basePath + "/";

      {
        stringstream ss ;
        ss <<  _hitHistoLocalName << "-" << iDet ;
        tempHistoName = ss.str();
      }



      double xMin =  0;
      double xMax =  _siPlanesLayerLayout->getSensitiveSizeX ( iDet );

      double yMin =  0;
      double yMax =  _siPlanesLayerLayout->getSensitiveSizeY ( iDet );

      int xNBin =  _siPlanesLayerLayout->getSensitiveNpixelX( iDet );
      int yNBin =  _siPlanesLayerLayout->getSensitiveNpixelY( iDet );


      AIDA::IHistogram2D * hitHistoLocal = AIDAProcessor::histogramFactory(this)->createHistogram2D( (basePath + tempHistoName).c_str(),
                                                                                                     xNBin, xMin, xMax, yNBin, yMin, yMax );
      if ( hitHistoLocal ) {
        hitHistoLocal->setTitle("Hit map in the detector local frame of reference");
        _aidaHistoMap.insert( make_pair( tempHistoName, hitHistoLocal ) );
      } else {
        streamlog_out ( ERROR1 )  << "Problem booking the " << (basePath + tempHistoName) << ".\n"
                                  << "Very likely a problem with path name. Switching off histogramming and continue w/o" << endl;
        _histogramSwitch = false;
      }

      // 2 should be enough because it
      // means that the sensor is wrong
      // by all its size.
      double safetyFactor = 2.0;

      xMin = safetyFactor * ( _siPlanesLayerLayout->getSensitivePositionX( iDet ) -
                                     ( 0.5 * _siPlanesLayerLayout->getSensitiveSizeX ( iDet ) ));
      xMax = safetyFactor * ( _siPlanesLayerLayout->getSensitivePositionX( iDet ) +
                                     ( 0.5 * _siPlanesLayerLayout->getSensitiveSizeX ( iDet )));

      yMin = safetyFactor * ( _siPlanesLayerLayout->getSensitivePositionY( iDet ) -
                                     ( 0.5 * _siPlanesLayerLayout->getSensitiveSizeY ( iDet )));
      yMax = safetyFactor * ( _siPlanesLayerLayout->getSensitivePositionY( iDet ) +
                                     ( 0.5 * _siPlanesLayerLayout->getSensitiveSizeY ( iDet )) );

      xNBin = static_cast< int > ( safetyFactor ) * _siPlanesLayerLayout->getSensitiveNpixelX( iDet );
      yNBin = static_cast< int > ( safetyFactor ) * _siPlanesLayerLayout->getSensitiveNpixelY( iDet );

      {
        stringstream ss ;
        ss <<  _hitHistoTelescopeName << "-" << iDet ;
        tempHistoName = ss.str();
      }
      AIDA::IHistogram2D * hitHistoTelescope =
        AIDAProcessor::histogramFactory(this)->createHistogram2D( ( basePath + tempHistoName ).c_str(),
                                                                  xNBin, xMin, xMax, yNBin, yMin, yMax );

      if ( hitHistoTelescope ) {
        hitHistoTelescope->setTitle("Hit map in the telescope frame of reference");
        _aidaHistoMap.insert( make_pair ( tempHistoName, hitHistoTelescope ) );
      } else {
        streamlog_out ( ERROR1 )  << "Problem booking the " << (basePath + tempHistoName) << ".\n"
                                  << "Very likely a problem with path name. Switching off histogramming and continue w/o" << endl;
        _histogramSwitch = false;
      }
 
     // Pixel map to show generated cluster shapes

      {
        stringstream ss ;
        ss <<  _pixelHistoName << "-" << iDet ;
        tempHistoName = ss.str();
      }

      xNBin =  _siPlanesLayerLayout->getSensitiveNpixelX( iDet );
      yNBin =  _siPlanesLayerLayout->getSensitiveNpixelY( iDet );

      xMin =  0;
      xMax =  xNBin;

      yMin =  0;
      yMax =  yNBin;

      AIDA::IHistogram2D * pixelHisto = AIDAProcessor::histogramFactory(this)->createHistogram2D( (basePath + tempHistoName).c_str(),
                                                                                                     xNBin, xMin, xMax, yNBin, yMin, yMax );
      if ( pixelHisto ) {
        pixelHisto->setTitle("Pixel map");
        _aidaHistoMap.insert( make_pair( tempHistoName, pixelHisto ) );
      } else {
        streamlog_out ( ERROR1 )  << "Problem booking the " << (basePath + tempHistoName) << ".\n"
                                  << "Very likely a problem with path name. Switching off histogramming and continue w/o" << endl;
        _histogramSwitch = false;
      }

   }

  } catch (lcio::Exception& e ) {

    streamlog_out ( ERROR1 ) << "No AIDAProcessor initialized. Type q to exit or c to continue without histogramming" << endl;
    string answer;
    while ( true ) {
      streamlog_out ( ERROR1 ) <<  "[q]/[c]" << endl;
      cin >> answer;
      transform( answer.begin(), answer.end(), answer.begin(), ::tolower );
      if ( answer == "q" ) {
        exit(-1);
      } else if ( answer == "c" )
        _histogramSwitch = false;
      break;
    }
  }
#endif
}
//
// ===============================================================================
//
//  Private function members
//


double EUTelMAPSdigi::CheckPathLimits()
{
  // For most tracks the whole path should fit inside a sensor

  double pathStart=0.;
  double pathEnd=1.;

  // Check all sensor dimensions ! 
  // (Probably checking Z would be enough :-)

  for(int idim=0; idim<3 ; idim++)
    {
    if(_localPosition[idim]+(pathStart-0.5)*_mokkaPath*_localDirection[idim] < 0)
                 pathStart= 0.5 - _localPosition[idim]/(_mokkaPath*_localDirection[idim]);

    if(_localPosition[idim]+(pathStart-0.5)*_mokkaPath*_localDirection[idim] > _localSize[idim])
                 pathStart= 0.5 + (_localSize[idim]- _localPosition[idim])/(_mokkaPath*_localDirection[idim]);

    if(_localPosition[idim]+(pathEnd-0.5)*_mokkaPath*_localDirection[idim] > _localSize[idim])
                 pathEnd= 0.5 + (_localSize[idim]- _localPosition[idim])/(_mokkaPath*_localDirection[idim]);

    if(_localPosition[idim]+(pathEnd-0.5)*_mokkaPath*_localDirection[idim] < 0)
                 pathEnd= 0.5 - _localPosition[idim]/(_mokkaPath*_localDirection[idim]);

    }

  // 

  if(  pathStart < 0. || pathStart > 1. || pathEnd < 0. || pathEnd > 1.)
       streamlog_out ( WARNING0 ) << "Warning in checking path limits: out of range " << endl
                    << "Start point = " << pathStart << "  End point = " << pathEnd << endl;

  if(  pathStart < 0. ) pathStart = 0.;
  if(  pathStart > 1. ) pathStart = 1.;
  if(  pathEnd < 0. ) pathEnd = 0.;
  if(  pathEnd > 1. ) pathEnd = 1.;

  // Apply corrections to track length and track position

  double midShift=(pathEnd+pathStart-1.)/2.;

  for(int idim=0; idim<3 ; idim++)
     _localPosition[idim]+=_mokkaPath*_localDirection[idim]*midShift;

  _mokkaPath*=(pathEnd - pathStart);

  return pathEnd - pathStart;

}


#endif
