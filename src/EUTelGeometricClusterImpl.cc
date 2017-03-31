#include "EUTelGeometricClusterImpl.h"
#include "EUTelGeometricPixel.h"

using namespace eutelescope;

EUTelGeometricClusterImpl::EUTelGeometricClusterImpl(IMPL::TrackerDataImpl* data): EUTelGenericSparseClusterImpl<EUTelGeometricPixel>(data)
{
  /*nothing else to do*/
} 

EUTelGeometricClusterImpl::~EUTelGeometricClusterImpl()
{
  
}

void EUTelGeometricClusterImpl::getClusterGeomInfo(float& xPos, float& yPos, float& xSize, float& ySize) const 
{
  float xMin = std::numeric_limits<float>::max(); 	//stores the largest possible value every pixel will be lower, 
  float yMin = xMin;					//so its OK for max 
  float xMax = -xMin;								
  float yMax = xMax;
  float xMinBoundary = 0;
  float xMaxBoundary = 0;   //boundary is +- 0.018
  float yMinBoundary = 0;
  float yMaxBoundary = 0;
  float old_xSize=0;
  float old_ySize=0;
  
  int xmin_int=0;
  int ymin_int=0;
  int xmax_int=0;
  int ymax_int=0;
  int x_tot=0;
  int y_tot=0;
  static int call_count=0;
  //int kenx=0;
  //int keny=0;
  
  EUTelGeometricPixel* pixel = new EUTelGeometricPixel;
  //Loop over all the pixels in the cluster
  for( unsigned int index = 0; index < size() ; index++ ){ 
    
    //Get the pixel
    getSparsePixelAt( index , pixel);
    
    //And its position
    float xCur = pixel->getPosX();
    float yCur = pixel->getPosY();

    //if(xCur<xMin||xCur>xMax){ kenx++; }
    //if(yCur<yMin||yCur>yMax){ keny++; }
    
    if( xCur < xMin )
      {
	xMin = xCur;
	xMinBoundary = pixel->getBoundaryX();
	xmin_int++;
      }
    if ( xCur > xMax ) 
      {
	xMax = xCur;
	xMaxBoundary = pixel->getBoundaryX();
	xmax_int++;
      }
    if ( yCur < yMin )
      {
	yMin = yCur;
	yMinBoundary = pixel->getBoundaryY();
	ymin_int++;
      }
    if ( yCur > yMax ) 
      {
	yMax = yCur;
	yMaxBoundary = pixel->getBoundaryY();
	ymax_int++;
      }
  }
  x_tot = xmin_int + xmax_int;
  y_tot = ymin_int + ymax_int;
  
  xSize = x_tot-1;
  ySize = y_tot-1;
  old_xSize = xMax + xMaxBoundary - xMin + xMinBoundary;
  old_ySize = yMax + yMaxBoundary - yMin + yMinBoundary;
  xPos = xMax + xMaxBoundary - 0.5 *old_xSize;
  yPos = yMax + yMaxBoundary - 0.5 *old_ySize;
  
  if(call_count%100==0){
    std::cout << " no. pixels in x and y: " << xSize << "," << ySize << "\n";
    std::cout << " (xMin, yMin): (" << xMin << ", " << yMin << ") and (xMax, yMax): (" << xMax << ", " << yMax << ")\n";
    std::cout << " MinBoundaries: (" << xMinBoundary << ", " << yMinBoundary << ") MaxBoundaries: (" << xMaxBoundary << ", " << yMaxBoundary << ")\n";
    std::cout << " oldSize: (" << old_xSize << ", " << old_ySize << ") Pos: (" << xPos << ", " << yPos << ")\n";
    std::cout << "\n";}
  
  call_count++;
  delete pixel;
}

  
void EUTelGeometricClusterImpl::getGeometricCenterOfGravity(float& xCoG, float& yCoG) const
{
  xCoG = 0;
  yCoG = 0;
    
  double totalCharge = 0;
    
  EUTelGeometricPixel* pixel = new EUTelGeometricPixel;
  for( unsigned int index = 0; index < size() ; index++ ) 
    {
      getSparsePixelAt( index , pixel);
	
      double curSignal = pixel->getSignal();
      xCoG += (pixel->getPosX())*curSignal;
      yCoG += (pixel->getPosY())*curSignal;
      totalCharge += curSignal;
    } 
    
  xCoG /= totalCharge;
  yCoG /= totalCharge;
  delete pixel;
}
  

