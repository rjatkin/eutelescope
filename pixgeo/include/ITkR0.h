#ifndef ITKR0_H
#define	ITKR0_H

  /** @class ITkR0
	* This class is the implementation of  @class EUTelGenericPixGeoDescr
	* for the ITkR0 which is the ATLAS ITk strip end-cap sensor closest to
	* the beam pipe.
	* The geoemtry is as following: Four rows of strips, each a different length.
	* The two inner rows have 1026 strips while the outer two rows have 1054 strips.
	* The strips are defined as annuli and focus on a point not at the origin
	* while they are arranged in a curve centered on the origin, due to the stereo angle. 
    */

//STL
#include <string> //std::string
#include <utility> //std::pair

//EUTELESCOPE
#include "EUTelGenericPixGeoDescr.h"

//ROOT
#include "TGeoMaterial.h"
#include "TGeoMedium.h"
#include "TGeoVolume.h"

namespace eutelescope {
namespace geo {

class ITkR0 : public EUTelGenericPixGeoDescr {
	
	public:
		ITkR0();
		~ITkR0();

		void createRootDescr(char const *);
		std::string getPixName(int, int);
		std::pair<int, int> getPixIndex(char const *);

	protected:
		TGeoMaterial* matSi;
		TGeoMedium* Si;
		TGeoVolume* plane,*row1strip,*row2strip,*row3strip,*row4strip;
		Double_t ephi,sphi;
		Double_t pi=3.14159265358979,deg=180/pi;
		Double_t dz=0.155,dphi=0.198306302808;
		Double_t l1=18.981,l2=23.981,l3=28.98,l4=31.981;
		Double_t r0=384.5,r1=r0+l1,r2=r1+l2,r3=r2+l3,r4=488.423;
		Double_t stereo=0.02, R=438.614;
		Int_t n_phi1=1026,n_phi2=1154;
		Double_t dphi1=0.0001932745,dphi2=0.0001718368;
		Double_t phi_i,b,r,c,x1,x,y1,y,gradient,theta;
		//corners of sensor
		Double_t Ax=-39.565,Ay=-56.608,Bx=-47.856,By=47.920,Cx=49.943,Cy=47.711,Dx=37.479,Dy=-56.398;
		//focus point of strips
		Double_t Fx=-8.771695,Fy=0.0877198;
		
};

extern "C"
{
	EUTelGenericPixGeoDescr* maker();
}

} //namespace geo
} //namespace eutelescope

#endif	//ITKR0_H
