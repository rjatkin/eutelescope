#include "MimosaR02.h"

namespace eutelescope {
  namespace geo {
    
    std::string GetPixName(int x, int y);
    
    //EUTelGenericPixGeoDescr still needs to be altered to take proper values for the R0 sensor dimensions. Strips are unaffected though
    MimosaR02::MimosaR02(): EUTelGenericPixGeoDescr(97.769, 106.417, 0.00031,//size X, Y, Z (size in mm of sensor) (taken using points Bx & Cx, Ay & By)
						  0, 1153, 0, 3,	//minX maxX minY maxY number of pixels in x and y  
						  9.3660734 )		//rad length
    {
      
      //Create the material for the sensor
      matSi = new TGeoMaterial( "Si", 28.0855 , 14.0, 2.33, -_radLength, 45.753206 );
      Si = new TGeoMedium("MimosaSilicon",1, matSi);

      plane = _tGeoManager->MakeBox("sensarea_box",Si,150,500,1);
      //Define the four strip volumes
      row1strip = _tGeoManager->MakeTubs( "sensarea_row1" , Si, r0, r1, dz, 90+(-dphi1/2)*deg,90+(dphi1/2)*deg);
      row2strip = _tGeoManager->MakeTubs( "sensarea_row2" , Si, r1, r2, dz, 90+(-dphi1/2)*deg,90+(dphi1/2)*deg);
      row3strip = _tGeoManager->MakeTubs( "sensarea_row3" , Si, r2, r3, dz, 90+(-dphi2/2)*deg,90+(dphi2/2)*deg);
      row4strip = _tGeoManager->MakeTubs( "sensarea_row4" , Si, r3, r4, dz, 90+(-dphi2/2)*deg,90+(dphi2/2)*deg);

      //The following formula is defined in the ATLAS12EC Technical Specs
      //get angle of first strips in first two rows
      phi_i=(-float(n_phi1)/2+0.5)*dphi1;
      b=-2*(2*R*sin(stereo/2))*sin(stereo/2+phi_i);    
      c=pow((2*R*sin(stereo/2)),2)-pow(r0,2);
      r=0.5*(-b+sqrt(pow(b,2)-4*c));
      x1=r*cos(phi_i+stereo) - R*cos(stereo);//x and y are defined as if sensor centre is at the origin and strips are in the x direction
      y1=r*sin(phi_i+stereo) - R*sin(stereo);
      c=pow((2*R*sin(stereo/2)),2)-pow(r1,2);
      r=0.5*(-b+sqrt(pow(b,2)-4*c));
      x =r*cos(phi_i+stereo) - R*cos(stereo);
      y =r*sin(phi_i+stereo) - R*sin(stereo);
      gradient=-(x-x1)/(y-y1);
      theta=atan(gradient)-dphi1;
  
      //placement of first two rows
      for( int i = 0; i < n_phi1; i++ ){
	TGeoCombiTrans* transform=new TGeoCombiTrans(0,0,0,new TGeoRotation("rot",0,0,0));
	TGeoCombiTrans* transform1=new TGeoCombiTrans(0,0,0,new TGeoRotation("rot",0,0,0));
	//get angle of each strip. will be same for both rows
	theta+=dphi1;
	//get position of each strip for first row
	phi_i=(i-float(n_phi1)/2+0.5)*dphi1;
	b=-2*(2*R*sin(stereo/2))*sin(stereo/2+phi_i);    
	c=pow((2*R*sin(stereo/2)),2)-pow(r0,2);
	r=0.5*(-b+sqrt(pow(b,2)-4*c));
	x=r*cos(phi_i+stereo) - R*cos(stereo);
	y=r*sin(phi_i+stereo) - R*sin(stereo);
	//create first transform
	//rotate to get correct angle of strip
	//position is rotated by 90 degrees to get it on y axis, hence (-y,x)
	transform->RotateZ(theta*deg-90);
	transform->SetTranslation(-y-r0*cos(theta),x-r0*sin(theta),0);
	//get each position of the strips for second row 
	c=pow((2*R*sin(stereo/2)),2)-pow(r1,2);
	r=0.5*(-b+sqrt(pow(b,2)-4*c));
	x=r*cos(phi_i+stereo) - R*cos(stereo);
	y=r*sin(phi_i+stereo) - R*sin(stereo);
	//create second transform
	transform1->RotateZ(theta*deg-90);
	transform1->SetTranslation(-y-r1*cos(theta),x-r1*sin(theta),0);
	//add the nodes
	plane->AddNode(row1strip,i+1,transform);
	plane->AddNode(row2strip,i+1,transform1);
      }
  
      //get angle of first strips in outer two rows
      phi_i=(-float(n_phi2)/2+0.5)*dphi2;
      b=-2*(2*R*sin(stereo/2))*sin(stereo/2+phi_i);    
      c=pow((2*R*sin(stereo/2)),2)-pow(r2,2);
      r=0.5*(-b+sqrt(pow(b,2)-4*c));
      x1=r*cos(phi_i+stereo) - R*cos(stereo);
      y1=r*sin(phi_i+stereo) - R*sin(stereo);
      c=pow((2*R*sin(stereo/2)),2)-pow(r3,2);
      r=0.5*(-b+sqrt(pow(b,2)-4*c));
      x =r*cos(phi_i+stereo) - R*cos(stereo);
      y =r*sin(phi_i+stereo) - R*sin(stereo);
      gradient=-(x-x1)/(y-y1);
      theta=atan(gradient)-dphi2;
  
      //placement of second two rows
      for( int i = 0; i < n_phi2; i++ ){
	TGeoCombiTrans* transform=new TGeoCombiTrans(0,0,0,new TGeoRotation("rot",0,0,0));
	TGeoCombiTrans* transform1=new TGeoCombiTrans(0,0,0,new TGeoRotation("rot",0,0,0));
	//get angle of each strip. will be same for both rows
	theta+=dphi2;
	//get each position of the strips for first row
	phi_i=(i-float(n_phi2)/2+0.5)*dphi2;
	b=-2*(2*R*sin(stereo/2))*sin(stereo/2+phi_i);    
	c=pow((2*R*sin(stereo/2)),2)-pow(r2,2);
	r=0.5*(-b+sqrt(pow(b,2)-4*c));
	x=r*cos(phi_i+stereo) - R*cos(stereo);
	y=r*sin(phi_i+stereo) - R*sin(stereo);
	//create first transform
	transform->RotateZ(theta*deg-90);
	transform->SetTranslation(-y-r2*cos(theta),x-r2*sin(theta),0);
	//get each position of the strips for second row 
	c=pow((2*R*sin(stereo/2)),2)-pow(r3,2);
	r=0.5*(-b+sqrt(pow(b,2)-4*c));
	x=r*cos(phi_i+stereo) - R*cos(stereo);
	y=r*sin(phi_i+stereo) - R*sin(stereo);
	//create second transform
	transform1->RotateZ(theta*deg-90);
	transform1->SetTranslation(-y-r3*cos(theta),x-r3*sin(theta),0);
	//add the nodes
	plane->AddNode(row3strip,i+1,transform);
	plane->AddNode(row4strip,i+1,transform1);//transform);
      }
  
    }
  
    MimosaR02::~ MimosaR02()
    {
      //It appears that ROOT will take ownership and delete that stuff! 
      //delete matSi,
      //delete Si;
    }
  
    void  MimosaR02::createRootDescr(char const * planeVolume)
    {
      //Get the plane as provided by the EUTelGeometryTelescopeGeoDescription
      TGeoVolume* topplane =_tGeoManager->GetVolume(planeVolume);
      //Add the sensitive area to the plane
      topplane->AddNode(plane, 1,new TGeoTranslation(0,0,0) );
       
    }
  
    std::string MimosaR02::getPixName(int x , int y){
      char buffer [100];
      //return path to the pixel, don't forget to shift indices by +1+
      snprintf( buffer, 100, "/sensarea_box_1/sensarea_row%d_%d",y+1,x+1);
    
      return std::string( buffer ); 
    }
  
  
  
    /*TODO*/ std::pair<int, int>  MimosaR02::getPixIndex(char const*){return std::make_pair(0,0); }
  
    EUTelGenericPixGeoDescr* maker()
    {
      MimosaR02* mPixGeoDescr = new MimosaR02();
      return dynamic_cast<EUTelGenericPixGeoDescr*>(mPixGeoDescr);
    }
  
  } //namespace geo
} //namespace eutelescope

