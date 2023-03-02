#ifndef PROBA_POINT_TYPE_H
#define PROBA_POINT_TYPE_H

#define PCL_NO_PRECOMPILE
#include <pcl/pcl_macros.h>
#include <pcl/point_types.h>
#include <pcl/point_cloud.h>
#include <pcl/io/pcd_io.h>

struct EIGEN_ALIGN16 _PointXYZP      // enforce SSE padding for correct memory alignment
{
  PCL_ADD_POINT4D;                  // preferred way of adding a XYZ+padding
  float proba ;
  //PCL_MAKE_ALIGNED_OPERATOR_NEW     // make sure our new allocators are aligned
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW     // make sure our new allocators are aligned
} ;                    
struct EIGEN_ALIGN16 PointXYZP : public _PointXYZP
{
    inline PointXYZP(const _PointXYZP &p) 
    {
        x = p.x;
        y = p.y;
        z = p.z;
        proba = p.proba;
    }
    
    inline PointXYZP()
    {
        x = 0;
        y = 0;
        z = 0;
        proba = 0;
    }
    
    inline PointXYZP(float _x, float _y, float _z, 
        float _proba)
    {
        x = _x ;
        y = _y;
        z = _z;
        proba = _proba;
    }
    //PCL_MAKE_ALIGNED_OPERATOR_NEW
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
} ;

POINT_CLOUD_REGISTER_POINT_STRUCT (PointXYZP,           
                                   (float, x, x)
                                   (float, y, y)
                                   (float, z, z)
                                   (float, proba, proba)
)

POINT_CLOUD_REGISTER_POINT_WRAPPER(PointXYZP, _PointXYZP)

#endif // PROBA_POINT_TYPE_H
