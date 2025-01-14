#ifndef GLOBAL_H
#define GLOBAL_H
#define GL_GLEXT_PROTOTYPES
#define GLM_FORCE_RADIANS

#define QT_NO_OPENGL_ES_2
#include <QGLWidget>
#include <QGLFunctions>

#include "GL/gl.h"
#include "GL/glu.h"
#include "GL/glext.h"
#include "glm/glm.hpp"

#include <glog/logging.h>

#include <iostream>
#include<numeric>
#include <vector>
#include <cmath>
#include <math.h>
#include <algorithm>
#include <fstream>
#include <iostream>
#define MESH_AMOUNT 50 //amount of
#define CAUSTIC_DOMAIN 0.5

/*Physical model*/
#define MATERIAL_REFRACTIV_INDEX 1.5//value for acrylic used in the paper to do test experimently

/*ceres solver settings*/
#define CONVERGENCE_LIMIT 9e-2001
#define AIR_REFRACTIV_INDEX 1
#define EBAR_DETH 39
#define EINT_WEIGHT 1.0
#define EBAR_WEIGHT 1.0
#define EDIR_WEIGHT 0.6
#define EREG_WEIGHT 10.0
#define MAX_Y 1
#define MAX_Z 1
#define RAY_SHOOT 50

// if set to true, it will solve for the reflective caustics instead of refractieve
#define REFLECTIVE_CAUSTICS false 

#if (REFLECTIVE_CAUSTICS == false)
    #define INCIDENT_RAY_X -1.0f 
    #define INCIDENT_RAY_Y 0.0f
    #define INCIDENT_RAY_Z 0.0f
#else
    #define INCIDENT_RAY_X 1.0f 
    #define INCIDENT_RAY_Y 0.0f
    #define INCIDENT_RAY_Z 0.0f
#endif

//typedef float (n_array)[961]; // size of vertex in the front face 961 without the edges 1089 with edges




#endif // GLOBAL_H
