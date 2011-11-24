// Camera.hpp
// Template Fluids Simulation
// N. Dommanget dommange@univ-mlv.fr


#ifndef __CAMERA_HPP__
#define __CAMERA_HPP__

#include "GLHeaders.hpp"

// A camera to watch the scene
class Camera
{
    public:

        // View Data
        GLfloat c[3];               // Camera position 
        GLfloat x[3];               // Camera axis x : right side
        GLfloat y[3];               // Camera axis y : up
        GLfloat z[3];               // Camera axis z : backward
        GLfloat view[16];           // View matrix

        // Projection data
        bool perspectiveProjection; // persepective projection:true / orthographic projection:false
        GLfloat left;               // x coord from center to left plane of frustum
        GLfloat right;              // x coord from center to right plane of frustum
        GLfloat bottom;             // y coord from center to bottom plane of frustum
        GLfloat top;                // y coord from center to top plane of frustum
        GLfloat Near;               // z coord from c to near plane of frustum
        GLfloat Far;                // z coord from c to far plane of frustum
        GLfloat projection[16];     // Projection matrix


        Camera();
        ~Camera();
        
        void init();

        void setProjectionData(GLfloat left, GLfloat right, GLfloat bottom, GLfloat top, GLfloat n, GLfloat f);
        void switchCameraProjection();

        void lookAt(GLfloat * c, GLfloat * aim, GLfloat * up);
        void setPerspectiveFromAngle(GLfloat fovy, GLfloat aspectRatio);

        void updateView();
        void updateProjection();
};

#endif // __CAMERA_HPP__
