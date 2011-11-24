// Camera.cpp
// Template Fluids Simulation
// N. Dommanget dommange@univ-mlv.fr


#include "Camera.hpp"
#include "Tools.hpp"

#include <iostream>
#define _USE_MATH_DEFINES
#include <math.h>



// Default constructor
Camera::Camera()
{
    this->init();
}

// Cleans memory for Camera
Camera::~Camera()
{

}


// Inits parameters for the camera
void Camera::init()
{
    // View Data
    // Camera position and orientation  
    GLfloat c[]={0.0f, 0.0f, 1.0f}; // Camera position    
    GLfloat aim[]={0.0f, 0.0f, 0.0f}; // Where we look
    GLfloat up[]={0.0f, 1.0f, 0.0f}; // Vector pointing over the camera
    lookAt(c, aim, up);

    // Projection data
    // Projection type : perspective:true / orthographic:false
    this->perspectiveProjection=true;
    // Projection frustum data
    GLfloat l=1.0;
    this->left=-l; 
    this->right=l;
    this->bottom=-l;
    this->top=l; 
    this->Near=0.1;
    this->Far=100.0;
    setToIdentity(this->projection);
    updateProjection();
}


// Sores the data necassary to evaluate the porjection matrix
void Camera::setProjectionData(GLfloat left, GLfloat right, GLfloat bottom, GLfloat top, GLfloat n, GLfloat f)
{
    this->left=left;
    this->right=right;
    this->bottom=bottom;
    this->top=top;
    this->Near=n;
    this->Far=f;
}


// Turns camera projection from perspective to ortho or inverse
void Camera::switchCameraProjection()
{
    std::cout<<"Camera projection : ";
	if (this->perspectiveProjection)
    {
        std::cout<<"orthographic"<<std::endl;
		this->perspectiveProjection=false;
    }
	else
    {
        std::cout<<"perspective"<<std::endl;
		this->perspectiveProjection=true;
    }

    // Changes the matrix accordingly
    updateProjection();
}


// Builds the camera axis from position c, aim (focus) of the camera, and up vector
void Camera::lookAt(GLfloat * c, GLfloat * aim, GLfloat * up)
{
    for (GLuint iCoord=0 ; iCoord<3; iCoord++)
    {
        this->c[iCoord]=c[iCoord]; // c : camera position
        this->y[iCoord]=up[iCoord]; // y : up vector
        this->z[iCoord]=c[iCoord]-aim[iCoord]; // z : from aim to camera position
    }
    // Opposite direction of axis z
    GLfloat minusZ[]={-this->z[0], -this->z[1], -this->z[2]};
    // axis x  : from axis y and axis z
    vectorProduct(this->y, this->z, this->x);
    // axis y  : from new axis x and opposite of axis z
    vectorProduct(this->x, minusZ, this->y);
    // Normalizes all axis vectors
    normalize(this->x);
    normalize(this->y);
    normalize(this->z);

    // Builds the new view matrix
    updateView();
}


// Updates view
void Camera::updateView()
{
    // Rotation to be aligned with correct camera axis
    GLfloat RcInv[]={this->x[0], this->y[0], this->z[0], 0.0f,
                     this->x[1], this->y[1], this->z[1], 0.0f,
                     this->x[2], this->y[2], this->z[2], 0.0f,
                     0.0f,          0.0f,          0.0f,          1.0f};

    // Translation to be at the right distance from the scene
    GLfloat TcInv[]={1.0f,           0.0f,           0.0f,           0.0f,
                     0.0f,           1.0f,           0.0f,           0.0f,
                     0.0f,           0.0f,           1.0f,           0.0f,
                     -this->c[0], -this->c[1], -this->c[2], 1.0f};

    // Initializes
    setToIdentity(this->view);
    // Rotates
    multMatrixBtoMatrixA(this->view, RcInv);  
    // Translates
    multMatrixBtoMatrixA(this->view, TcInv); 
}

// Set the perspective like gluPerspective
void Camera::setPerspectiveFromAngle(GLfloat fovy, GLfloat aspectRatio)
{
    this->top = this->Near * tan(fovy/2.0f);
    this->bottom = - this->top;
    this->left = this->bottom * aspectRatio;
    this->right = this->top * aspectRatio;
    
    updateProjection();
}


// Updates the projection matrix from the data
void Camera::updateProjection()
{  
    GLfloat l=this->left;
    GLfloat r=this->right;
    GLfloat b=this->bottom;
    GLfloat t=this->top;
    GLfloat n=this->Near;
    GLfloat f=this->Far;

    if (this->perspectiveProjection) // Perspective projection
    {
        GLfloat P[]={ (2.0f*n)/(r-l), 0.0f,           0.0f,              0.0f,
                      0.0f,           (2.0f*n)/(t-b), 0.0f,              0.0f,
                      (r+l)/(r-l),   (t+b)/(t-b),   -(f+n)/(f-n),    -1.0f,
                      0.0f,           0.0f,           -(2.0f*f*n)/(f-n), 0.0f}; // Perspective projection
        for (int iMatrixCoord=0 ; iMatrixCoord<16 ; iMatrixCoord++)
            this->projection[iMatrixCoord]=P[iMatrixCoord];
    }
    else // Orthographic projection
    { 
        GLfloat P[]={ 2.0f/(r-l),   0.0f,         0.0f,         0.0f,
                      0.0f,         2.0f/(t-b),   0.0f,         0.0f,
                      0.0f,         0.0f,         -2.0f/(f-n),  0.0f,
                    -(r+l)/(r-l), -(t+b)/(t-b), -(f+n)/(f-n), 1.0f}; // Orthographic projection
        for (int iMatrixCoord=0 ; iMatrixCoord<16 ; iMatrixCoord++)
            this->projection[iMatrixCoord]=P[iMatrixCoord];
    }
}

