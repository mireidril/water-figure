// Scene.hpp
// Template Fluids Simulation
// N. Dommanget dommange@univ-mlv.fr


#ifndef __SCENE_HPP__
#define __SCENE_HPP__

#include "GLHeaders.hpp"

// Forward declaration
class Camera;
class Object;

// A container for objects
class Scene
{
    public:
        GLuint maxStoredObjects;         // An initial limit of storable objects
        GLuint maxDrawnObjects;          // An initial limit of drawable objects
        
        Object ** storedObjects;         // Library of Objects to use from GPU
        GLuint nbStoredObjects;          // Number of stored objects
        
        GLuint * drawnObjects;           // Indices of objects from library to really draw
        GLfloat * drawnObjectsColors;    // Color for each drawn object
        GLfloat * drawnObjectsModels;    // Transformation matrix for each drawn object
        GLuint * drawnObjectsShaderIDs;  // ID of the shader to use for each drawn object  
        GLuint * drawnObjectsTexture0IDs; // ID of the texture to use on each drawn object for unit 0
        GLuint * drawnObjectsTexture1IDs; // ID of the texture to use on each drawn object for unit 1
        GLuint nbDrawnObjects;

        GLfloat defaultColor[4];         // Default color for drawn elements
        GLfloat defaultModel[16];        // Default transformation matrix for drawn elements
        GLuint defaultShaderID;          // Default shaderID for drawn elements        
        GLuint defaultTextureID;         // Default textureID for drawn elements
        
        GLfloat lightPosition[4];        // Position of the light used in shader
        GLfloat lightPower;              // Power of the light used in shader
        
        Camera * camera;                 // Camera used to watch the scene

        Scene();
        ~Scene();

        void init();

        GLuint storeObject(Object * object);
        GLuint addObjectToDraw(GLuint indexStoredObject);

        void setDrawnObjectColor(GLuint indexDrawnObject, GLfloat * color);
        void setDrawnObjectModel(GLuint indexDrawnObject, GLfloat * model);
        void setDrawnObjectShaderID(GLuint indexDrawnObject, GLuint shaderID);
        void setDrawnObjectTextureID(GLuint indexDrawnObject, GLuint textureUnit, GLuint textureID);    
        
        void setDefaultColor(GLfloat * defaultColor);
        void setDefaultModel(GLfloat * defaultModel);
        void setDefaultShaderID(GLuint defaultShaderID);
        void setDefaultTextureID(GLuint defaultTextureID);        
        
        void setLight(GLfloat * position, GLfloat power);

        void setAppearance(GLuint indexDrawnObject);
        void drawObjectsOfScene();
};

#endif // __SCENE_HPP__
