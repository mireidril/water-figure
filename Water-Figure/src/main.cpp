// main.c
// Template Fluids Simulation
// N. Dommanget dommange@univ-mlv.fr


#include "GLHeaders.hpp"
#include "Application.hpp"
#include "Builders.hpp"
#include "Tools.hpp"
#include "Scene.hpp"
#include "Simulation.hpp"
#include "Camera.hpp"
#include "Object.hpp"

#include <iostream>
#include <vector>
#include <string>

// Entry point in the program
int main(int argc, char **argv)
{
    std::cout<<std::endl<<std::endl<<"________Template Fluids Simulation________"<<std::endl<<std::endl;  
    //__________________________________________________________________________ 
    
    // Application creation
    //std::cout<<"  Application creation"<<std::endl;
    Application * application=new Application();
    //__________________________________________________________________________ 

	// Elements creation
    //std::cout<<"    Scene creation"<<std::endl;
    Scene * scene=new Scene();
    application->setScene(scene);
    //__________________________________________________________________________

    // Camera position and orientation
    //std::cout<<"    Camera settings"<<std::endl;
    GLfloat c[]={0.0, 0.0, 2.0};
	GLfloat aim[]={0.0, 0.0, 0.0};
	GLfloat up[]={0.0, 1.0, 0.0};
    scene->camera->lookAt(c, aim, up);       
    //__________________________________________________________________________

    // Creation of the shaders
    //std::cout<<"    Shaders reading, compiling and linking"<<std::endl;
    
    // Compilation and storage of the program
    std::vector<std::string> files;

    // The special shaders used for mac os are only necesary until Snow Leopard
    // From Leon on, mac os drivers version for Opengl is 3.2
    // Therefore the shaders in folder oldGLSL become out-of-date
	#ifdef __APPLE__
        files.push_back("../shaders/oldGLSL/shaderTools.glsl");
    #else
        files.push_back("../shaders/shaderTools.glsl");
    #endif
    
    // One regular Phong shader
	#ifdef __APPLE__
	    files.push_back("../shaders/oldGLSL/lightingShader.glsl");
    #else
	    files.push_back("../shaders/lightingShader.glsl");
    #endif
    GLuint lightingShaderID=loadProgram(files); files.pop_back();
    scene->setDefaultShaderID(lightingShaderID);
    // Uniforms filling
    glUseProgram(lightingShaderID);
        // Material creation and to shader
    GLfloat ambient[]={1.0f, 1.0f, 1.0f, 1.0f}; GLfloat ka=0.01f;
    GLfloat diffuse[]={1.0f, 1.0f, 1.0f, 1.0f}; GLfloat kd=1.0f;
    GLfloat specular[]={1.0f, 1.0f, 1.0f, 1.0f}; GLfloat ks=2.0f; GLfloat shininess=5.0f;
    setMaterialInShader(lightingShaderID, ambient, diffuse, specular, ka, kd, ks, shininess);
    
    // One texture and Phong shader
	#ifdef __APPLE__
	    files.push_back("shaders/oldGLSL/lightingTexturingShader.glsl");
    #else
	    files.push_back("../shaders/lightingTexturingShader.glsl");
    #endif
    GLuint lightingTextureShaderID=loadProgram(files); files.pop_back();
    glUseProgram(lightingTextureShaderID);
    setMaterialInShader(lightingTextureShaderID, ambient, diffuse, specular, ka, kd, ks, shininess);
    setTextureUnitsInShader(lightingTextureShaderID);
    
    // One sprite texture and Phong shader
	#ifdef __APPLE__
	    files.push_back("shaders/oldGLSL/spriteShader.glsl");
    #else
	    files.push_back("../shaders/spriteShader.glsl");
    #endif
    GLuint spriteShaderID=loadProgram(files); files.pop_back();
    glUseProgram(spriteShaderID);
    setMaterialInShader(spriteShaderID, ambient, diffuse, specular, ka, kd, ks, shininess);
    setTextureUnitsInShader(spriteShaderID);
    
    //__________________________________________________________________________ 

    // Creation of the objects to store
    //std::cout<<"    Generic objects creation and storage :"<<std::endl;
    Object * objectL0=new Object(GL_LINES);     GLuint storedObjectL0ID=scene->storeObject(objectL0);
    /*Object * objectT0=new Object(GL_TRIANGLES); GLuint storedObjectT0ID=scene->storeObject(objectT0); 
    Object * objectT1=new Object(GL_TRIANGLES); GLuint storedObjectT1ID=scene->storeObject(objectT1);
    Object * objectT2=new Object(GL_TRIANGLES); GLuint storedObjectT2ID=scene->storeObject(objectT2);
    Object * objectT3=new Object(GL_TRIANGLES); GLuint storedObjectT3ID=scene->storeObject(objectT3);*/
    //__________________________________________________________________________
    
    // Object building
    //std::cout<<"    Generic objects building :"<<std::endl;
    // Axis
    //buildAxis(objectL0);
    
	// A procedural plane
    /*buildEarthPlane(objectT0, 200.0, 200.0);
    
    // A procedural sphere
    GLfloat radius=0.2; GLint discLong=20; GLint discLat=10;
    buildSphere_TrSmoothNonRed(objectT1, radius, discLat, discLong);
    
    // An obj monkey head (smooth edges)
    bool smoothObjectFlag=true;
    std::string fileName="../objs/monkeyHead.obj";
    buildObjectGeometryFromOBJ(objectT2, fileName, smoothObjectFlag);
    
    // An obj textured house  (hard edges)
    smoothObjectFlag=false;
    fileName="../objs/house.obj";
    buildObjectGeometryFromOBJ(objectT3, fileName, smoothObjectFlag); */
    //__________________________________________________________________________
              
    // Objects we want to see
    //std::cout<<"    Objects to draw setting"<<std::endl;
    
    // Axis
    //GLuint axisID=scene->addObjectToDraw(storedObjectL0ID); 
       
    /*// One plane (grey by default)
    GLuint planeID=scene->addObjectToDraw(storedObjectT0ID);
    GLfloat grass[]={0.35, 0.74, 0.25, 1.0};
    scene->setDrawnObjectColor(planeID, grass);
    
    // One yellow sphere
    GLuint sphereID=scene->addObjectToDraw(storedObjectT1ID);
    GLfloat T1[16]; GLfloat t1[3]={-0.5, 1.0, -0.2}; setToTranslate(T1, t1);
    scene->setDrawnObjectModel(sphereID, T1);
    GLfloat yellow[]={1.0, 1.0, 0.0, 1.0};
    scene->setDrawnObjectColor(sphereID, yellow);
    
    // One Orange monkey head
    GLuint monkeyHeadID=scene->addObjectToDraw(storedObjectT2ID);
    GLfloat orange[]={0.7, 0.4, 0.1, 1.0};
    GLfloat T2[16]; GLfloat t2[3]={0.5, 0.5, 0.5}; setToTranslate(T2, t2);
    GLfloat S2[16]; GLfloat s2[3]={0.1, 0.1, 0.1}; setToScale(S2, s2);
    multMatrixBtoMatrixA(T2, S2);
    scene->setDrawnObjectModel(monkeyHeadID, T2);
    scene->setDrawnObjectColor(monkeyHeadID, orange);
    
    // A field of textured houses
    GLuint houseTextureDiffuseID=loadTexture("../textures/house_diffuse.ppm");
    GLuint houseTextureSpecularID=loadTexture("../textures/house_spec.ppm");
    GLuint fieldObjectID; GLfloat T3[16]; GLfloat t3[3];
    GLint disc=5; GLfloat step=4.0/GLfloat(disc); GLfloat x=-2.0; GLfloat z=-2.0;
    for (int iX=-disc ; iX<=disc ; iX++)
    {
        x=iX*step;
        for (int iZ=-disc ; iZ<=disc ; iZ++)
        {
            z=iZ*step;
            fieldObjectID=scene->addObjectToDraw(storedObjectT3ID);
            t3[0]=iX; t3[1]=0.5; t3[2]=iZ; setToTranslate(T3, t3);
            scene->setDrawnObjectModel(fieldObjectID, T3);
            scene->setDrawnObjectShaderID(fieldObjectID, phongTextureShaderID);
            scene->setDrawnObjectTextureID(fieldObjectID, 0, houseTextureDiffuseID);
            scene->setDrawnObjectTextureID(fieldObjectID, 1, houseTextureSpecularID);
        }
    }*/
    //__________________________________________________________________________
    
    // Other informations
    // Background color creation
    //GLfloat skyColor[]={0.57, 0.80, 0.93, 1.0};
    GLfloat black[]={0.0, 0.0, 0.0, 1.0};
    application->setBackgroundColor(black);
    // Light creation
    GLfloat lightPosition[]={0.0, 5.0, 0.0, 1.0}; GLfloat lightPower=1.0;
    scene->setLight(lightPosition, lightPower);
    
    //__________________________________________________________________________
        
    // Simulation settings
    
    bool surface=true;
    //bool surface=false;     
    GLfloat size=2.0;					//largeur de la grille dans l'espace
    float density=1000.0;
    float viscosity=0.0;
    GLuint nbSamplesOneDir=10; 
    GLuint nbSamplesX=nbSamplesOneDir; 
    GLuint nbSamplesY=nbSamplesOneDir;
    GLuint nbSamplesZ=nbSamplesOneDir;
    GLuint nbParticlesCoef=2; // On surface, nbSamples*4^nbParticlesCoefs particles //particles density
                              // In volume,  nbSamples*6^nbParticlesCoefs particles
    bool solidWalls=true;
    //bool solidWalls=false;
  
    application->simulation=new Simulation(scene,
                                           surface,
                                           size,
                                           density,
                                           viscosity,
                                           nbSamplesX,
                                           nbSamplesY,
                                           nbSamplesZ,
                                           nbParticlesCoef,
                                           lightingShaderID,
                                           spriteShaderID,
                                           solidWalls);
    
	//application->simulation->drawSamples();					// Centre des carrés
    //application->simulation->drawGrid();						// Grille
    //application->simulation->drawVelocitiesBorders();			// 4 bords des carrés
    //application->simulation->drawVelocitiesCenters();			// Vitesse des points centraux
    application->simulation->drawParticles();					// Particules
    //application->simulation->drawParticlesVelocities();		// Vitesse des particules
	//application->simulation->drawForces();					// Forces appliquées
	//application->simulation->drawPressures();					// Pressions par carrés
	//application->simulation->drawTypes();						// Types des carrés
    
    
    //__________________________________________________________________________
    
    // Errors checking
    //std::cout<<"    Errors cheking before the loop"<<std::endl;
    printGlErrors(); std::cout<<std::endl;
    //__________________________________________________________________________  
      
    // Loop start !
    //std::cout<<"    Looping !! "<<std::endl;
    application->lastStartTime=getTime();
    application->initTimers();
    application->loop();
    printGlErrors(); std::cout<<std::endl;
    //__________________________________________________________________________

    // Cleaning and finishing
    //std::cout<<"    Cleaning before leaving"<<std::endl; 
    delete application;
   
	return 0;
}

