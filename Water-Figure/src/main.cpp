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
    //__________________________________________________________________________
        
    // Simulation settings
    
    bool surface=true;
    //bool surface=false;     
    GLfloat size=2.0;					//largeur de la grille dans l'espace
    float density=20.0;
    float viscosity=0.0;
    GLuint nbSamplesOneDir=80; 
    //GLuint nbSamplesX=nbSamplesOneDir; 
    GLuint nbSamplesX=120;
    GLuint nbSamplesY=120;
    // Proportion de l'écran
    //nbSamplesY=nbSamplesOneDir*application->height/application->width; 
    GLuint nbSamplesZ=nbSamplesOneDir;
    GLuint nbParticlesCoef=0; // On surface, nbSamples*4^nbParticlesCoefs particles //particles density
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

