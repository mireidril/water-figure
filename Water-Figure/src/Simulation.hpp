// Simulation.hpp
// Template Fluids Simulation
// N. Dommanget dommange@univ-mlv.fr


#ifndef __SIMULATION_HPP__
#define __SIMULATION_HPP__

#include "GLHeaders.hpp"

#include <SDL.h>

#include <vector>
#include <algorithm>

// Mersenne Twister random generation
#include "../api/random/MersenneTwister.h"
#include "../api/eigen-eigen-3.0.3/Eigen/Dense"
#include "../api/TinyXML/tinyxml.h"

// Enum of types of cells
enum { FLUID, AIR, SOLID };

// Forward declaration
class Scene;
class Object;

// A container for objects
class Simulation
{
     public:
        GLuint ppmPhase;						// Which phase of the PPM are we ?
     
        // About the data stored on centers of the MAC grid
        GLuint nbSamplesX;       // Samples number on X axis
        GLuint nbSamplesY;       // Samples number on Y axis
        GLuint nbSamplesZ;       // Samples number on Z axis
        GLuint nbSamples;        // Total samples number
        GLfloat * samples;       // Samples positions
        GLfloat * pressures;     // Samples pressures
        GLfloat * colors;        // Samples colors
		GLfloat * forces;		 // Samples forces
		GLuint * types;			 // Types (0 : cell is fluid, 1 : cell is empty, 2 : cell is solid)
           
        // About the data stored on borders of the MAC grid
        GLuint nbBordersX;       // Number of left/right cell borders
        GLuint nbBordersY;       // Number of bottom/top cell borders
        GLuint nbBordersZ;       // Number of back/front cell borders
        GLfloat * velocitiesX;   // X velocity components, sampled on left/right borders
        GLfloat * velocitiesY;   // Y velocity components, sampled on bottom/top borders
        GLfloat * velocitiesZ;   // Z velocity components, sampled on back/front borders
  
        // About the particles used for visualization (advected on the whole grid)
        GLuint nbParticlesCoef;   // Coefficient for density of particles per cell
        GLuint nbParticles;       // Number of rendered particles          
                                  // (=nbSamples*[4,8]^nbParticlesCoef)
        GLfloat * particles;      // Particle positions
        GLfloat * particleColors; // Grid interpolated particles colors
        GLfloat * particleVelocities; // Grid interpolated particles velocities
        
        // Geometric parameters
        bool surface;             // True if 2D simulation
        GLfloat size;             // Width of the grid
        GLfloat h;                // Cell side length (square or cubes)
        GLfloat offset[3];        // Start point left/bottom/back of the grid
        GLfloat offsetInCell;     // Offset from left/bottom/back of grid cell
        bool solidWalls;          // True if grid boundary is expected to be solid
        GLfloat vectorScale;  // Scale coefficient to display vector lengths
        
        // Simulation parameters
        GLfloat density;          // Expected density for the fluid
        GLfloat viscosity;        // Expected viscosity for the fluid
        GLfloat g[3];             // Vector of acceleration due to gravity
        
        // Additional required data
        Scene * scene;            // Scene in which to draw the simulation
        GLuint defaultShaderID;   // Shader used for test rendering
        GLuint spriteShaderID;    // Shader used for particles
        MTRand * rand;            // Mersenne-Twister random generator
        
        // Objects for OpenGL visualization
        Object * objectSamples;    // For rendering samples (position/color checking)
        Object * objectVelocitiesCenters; // For rendering velocities interpolated on samples
        Object * objectVelocitiesBorders; // For rendering velocities components
        Object * objectParticles;  // For rendering particles (position/color checking)
        Object * objectParticleVelocities; // For rendering particles interpolated velocities
		Object * objectPressures; //For rendering pressures (pressure checking)
		Object * objectForces; //For rendering pressures (pressure checking)
        
        Simulation(Scene * scene, bool surface, GLfloat size, GLfloat density, GLfloat viscosity, GLuint nbSamplesX, GLuint nbSamplesY, GLuint nbSamplesZ, GLuint nbParticlesCoef, GLuint defaultShaderID, GLuint spriteShaderID, bool solidWalls);
        ~Simulation();

		void getCell(GLfloat * pos, int * iX, int * iY, int * iZ);

        void initSimulation();
        void initVisualization();
        
        void initSamples();
        void buildSamples(Object * object);
        void drawSamples();
        
        void buildGrid(Object * object);
        void drawGrid();
        
        void initVelocitiesBorders();
        void buildVelocitiesBorders(Object * object);
		void enforceVelocitiesBorders();
        void drawVelocitiesBorders();
        
        void updateVelocitiesFromCenters(GLfloat * centeredVel);
        void updateVelocitiesFromBorders();
        
        void buildVelocitiesCenters(Object * object);
        void drawVelocitiesCenters();
        
        void initParticles();
        //void addParticles(GLuint & iSample);
        void buildParticles(Object * object);
        void drawParticles();
        
        void buildParticleVelocities(Object * object);
        void drawParticlesVelocities();

		void buildForces(Object * object);
		void drawForces();
		void applyForces(GLfloat dt);

		void buildPressures(Object * object);
		void drawPressures();
		void setDivergences(double * divergences);

		GLuint type(int iX, int iY, int iZ);
		void buildTypes(Object* object);
		void drawTypes();
        
        void interpolateFromCenters(GLfloat * data, GLfloat * position, GLfloat * result);
        void interpolateFromBorders(GLfloat * dataX, GLfloat * dataY, GLfloat * dataZ, GLfloat * position, GLfloat * result);
        void forwardEuler1stOdrTimeIntegration(GLfloat direction, GLfloat * knownValue, GLfloat * knownVariation, GLfloat dt, GLfloat * resultValue);
        void rungeKutta2ndOdrTimeIntegration(GLfloat direction, GLfloat * knownValue, GLfloat * knownVariation, GLfloat dt, GLfloat * resultValue);
        void rungeKutta3rdOdrTimeIntegration(GLfloat direction, GLfloat * knownValue, GLfloat * knownVariation, GLfloat dt, GLfloat * resultValue);
        void integrate(GLfloat direction, GLfloat * knownValue, GLfloat * knownVariation, GLfloat dt, GLfloat * resultValue);
        
        void advectParticles(GLfloat dt); 
        void advectColors(GLfloat dt);
        void advectVelocities(GLfloat dt);
		void project(GLfloat dt);
		
		void multiplySparseMatrix(double * Adiag, double * Aright, double * Atop, Eigen::VectorXd v, Eigen::VectorXd * result);
		void conjugateGradient(double * Adiag, double * Aright, double * Atop, double * b);
		void MICPreconditioner(double * Adiag, double * Aright, double * Atop, Eigen::VectorXd* precon);
		void applyPreconditioner(double * Aright, double * Atop, Eigen::VectorXd precon, Eigen::VectorXd r, Eigen::VectorXd * z);
		   
		void threeColorImageHandler(unsigned char *image);
		//All the samples which are initially FLUID
		//std::vector<GLuint> infiniteFluidsCase;
		//void isFluidCaseEmpty();
		
		void loadForcesFrom(const char* pFilename);
		void parse( TiXmlNode* pParent, unsigned int indent);
		void applyForcesLoaded( int x, int y, int x2, int y2, double valueX, double valueY );
		
		void resetForces();
		
        void update();
        
        void render();
};

#endif // __SIMULATION_HPP__

