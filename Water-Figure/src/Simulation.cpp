// Simulation.cpp
// Template Fluids Simulation
// N. Dommanget dommange@univ-mlv.fr

#include "Simulation.hpp"
#include "Scene.hpp"
#include "Tools.hpp"
#include "Camera.hpp"
#include "Object.hpp"

#include <iostream>
#define _USE_MATH_DEFINES
#include <math.h>



//-------------------------------------
//--------PPM interpretation-----------
//-------------------------------------
///Types : 	black : solid
///			white : air
///			color : fluid of this color

void Simulation::threeColorImageHandler(unsigned char *image)
{
	std::cout<<"threeColorImageHandler nbSamplesX:"<<nbSamplesX<<" nbSamplesY"<<nbSamplesY<<std::endl;
	for (GLuint iY=0 ; iY<nbSamplesY ; iY++)
	{
		for (GLuint iX=0 ; iX<nbSamplesX ; iX++)
		{
            GLuint iSample = (iY*(nbSamplesX)+ iX);
            GLuint iImage = (iY*(nbSamplesX)+ iX)*3;
            if((image[iImage] == 0)&& (image[iImage+1] == 0) && (image[iImage+2] == 0) ){
                types[iSample] = SOLID;
                //std::cout<<"SOLID";
            }
            else if((image[iImage] == 255)&& (image[iImage+1] == 255) && (image[iImage+2] == 255) ){
                types[iSample] = AIR;
                //std::cout<<"AIR";
            }
            else {
                types[iSample] = FLUID;
                colors[iSample*4+0]= int(image[iImage]);
                colors[iSample*4+1]= int(image[iImage+1]);
                colors[iSample*4+2]= int(image[iImage+2]);
                colors[iSample*4+3]=1.0f;
				//std::cout<<"FLUID "<<int(image[iImage])<<"  ,  "<<int(image[iImage+1])<<"  ,  "<<int(image[iImage+2])<<std::endl; 
				//TODO this actually init the CELL's COLOR. 
				//To colorate the particles, use interpolateFromCenters(colors, &(particles[iParticles*4+0]), &(particleColors[iParticles*4+0])); as in initParticles()
            }
            //std::cout<<"  X:"<<iX<<" Y:"<<iY<<" iSample: "<<iSample<<"  iImage :"<<iImage<<std::endl;
        }
    }
}

//-------------------------------------
//--------XML interpretation-----------
//-------------------------------------
// load the named file and parse it
void Simulation::loadForcesFrom(const char* pFilename)
{
	TiXmlDocument doc(pFilename);
	bool loadOkay = doc.LoadFile();
	if (loadOkay)
	{
		printf("\n%s:\n", pFilename);
		//Before parsing, reset forces
		resetForces();
		parse( &doc, 0 );
	}
	else
	{
		printf("Failed to load file \"%s\"\n", pFilename);
	}
}

//Parse : for each force, parse the attributes, then call applyForcesLoaded.
void Simulation::parse( TiXmlNode* pParent, unsigned int indent = 0 )
{
	if ( !pParent ) return;
	
	TiXmlElement* child = pParent->FirstChild( "forces" )->FirstChild( "force" )->ToElement();
	if( child )
		for( ; child; child=child->NextSiblingElement() )
		{
			int x, y , x2, y2;
			double valueX, valueY;
			child->QueryIntAttribute("x", &x);
			child->QueryIntAttribute("y", &y);
			child->QueryIntAttribute("x2", &x2);
			child->QueryIntAttribute("y2", &y2);
			child->QueryDoubleAttribute("valueX", &valueX);
			child->QueryDoubleAttribute("valueY", &valueY);
			
			applyForcesLoaded( x, y, x2, y2, valueX, valueY );
		}
	
}

// Apply the force given.
// The points (x,y) and (x2,y2) form a rectangle where the force (valueX, valueY) is going to be applied.
void Simulation::applyForcesLoaded( int x, int y, int x2, int y2, double valueX, double valueY )
{
	printf( "applyForcesLoaded %d %d %d %d %f %f\n", x, y, x2, y2, valueX, valueY );
	
	if( x > x2 )
		std::swap( x, x2 );
	if( y > y2 )
		std::swap( y, y2 );
	
	//In the rectangle
	for( int iX = x; iX <= x2; ++iX)
	{
		for( int iY = y; iY <= y2; ++iY)
		{
			GLuint iSample = (iY*(nbSamplesX)+ iX);
			//std::cout<<"iX : "<<iX<<" \tiY : "<<iY<<" \tindex : "<<index<<std::endl;
			forces[iSample*4+0]+= valueX ;
			forces[iSample*4+1]+= valueY ;
		}
	}
}

//Reset all forces 
void Simulation::resetForces()
{
	for( int iX = 0; iX < nbSamplesX; ++iX)
	{
		for( int iY = 0; iY < nbSamplesY; ++iY)
		{
			GLuint iSample = (iY*(nbSamplesX)+ iX);
			forces[iSample*4+0] = 0.0 ;
			forces[iSample*4+1] = 0.0;
		}
	}
}


// Constructor
Simulation::Simulation(Scene * scene, bool surface=true, GLfloat size=2.0, GLfloat density=1000.0, GLfloat viscosity=0.0, GLuint nbSamplesX=10, GLuint nbSamplesY=10, GLuint nbSamplesZ=1, GLuint nbParticlesCoef=1, GLuint defaultShaderID=1, GLuint spriteShaderID=1, bool solidWalls=true)
{
    // About the data stored on centers of the MAC grid
    
        // Samples number on X axis
    this->nbSamplesX=nbSamplesX;
        // Samples number on Y axis
    this->nbSamplesY=nbSamplesY;
           // Samples number on Z axis
    this->nbSamplesZ=nbSamplesZ;
    if (surface) this->nbSamplesZ=1;
        // Total samples number
    this->nbSamples=this->nbSamplesX*this->nbSamplesY*this->nbSamplesZ;
        // Samples positions
    this->samples=new GLfloat[this->nbSamples*4];
        // Samples pressures
    this->pressures=new GLfloat[this->nbSamples];
        // Samples colors
    this->colors=new GLfloat[this->nbSamples*4];
           
           
    // About the data stored on borders of the MAC grid
    
        // Number of left/right cell borders
    this->nbBordersX=(this->nbSamplesX+1)*this->nbSamplesY*this->nbSamplesZ;
        // Number of bottom/top cell borders
    this->nbBordersY=this->nbSamplesX*(this->nbSamplesY+1)*this->nbSamplesZ;
        // Number of back/front cell borders
    this->nbBordersZ=this->nbSamplesX*this->nbSamplesY*(this->nbSamplesZ+1);
    if (surface) this->nbBordersZ=0;
        // X velocity components, sampled on left/right borders
    this->velocitiesX=new GLfloat[this->nbBordersX];
        // Y velocity components, sampled on bottom/top borders
    this->velocitiesY=new GLfloat[this->nbBordersY];
        // Z velocity components, sampled on back/front borders
    this->velocitiesZ=new GLfloat[this->nbBordersZ];
    
    
    // About the particles used for visualization (advected on the whole grid)
    
        // Coefficient for density of particles per cell
    this->nbParticlesCoef=nbParticlesCoef;
        // Number of rendered particles (=nbSamples*[4,8]^nbParticlesCoef)
    this->nbParticles=this->nbSamples*(GLuint)pow(8.0, (int) this->nbParticlesCoef);
    if (surface) this->nbParticles=this->nbSamples*(GLuint)pow(4.0, (int) this->nbParticlesCoef);
        // Particle positions
    this->particles=new GLfloat[this->nbParticles*4];
        // Grid interpolated particles colors
    this->particleColors=new GLfloat[this->nbParticles*4];
        // Grid interpolated particles velocities
    this->particleVelocities=new GLfloat[this->nbParticles*4*2];      
        
        
    // Geometric parameters
    
        // True if 2D simulation
    this->surface=surface;
        // Width of the grid
    this->size=size;
        // Cell side length (square or cubes)
    this->h=size/float(nbSamplesX);
        // Start point left/bottom/back of the grid
    //this->offset[0]=-size/2.0f; this->offset[1]=-size/2.0f; this->offset[2]=-size/2.0f; OLD
    // Start point left/bottom/back of the grid
	this->offset[0]=-size/2.0;
	this->offset[1]=-size*(GLfloat(nbSamplesY)/GLfloat(nbSamplesX))/2.0;
	this->offset[2]=-size*(GLfloat(nbSamplesZ)/GLfloat(nbSamplesX))/2.0;

    if (surface) this->offset[2]=size/2.0f;
        // Offset from left/bottom/back of grid cell
    this->offsetInCell=this->h/2.0f;
        // True if grid boundary is expected to be solid
    this->solidWalls=solidWalls;
        // Scale coefficient to display vector lengths
    this->vectorScale=1.0f/3.0f;
        
    // Simulation parameters
    
        // Expected density for the fluid
    this->density=1000.0f;
        // Expected viscosity for the fluid
    this->viscosity=0.0f;
        // Vector of acceleration due to gravity
    this->g[0]= 0.0f; this->g[1]=-9.81f; this->g[2]=0.0f;
        
        
    // Additional required data
    
        // Scene in which to draw the simulation
    this->scene=scene;
        // Shader used for test rendering
    this->defaultShaderID=defaultShaderID;
        // Shader used for particles
    this->spriteShaderID=spriteShaderID;
        // Mersenne-Twister random generator
    this->rand=new MTRand(1995);


    // Objects for OpenGL visualization
    
        // For rendering samples (position/color checking)
    this->objectSamples=NULL;  
        // For rendering velocities interpolated on samples
    this->objectVelocitiesCenters=NULL;  
        // For rendering velocities components
    this->objectVelocitiesBorders=NULL; 
        // For rendering particles (position/color checking)
     this->objectParticles=NULL;   
        // For rendering particles interpolated velocities
    this->objectParticleVelocities=NULL;
		// For rendering pressures (pressure checking)
	this->objectPressures=NULL;
    
	//Sample forces
	this->forces = new GLfloat[this->nbSamples*4];

	//Types
	this->types = new GLuint[this->nbSamples];

    // Printing a few informations
    if (surface) std::cout<<"Simulation on a surface ("<<this->nbSamplesX<<"*"<<this->nbSamplesY<<") !"<<std::endl;
    else         std::cout<<"Simulation in a volume ("<<this->nbSamplesX<<"*"<<this->nbSamplesY<<"*"<<this->nbSamplesZ<<") !"<<std::endl;
    std::cout<<"nbSamples="<<this->nbSamples<<std::endl;
    std::cout<<"nbParticles="<<this->nbParticles<<std::endl;


    // Initialisation of buffers

    // Initialisation of required positions colors and velocities for simulation
    this->initSimulation();

    // Initialisation of particles for visualization
    this->initVisualization();
}


// Cleans memory for Simulation
Simulation::~Simulation()
{
    delete [] samples;
    delete [] pressures;
    delete [] colors;
	delete [] forces;
	delete [] types;
    
    delete [] velocitiesX;
    delete [] velocitiesY;
    delete [] velocitiesZ;
    
    delete [] particles;
    delete [] particleColors;
    delete [] particleVelocities;
}

// Sets integer indices of the cell position pos overlaps
// 3D not implemented
void Simulation::getCell(GLfloat * pos, int * iX, int * iY, int * iZ)
{
	(*iX)=int(floor((pos[0]-offset[0])/h));
	(*iY)=int(floor((pos[1]-offset[1])/h));
}


// Inits the data sampled on centers and borders of the MAC grid
void Simulation::initSimulation()
{
    // Inits the samples positions and colors
    this->initSamples();
    
    // Inits the cell borders velocity components
    this->initVelocitiesBorders();
    std::cout<<"INIT"<<std::endl;
    
    
    //---------------LOAD PPM--------------------
    GLuint height;
    GLuint width;
    unsigned char *	image_ppm = loadPPM("../ppm/image1.ppm", &height, &width);
    
    //-------Get informations from ppm ----------
    this->threeColorImageHandler(image_ppm);
    
    //----------Get forces from XML -------------
    this->loadForcesFrom("../xml/forces1.xml");
}


// Inits the particles spread over the grid for visualization
void Simulation::initVisualization()
{
    this->initParticles();
}

// Inits the samples positions and marks 4 zones with different colors
// 3D working
void Simulation::initSamples()
{
	for (GLuint iZ=0 ; iZ<nbSamplesZ ; iZ++)
	{
	    for (GLuint iY=0 ; iY<nbSamplesY ; iY++)
	    {
	        for (GLuint iX=0 ; iX<nbSamplesX ; iX++)
	        {
	            GLuint index=iZ*(nbSamplesY*nbSamplesX)+iY*nbSamplesX+iX;
	            samples[index*4+0]=this->offset[0]+this->offsetInCell+iX*this->h;
	            samples[index*4+1]=this->offset[1]+this->offsetInCell+iY*this->h;
	            if (surface) samples[index*4+2]=this->offset[2]+iZ*this->h;
	            else         samples[index*4+2]=this->offset[2]+this->offsetInCell+iZ*this->h;
	            samples[index*4+3]=1.0f;
	            
	            
	            colors[index*4+0]=0.0f;
	            colors[index*4+1]=0.0f;
	            colors[index*4+2]=0.0f;
	            colors[index*4+3]=1.0f;
	            
	            
	            // Right : red
	            if (samples[index*4+0]>0.0f) colors[index*4+0]=1.0f;
	            // Top : Green
	            if (samples[index*4+1]>0.0f) colors[index*4+1]=1.0f;
	            // Bottom-left : blue / Top-right : white
	            if ((int(colors[index*4+0]+colors[index*4+1])%2)==0)
	                colors[index*4+2]=1.0f;
				

				// Forces
	            forces[index*4+0]=0.0f;
				forces[index*4+1]=0.0f;//-2.0f;
				forces[index*4+2]=0.0f;
				forces[index*4+3]=0.0f;
				
				// Pour l'image1.ppm ! Forces qui poussent vers l'arbre.
				/*
				if (iX>=0 && iX<= 58)
					forces[index*4+0]=4.0f;
					
				if (iX>=64 && iX<= 120)
					forces[index*4+0]=-4.0f;
					
				if (iX>=58 && iX<= 64)
					forces[index*4+1]=4.0f;
				*/
				/*
				if (iX>=(nbSamplesX/2) && iX<= ((nbSamplesX/2) + 2) )
					forces[index*4+1]=10.0f;
					
				if (iX % ((nbSamplesX)/4) == 0 )
					forces[index*4+1]=-0.1f;
				*/
				// Pressures
				pressures[index]=1.0;

				// Types
				
				types[index] = AIR;
				/*
				//if(iY <= (nbSamplesY/2 + 5) && iY >= (nbSamplesY/2) && iX >= (nbSamplesX/2 +1) && iX <= (3*nbSamplesX/4 + 3))
				
				if(iX >= (nbSamplesX/2) && iY >= (nbSamplesY/2))
				{
					types[index]=FLUID;
				}*/
	        }
	    }
	}
}


// Inits velocity directional components (sampled on cell borders)
// 3D not implemented
void Simulation::initVelocitiesBorders()
{
	GLuint iBordersX, iBordersY;
    for (GLuint iY=0 ; iY<nbSamplesY ; iY++)
	{
	    for (GLuint iX=0 ; iX<(nbSamplesX+1) ; iX++)
	    {   
	        iBordersX=iY*(nbSamplesX+1)+iX;
	        
	        // Inits circular velocity field (velocityX=-borderPosition[1])
            //velocitiesX[iBordersX]=-(offset[1]+offsetInCell+iY*h)/size;
            
            // Inits diagonally growing velocity field (velocityX=k*borderPosition[0])
            velocitiesX[iBordersX]=0;//((float)iX/(float)nbSamplesX)/size;
            
            // If boundaries are solid, normal velocity component must be null
            /*if ( (solidWalls) && ((iX==0) || (iX==nbSamplesX)) )
            {
           		velocitiesX[iBordersX] = 0.0;
            }*/
        }
    }
    for (GLuint iY=0 ; iY<(nbSamplesY+1) ; iY++)
	{
	    for (GLuint iX=0 ; iX<nbSamplesX ; iX++)
	    {   
	        iBordersY=iY*nbSamplesX+iX;
	        
	        // Inits circular velocity field (velocityY=borderPosition[0])
            //velocitiesY[iBordersY]=(offset[0]+offsetInCell+iX*h)/size;
            
            // Inits diagonally growing velocity field (velocityY=k*borderPosition[1])
           	velocitiesY[iBordersY]=0;//((float)iY/(float)nbSamplesY)/size;
            
            // If boundaries are solid, normal velocity component must be null
            /*if ( (solidWalls) && ((iY==0) || (iY==nbSamplesY)) )
            {
            	velocitiesY[iBordersY] = 0.0;
            }*/
        }
    }

	enforceVelocitiesBorders();
}

// Init the particles positions, and interpolates colors and velocities from grid
// 3D not implemented
void Simulation::initParticles()
{
    //GLuint nb=(GLuint) pow(4.0, (int)nbParticlesCoef);
    GLuint nbAxis=(GLuint)pow(2.0, (int)nbParticlesCoef);  
    
    GLfloat halfCell=h/2.0f;
    GLfloat step=h/float(nbAxis);
    GLfloat sideOffsetX=h/(4.0f*nbAxis);
    	
	GLuint iParticles=0;
	GLfloat * particlesTmp=new GLfloat[this->nbParticles*4]; //new
	// Inits 4 particles in 2D cells (optimal distribution)
	// if (nbParticlesCoef==1) : corners of a losange in the square cell
	// if (nbParticlesCoef==2) : twice as more ...
	for (GLuint iSamples=0 ; iSamples<nbSamples ; iSamples++)
	{
		if (types[iSamples]==FLUID)
		{
			for (GLuint iY=0 ; iY<nbAxis ; iY++)
			{
				GLfloat offsetY=-halfCell + step/2.0 + iY*step;
				for (GLuint iX=0 ; iX<nbAxis ; iX++)
				{
					GLfloat offsetX=-halfCell + sideOffsetX + iX*step;
					if (iY%2==1) offsetX+=step/2.0;
					particlesTmp[iParticles*4+0]=samples[iSamples*4+0]+offsetX;
					particlesTmp[iParticles*4+1]=samples[iSamples*4+1]+offsetY;
					particlesTmp[iParticles*4+2]=samples[iSamples*4+2];
					particlesTmp[iParticles*4+3]=samples[iSamples*4+3];
					iParticles++;
				}
			}
		}
	}
	nbParticles=iParticles;
	delete [] particles;
	delete [] particleColors;
	delete [] particleVelocities;
	this->particles=new GLfloat[this->nbParticles*4];
	this->particleColors=new GLfloat[this->nbParticles*4];
	this->particleVelocities=new GLfloat[this->nbParticles*4*2];
	
	for (GLuint iParticlesCoord=0 ; iParticlesCoord<4*this->nbParticles ; iParticlesCoord++)
		this->particles[iParticlesCoord]=particlesTmp[iParticlesCoord];
		
	delete [] particlesTmp;
	for (GLuint iParticles=0 ; iParticles<this->nbParticles ; iParticles++)
	{
	    // Moves particles randomly and at close range around initial position
	    GLfloat radius= (float)getRand(*(this->rand), true, 0.0, 2.0*sideOffsetX);
	    GLfloat angle=(float)getRand(*(this->rand), true, 0.0, 2.0*M_PI);
	    particles[iParticles*4+0]+=radius*cos(angle);
	    particles[iParticles*4+1]+=radius*sin(angle);

        // Interpolates color and velocity bilinearly for each particle
        for (GLuint iCoord=0 ; iCoord<4 ; iCoord++)
        {
            particleColors[iParticles*4+iCoord]=0.0;
            if (iCoord==3) particleColors[iParticles*4+iCoord]=1.0;
            particleVelocities[iParticles*4+iCoord]=0.0;
        }
        interpolateFromCenters(colors, &(particles[iParticles*4+0]), &(particleColors[iParticles*4+0]));
        
        /*
        // Couleurs de Nadine : blanc jaune rouge
        particleColors[iParticles*4+0] = 1.0;
        particleColors[iParticles*4+3] = 1.0;
        if(iParticles < 2 * nbParticles/3) particleColors[iParticles*4+1] = 1.0;
        if(iParticles < 1 * nbParticles/3) particleColors[iParticles*4+2] = 1.0;
        */
        interpolateFromBorders(velocitiesX, velocitiesY, velocitiesZ, &(particles[iParticles*4+0]), &(particleVelocities[iParticles*4+0]));
    }
}


// Builds an object to visualize the grid (borders of cells as lines)
// 3D bugged
void Simulation::buildGrid(Object * object)
{
    std::cout<<"Building on GPU : Grid (borders)"<<std::endl;
    GLuint tmpNbSamplesZ=nbSamplesZ;
    if (surface) tmpNbSamplesZ=0;

    object->nbVertices=(nbSamplesX+1)*(nbSamplesY+1)*(tmpNbSamplesZ+1);
    object->nbIndices=(nbBordersX+nbBordersY+nbBordersZ)*2;
    GLfloat * vertices = new GLfloat[object->nbVertices*4];
    GLuint * indices = new GLuint[object->nbIndices];
    
	GLuint iIndices=0;
	for (GLuint iZ=0 ; iZ<(tmpNbSamplesZ+1) ; iZ++)
	{
	    for (GLuint iY=0 ; iY<(nbSamplesY+1) ; iY++)
	    {
	        for (GLuint iX=0 ; iX<(nbSamplesX+1) ; iX++)
	        {
	            GLuint index=iZ*((nbSamplesY+1)*(nbSamplesX+1))
	                        +iY*                (nbSamplesX+1)
	                        +iX;
	            // Builds a vertex at each cell corner
	            vertices[index*4+0]=offset[0]+iX*h;
	            vertices[index*4+1]=offset[1]+iY*h;
	            vertices[index*4+2]=offset[2]+iZ*h;
	            vertices[index*4+3]=1.0;
	            
	            // Builds x oriented lines
	            if (iX<nbSamplesX) // x oriented segments
	            {
	                indices[iIndices++]= iZ   *((nbSamplesY+1)*(nbSamplesX+1))
	                                   + iY   *                (nbSamplesX+1)
	                                   + iX;
	                indices[iIndices++]= iZ   *((nbSamplesY+1)*(nbSamplesX+1))
	                                   + iY   *                (nbSamplesX+1)
	                                   +(iX+1);
	            }
	            // Builds y oriented lines
	            if (iY<nbSamplesY) // y oriented segments
	            {
	                indices[iIndices++]= iZ   *((nbSamplesY+1)*(nbSamplesX+1))
	                                   + iY   *                (nbSamplesX+1)
	                                   + iX;
	                indices[iIndices++]= iZ   *((nbSamplesY+1)*(nbSamplesX+1))
	                                   +(iY+1)*                (nbSamplesX+1)
	                                   + iX;
	            }
	            // Builds z oriented lines
	            if ((!surface) && (iZ<nbSamplesZ)) // z oriented segments
	            {
	                indices[iIndices++]= iZ   *((nbSamplesY+1)*(nbSamplesX+1))
	                                   + iY   *                (nbSamplesX+1)
	                                   + iX;
	                indices[iIndices++]=(iZ+1)*((nbSamplesY+1)*(nbSamplesX+1))
	                                   + iY*                   (nbSamplesX+1)
	                                   + iX;
	            }
	        }
	    }
	}
	// Sends the data into buffers on the GPU
    object->sendPrimitives(vertices, indices);
}


// Builds an object to visualize the samples (positions and colors)
// The object will then be accessed and updated each frame (advectColors(...))
// 3D working
void Simulation::buildSamples(Object * object)
{
    std::cout<<"Building on GPU : samples positions"<<std::endl;
    
    object->nbVertices=nbSamples;
    
    // No indices are necessary since we draw GL_POINTS
	GLuint * indices=NULL;
	
    // Sends the data into buffers on the GPU
    object->sendPrimitives(samples, indices);
    object->sendColors(colors);
    
    // The colors buffer will be updated each frame so is stored in dynamic memory
    glBindBuffer(GL_ARRAY_BUFFER, objectSamples->colorsVboId);
    glBufferData(GL_ARRAY_BUFFER, objectSamples->nbVertices*4*sizeof(GLfloat), colors, GL_DYNAMIC_DRAW);
}




// Builds an object to visualize the velocity components as colors (at borders)
// The object will then be accessed and updated each frame (advectVelocities(...))
// 3D not implemented
void Simulation::buildVelocitiesBorders(Object * object)
{
    std::cout<<"Building on GPU : cell borders velocities"<<std::endl;
    
    object->nbVertices=(nbBordersX+nbBordersY+nbBordersZ);
	GLuint * indices=NULL;
	
    GLfloat * vertices = new GLfloat[object->nbVertices*4];
    GLfloat * colors = new GLfloat[object->nbVertices*4];

	GLuint iBorders=0;
	GLuint iBordersX, iBordersY;
    for (GLuint iY=0 ; iY<nbSamplesY ; iY++)
	{
	    for (GLuint iX=0 ; iX<(nbSamplesX+1) ; iX++)
	    {   
	        iBordersX=iY*(nbSamplesX+1)+iX;
            // Velocity x components are sampled on at the middle 
            // of the yz oriented cell borders (y in 2D)
            vertices[iBorders*4+0]=offset[0]+iX*h;
            vertices[iBorders*4+1]=offset[1]+offsetInCell+iY*h;
            vertices[iBorders*4+2]=offset[2];
            vertices[iBorders*4+3]=1.0;
            
            // Velocity x components are black to red
            GLfloat normVelocityX=(velocitiesX[iBordersX]);
            if (normVelocityX<0.0) normVelocityX=-normVelocityX;
            colors[iBorders*4+0]=normVelocityX;
            colors[iBorders*4+1]=0.0;
            colors[iBorders*4+2]=0.0;
            colors[iBorders*4+3]=1.0;
            
            iBorders++;
        }
    }
    
    for (GLuint iY=0 ; iY<(nbSamplesY+1) ; iY++)
	{
	    for (GLuint iX=0 ; iX<nbSamplesX ; iX++)
	    {   
            iBordersY=iY*nbSamplesX+iX;
            // Velocity y components are sampled on at the middle 
            // of the xz oriented cell borders (x in 2D)
            vertices[iBorders*4+0]=offset[0]+offsetInCell+iX*h;
            vertices[iBorders*4+1]=offset[1]+iY*h;
            vertices[iBorders*4+2]=offset[2];
            vertices[iBorders*4+3]=1.0;
            
            // Velocity y components are black to green
            GLfloat normVelocityY=(velocitiesY[iBordersY]);
            if (normVelocityY<0.0) normVelocityY=-normVelocityY;
            colors[iBorders*4+0]=0.0;
            colors[iBorders*4+1]=normVelocityY;
            colors[iBorders*4+2]=0.0;
            colors[iBorders*4+3]=1.0;

            iBorders++;
        }
    }

    // Sends the data into buffers on the GPU
    object->sendPrimitives(vertices, indices);
    object->sendColors(colors);
    
    // The colors buffer will be updated each frame so is stored in dynamic memory
    glBindBuffer(GL_ARRAY_BUFFER, object->colorsVboId);
    glBufferData(GL_ARRAY_BUFFER, object->nbVertices*4*sizeof(GLfloat), colors, GL_DYNAMIC_DRAW);
}


// Builds an object to visualize the velocity vectors interpolated at cell centers
// The object will then be accessed and updated each frame (advectVelocities(...))
// 3D not implemented
void Simulation::buildVelocitiesCenters(Object * object)
{
    std::cout<<"Building on GPU : cell centers velocities"<<std::endl;
    
    object->nbVertices=nbSamples*2;
    object->nbIndices=object->nbVertices;
    
    GLfloat * vertices = new GLfloat[object->nbVertices*4];
    GLfloat * colors = new GLfloat[object->nbVertices*4];
    GLuint * indices = new GLuint[object->nbIndices];
    
    for (GLuint iY=0 ; iY<nbSamplesY ; iY++)
    {
        for (GLuint iX=0 ; iX<nbSamplesX ; iX++)
        {
            GLuint iSamples=iY*nbSamplesX+iX;
            GLuint index=iSamples*2+0;
            // Position of vectors beginning
            vertices[index*4+0]=samples[iSamples*4+0];
            vertices[index*4+1]=samples[iSamples*4+1];
            vertices[index*4+2]=samples[iSamples*4+2];
            vertices[index*4+3]=1.0f;
            colors[index*4+0]=0.0f;
            colors[index*4+1]=0.0f;
            colors[index*4+2]=1.0f;
            colors[index*4+3]=1.0f;
            indices[index]=index;
                    
            index=iSamples*2+1;
            GLuint indexVelocitiesXLeft  =iY     *(nbSamplesX+1)+ iX;
            GLuint indexVelocitiesXRight =iY     *(nbSamplesX+1)+(iX+1);
            GLuint indexVelocitiesYBottom=iY     * nbSamplesX   + iX;
            GLuint indexVelocitiesYTop   =(iY+1) * nbSamplesX   + iX;

            GLfloat velocity[]={0.0f, 0.0f, 0.0f, 0.0f};
            // X component is interpolated from left and right
            velocity[0]=(velocitiesX[indexVelocitiesXLeft]
                        +velocitiesX[indexVelocitiesXRight])/2.0f;
            // Y component is interpolated from bottom and top           
            velocity[1]=(velocitiesY[indexVelocitiesYBottom]
                        +velocitiesY[indexVelocitiesYTop])/2.0f;
                        
            // The vector is printed inversed and transparent at the end (for effect)
            vertices[index*4+0]=samples[iSamples*4+0]-velocity[0]*vectorScale;
            vertices[index*4+1]=samples[iSamples*4+1]-velocity[1]*vectorScale;
            vertices[index*4+2]=samples[iSamples*4+2]-velocity[2]*vectorScale;
            vertices[index*4+3]=1.0;
            colors[index*4+0]=0.0;
            colors[index*4+1]=0.0;
            colors[index*4+2]=0.0;
            colors[index*4+3]=1.0;
            indices[index]=index;
        }  
    }
    
    // Sends the data into buffers on the GPU
    object->sendPrimitives(vertices, indices);
    object->sendColors(colors);
    
    // The vertices buffer will be updated each frame so is stored in dynamic memory
    glBindBuffer(GL_ARRAY_BUFFER, object->vboId);
    glBufferData(GL_ARRAY_BUFFER, object->nbVertices*4*sizeof(GLfloat), vertices, GL_DYNAMIC_DRAW);
}


// Builds an object to visualize the particles 
// The object will then be accessed and updated each frame (advectParticles(...))
// 3D working theorically
void Simulation::buildParticles(Object * object)
{
    std::cout<<"Building on GPU : particles"<<std::endl;
    
    object->nbVertices=nbParticles;
	GLuint * indices=NULL;
	
    // Sends the data into buffers on the GPU
    object->sendPrimitives(particles, indices);
    object->sendColors(particleColors);
    
    // The vertices buffer will be updated each frame so is stored in dynamic memory
    glBindBuffer(GL_ARRAY_BUFFER, object->vboId);
    glBufferData(GL_ARRAY_BUFFER, object->nbVertices*4*sizeof(GLfloat), particles, GL_DYNAMIC_DRAW);
    // The colors buffer will be updated each frame so is stored in dynamic memory
    glBindBuffer(GL_ARRAY_BUFFER, object->colorsVboId);
    glBufferData(GL_ARRAY_BUFFER, object->nbVertices*4*sizeof(GLfloat), particleColors, GL_DYNAMIC_DRAW);    
}


// Builds an object to visualize the particles interpolated velocities
// The object will then be accessed and updated each frame (advectParticles(...))
// 3D not implemented
void Simulation::buildParticleVelocities(Object * object)
{
    std::cout<<"Building on GPU : particle velocities"<<std::endl;
    
    object->nbVertices=nbParticles*2;
    object->nbIndices=object->nbVertices;
    
	GLfloat * colors = new GLfloat[object->nbVertices*2*4];
	GLuint * indices = new GLuint[object->nbIndices];
	
	for (GLuint iParticles=0 ; iParticles<nbParticles ; iParticles++)
	{
	    indices[iParticles*2+0]=iParticles*2+0;
	    indices[iParticles*2+1]=iParticles*2+1;
	    for (GLuint iCoord=0 ; iCoord<4 ; iCoord++)
	    {
	        // Vector start on particle location
	        particleVelocities[(iParticles*2+0)*4+iCoord]=particles[iParticles*4+iCoord];
	        // Vector starts with the particle interpolated color
	        colors[(iParticles*2+0)*4+iCoord]=particleColors[iParticles*4+iCoord];
	        
	        // Bilinear inteprolation of velocity
	        GLfloat particleVelocity[]={0.0, 0.0, 0.0, 0.0};
	        if (surface)
            {
                interpolateFromBorders(velocitiesX, velocitiesY, velocitiesZ, &(particles[iParticles*4+0]), particleVelocity);
            }else{}
            
            // The vector is printed inversed and transparent at the end (for effect)
	        particleVelocities[(iParticles*2+1)*4+iCoord]=particles[iParticles*4+iCoord]-(particleVelocity[iCoord]*vectorScale);
	        colors[(iParticles*2+1)*4+iCoord]=0.0;
	    }
	}
	
    // Sends the data into buffers on the GPU
    object->sendPrimitives(particleVelocities, indices);
    object->sendColors(colors);
    
    // The vertices buffer will be updated each frame so is stored in dynamic memory
    glBindBuffer(GL_ARRAY_BUFFER, object->vboId);
    glBufferData(GL_ARRAY_BUFFER, object->nbVertices*4*sizeof(GLfloat), particleVelocities, GL_DYNAMIC_DRAW);
    
    // The colors buffer will be updated each frame so is stored in dynamic memory
    glBindBuffer(GL_ARRAY_BUFFER, object->colorsVboId);
    glBufferData(GL_ARRAY_BUFFER, object->nbVertices*4*sizeof(GLfloat), colors, GL_DYNAMIC_DRAW);    
}


// Creates, builds and add to draw list an object for samples visualization
void Simulation::drawSamples()
{
    objectSamples=new Object(GL_POINTS);
    GLuint storedObjectSamples=scene->storeObject(objectSamples);
    this->buildSamples(objectSamples);
    GLuint samplesID=scene->addObjectToDraw(storedObjectSamples);
    scene->setDrawnObjectShaderID(samplesID, defaultShaderID);
}


// Creates, builds and add to draw list an object for grid visualization
void Simulation::drawGrid()  
{
    Object * objectGrid=new Object(GL_LINES);
    GLuint storedObjectGrid=scene->storeObject(objectGrid);
    this->buildGrid(objectGrid);
    GLuint gridID=scene->addObjectToDraw(storedObjectGrid);  
    scene->setDrawnObjectShaderID(gridID, defaultShaderID);
}
   
   
// Creates, builds and add to draw list an object for velocity components visualization
void Simulation::drawVelocitiesBorders() 
{   
    this->objectVelocitiesBorders=new Object(GL_POINTS);
    GLuint storedObjectVelocitiesBorders=scene->storeObject(objectVelocitiesBorders);
	this->buildVelocitiesBorders(objectVelocitiesBorders);
	GLuint velocitiesBordersID=scene->addObjectToDraw(storedObjectVelocitiesBorders);
	scene->setDrawnObjectShaderID(velocitiesBordersID, defaultShaderID);
}


// Creates, builds and add to draw list an object for velocity vectors visualization
void Simulation::drawVelocitiesCenters() 
{ 
    this->objectVelocitiesCenters=new Object(GL_LINES);
    GLuint storedObjectVelocitiesCenters=scene->storeObject(objectVelocitiesCenters);
    this->buildVelocitiesCenters(objectVelocitiesCenters);
    GLuint velocitiesCentersID=scene->addObjectToDraw(storedObjectVelocitiesCenters);
    scene->setDrawnObjectShaderID(velocitiesCentersID, defaultShaderID);
}

   
// Creates, builds and add to draw list an object for particles visualization
void Simulation::drawParticles()
{
    this->objectParticles=new Object(GL_POINTS);
    GLuint storedObjectParticles=scene->storeObject(objectParticles);
    this->buildParticles(objectParticles);
    GLuint particlesID=scene->addObjectToDraw(storedObjectParticles);
    
    // Particles are drawn with point sprites 
    // (textures aligned on screen for each point)
    // Black parts of sprites are transparent (blending enabled, depth test disabled)
    scene->setDrawnObjectShaderID(particlesID, spriteShaderID);
    GLuint starTextureID=loadTexture("../textures/star.ppm");
    scene->setDrawnObjectTextureID(particlesID, 0, starTextureID);
}

// Creates, builds and add to draw list an object for particles velocities visualization
void Simulation::drawParticlesVelocities()
{
    this->objectParticleVelocities=new Object(GL_LINES);
    GLuint storedObjectParticleVelocities=scene->storeObject(objectParticleVelocities);
    this->buildParticleVelocities(objectParticleVelocities);
    GLuint particleVelocitiesID=scene->addObjectToDraw(storedObjectParticleVelocities);
    scene->setDrawnObjectShaderID(particleVelocitiesID, defaultShaderID);
}



// Interpolates (bi/tri)linearly a data field sampled on cell centers
// 3D not implemented
void Simulation::interpolateFromCenters(GLfloat * data, GLfloat * position, GLfloat * result)
{
    GLfloat xCorner=(position[0]-offset[0])/h;
	GLfloat yCorner=(position[1]-offset[1])/h;

    if ((xCorner<0.0f) || (xCorner>nbSamplesX) || (yCorner<0.0f) || (yCorner>nbSamplesY))
        return;
        
    // Left-bottom closest sample
    GLfloat x=xCorner-0.5f;
    GLfloat y=yCorner-0.5f;
    
    // Distances to left-bottom closest sample
    GLfloat toLeft  =x-floor(x);
    GLfloat toBottom=y-floor(y);    
    
    // Corresponding indices on x/y axis
    int iX=int(floor(x));
    int iY=int(floor(y));
    
    // Four neighbours indices
    int x0y0= iY   *nbSamplesX+iX  ;
    int x1y0= iY   *nbSamplesX+iX+1;
    int x0y1=(iY+1)*nbSamplesX+iX  ;
    int x1y1=(iY+1)*nbSamplesX+iX+1;
    if (x<0.0)              { x0y0=x1y0; x0y1=x1y1;}
    if (x>(nbSamplesX-1.0)) { x1y0=x0y0; x1y1=x0y1;}
    if (y<0.0)              { x0y0=x0y1; x1y0=x1y1;}
    if (y>(nbSamplesY-1.0)) { x0y1=x0y0; x1y1=x1y0;}
    // Bilinear interpolation (4 components) from the four neighbours
    biLinearInterpolation(4, &(data[x0y0*4]), &(data[x1y0*4]), &(data[x0y1*4]), &(data[x1y1*4]), toLeft, toBottom, result);
}


// Interpolates (bi/tri)linearly a data field sampled on cell borders per components
// 3D not implemented
void Simulation::interpolateFromBorders(GLfloat * dataX, GLfloat * dataY, GLfloat * dataZ, GLfloat * position, GLfloat * result)
{
    GLfloat xCorner=(position[0]-offset[0])/h;
	GLfloat yCorner=(position[1]-offset[1])/h;

    if ((xCorner<0.0f) || (xCorner>nbSamplesX) || (yCorner<0.0f) || (yCorner>nbSamplesY))
        return;

    GLfloat x=xCorner;
    GLfloat y=yCorner-0.5f;
    // Distances to left closest sample and corresponding y  
    GLfloat toLeft  =x-floor(x);
    GLfloat toBottom=y-floor(y);
    
    // Corresponding indices on x/y axis
    int indexLeft  =int(floor(x));
    int indexBottom=int(floor(y)); 
       
    // Corresponding index in dataX
    int indexDataXLeft=indexBottom*(nbSamplesX+1)+indexLeft;
    
    // Four dataX neighbours indices
    int x0y0=indexDataXLeft;
    int x1y0=indexDataXLeft+1;
    int x0y1=indexDataXLeft +(nbSamplesX+1);
    int x1y1=indexDataXLeft+(nbSamplesX+1)+1;
    if (y<0.0)              { x0y0=x0y1; x1y0=x1y1;}
    if (y>(nbSamplesY-1.0)) { x0y1=x0y0; x1y1=x1y0;}
    // Bilinear interpolation (1 x components) from the four neighbours
    biLinearInterpolation(1, &(dataX[x0y0]), &(dataX[x1y0]), &(dataX[x0y1]), &(dataX[x1y1]), toLeft, toBottom, &(result[0]));
    
    x=xCorner-0.5f;
    y=yCorner;
    // Distances to bottom closest sample and corresponding x   
    toLeft  =x-floor(x);
    toBottom=y-floor(y);
    
    // Corresponding indices on x/y axis     
    indexBottom=int(floor(y));
    indexLeft  =int(floor(x));
    
    // Corresponding index in dataY
    int indexDataYBottom=indexBottom*nbSamplesX+indexLeft;
    // Four dataY neighbours indices   
    x0y0=indexDataYBottom;
    x1y0=indexDataYBottom+1;
    x0y1=indexDataYBottom+nbSamplesX;
    x1y1=indexDataYBottom+nbSamplesX+1;
    if (x<0.0)              { x0y0=x1y0; x0y1=x1y1;}
    if (x>(nbSamplesX-1.0)) { x1y0=x0y0; x1y1=x0y1;}
    // Bilinear interpolation (1 y components) from the four neighbours
    biLinearInterpolation(1, &(dataY[x0y0]), &(dataY[x1y0]), &(dataY[x0y1]), &(dataY[x1y1]), toLeft, toBottom, &(result[1]));
}


// Simple and intuitive (but inaccurate) time integration method (Forward Euler order 1)
// if direction is -1 : evaluates value is value at t-dt
void Simulation::forwardEuler1stOdrTimeIntegration(GLfloat direction, GLfloat * knownValue, GLfloat * knownVariation, GLfloat dt, GLfloat * resultValue)
{
	for (GLuint iCoord=0 ; iCoord<4 ; iCoord++)
    {
		resultValue[iCoord] = knownValue[iCoord] + direction * knownVariation[iCoord] * dt;
	}
}


// Sufficienty accurate time integration method (Runge-Kutta order 2)
// if direction is -1 : evaluates value is value at t-dt
void Simulation::rungeKutta2ndOdrTimeIntegration(GLfloat direction, GLfloat * knownValue, GLfloat * knownVariation, GLfloat dt, GLfloat * resultValue)
{
	GLfloat intermediateValue[4];
    // This time "variation" is approximated at what it would be at t+dt/2
	for (GLuint iCoord=0 ; iCoord<4 ; iCoord++)
    {
		intermediateValue[iCoord] = knownValue[iCoord] + direction * knownVariation[iCoord] * dt * 0.5f;
	}
	GLfloat intermediateVariation[]={0.f, 0.f, 0.f, 0.f};
    interpolateFromBorders(velocitiesX, velocitiesY, velocitiesZ, intermediateValue, intermediateVariation);
	
	
	for (GLuint iCoord=0 ; iCoord<4 ; iCoord++)
    {
		resultValue[iCoord] = knownValue[iCoord] + direction * intermediateVariation[iCoord] * dt;
	}
	
	
}

// Sophisticated (best trade-of accuracy/cost) time integration method (Runge-Kutta order 3)
// if direction is -1 : evaluates value is value at t-dt
void Simulation::rungeKutta3rdOdrTimeIntegration(GLfloat direction, GLfloat * knownValue, GLfloat * knownVariation, GLfloat dt, GLfloat * resultValue)
{
    // Variation is approximated at several instants between t and t+dt
    // Final value uses these successive variations instead of only the starting one as in Forward Euler
    
    GLfloat intermediateValue1[4];
    for (GLuint iCoord=0 ; iCoord<4 ; iCoord++)
        intermediateValue1[iCoord]=knownValue[iCoord]+direction*(1.0f/2.0f)*dt*knownVariation[iCoord];
    GLfloat intermediateVariation1[]={0.0f, 0.0f, 0.0f, 0.0f};
    interpolateFromBorders(velocitiesX, velocitiesY, velocitiesZ, intermediateValue1, intermediateVariation1);
    
    GLfloat intermediateValue2[4];
    for (GLuint iCoord=0 ; iCoord<4 ; iCoord++)
        intermediateValue2[iCoord]=knownValue[iCoord]+direction*(3.0f/4.0f)*dt*intermediateVariation1[iCoord];        
    GLfloat intermediateVariation2[]={0.0f, 0.0f, 0.0f, 0.0f};
    interpolateFromBorders(velocitiesX, velocitiesY, velocitiesZ, intermediateValue2, intermediateVariation2);
    
    for (GLuint iCoord=0 ; iCoord<4 ; iCoord++)
        resultValue[iCoord]=knownValue[iCoord]
                           +direction*(2.0f/9.0f)*dt*knownVariation[iCoord]
                           +direction*(3.0f/9.0f)*dt*intermediateVariation1[iCoord]
                           +direction*(4.0f/9.0f)*dt*intermediateVariation2[iCoord];
}


// Interpolation scheme choice for the whole simulation
// if direction is -1 : evaluates value is value at t-dt
void Simulation::integrate(GLfloat direction, GLfloat * knownValue, GLfloat * knownVariation, GLfloat dt, GLfloat * resultValue)
{
    for (GLuint iCoord=0 ; iCoord<4 ; iCoord++)
        resultValue[iCoord]=knownValue[iCoord];
        
    //forwardEuler1stOdrTimeIntegration(direction, knownValue, knownVariation, dt, resultValue);
    //rungeKutta2ndOdrTimeIntegration(direction, knownValue, knownVariation, dt, resultValue);
    rungeKutta3rdOdrTimeIntegration(direction, knownValue, knownVariation, dt, resultValue);
}


// Moves the particles by interpolating velocity and integrating the new position from it, changes the color as well by interpolating on the color field
// 3D working theorically
void Simulation::advectParticles(GLfloat dt)
{
    //GLfloat * velocityColors = new GLfloat[nbParticles*4*2];

	//Set all cells non-solid to AIR
	for (GLuint iSamples=0 ; iSamples<nbSamples ; iSamples++)
	{
		if (types[iSamples]!=SOLID)
		{
			types[iSamples] = AIR;
		}
	}
    
    for (GLuint iParticles=0 ; iParticles<nbParticles ; iParticles++)
	{
		//Set the cell to FLUID  
        int iX, iY, iZ;
        getCell( &(particles[iParticles*4+0]), &iX, &iY, &iZ);
        if (type(iX, iY, 0)!=SOLID)
        	types[ iY*nbSamplesX+iX ] = FLUID;
	
        GLfloat knownVelocity[]={0.0, 0.0, 0.0, 0.0};
        interpolateFromBorders(velocitiesX, velocitiesY, velocitiesZ, &(particles[iParticles*4+0]), knownVelocity);

        // What is the new position after dt ?
        integrate(1.0,                            // direction
                &(particles[iParticles*4+0]),     // startValue
                  knownVelocity,                  // startVariation 
                  dt,                             // dt
                &(particles[iParticles*4+0])); 		// resultValue     
        
        // What is the velocity at this new position ?
        GLfloat newVelocity[]={0.0, 0.0, 0.0, 0.0};
        interpolateFromBorders(velocitiesX, velocitiesY, velocitiesZ, &(particles[iParticles*4+0]), newVelocity);
        for (GLuint iCoord=0 ; iCoord<4 ; iCoord++)
        {
            particleVelocities[(iParticles*2+0)*4+iCoord]=particles[iParticles*4+iCoord];
            particleVelocities[(iParticles*2+1)*4+iCoord]=particles[iParticles*4+iCoord]-(newVelocity[iCoord]*vectorScale);
        }
        
        
        // What is the color at this new position ?
        /* ADVECT COLOR DESACTIVATED
        interpolateFromCenters(colors, &(particles[iParticles*4+0]), &(particleColors[iParticles*4+0]));
        for (GLuint iCoord=0 ; iCoord<4 ; iCoord++)
        {
            velocityColors[(iParticles*2+0)*4+iCoord]=particleColors[iParticles*4+iCoord];
	        velocityColors[(iParticles*2+1)*4+iCoord]=0.0;
	    }*/
    }
    if (objectPressures != NULL)
	{
		GLfloat * colorsType = new GLfloat[nbSamples * 4];
		for(GLuint i = 0; i < nbSamples; ++i)
		{
			if(types[i] == FLUID)
			{
				colorsType[i * 4 + 0] = 1.0;
				colorsType[i * 4 + 1] = 0.0; 
				colorsType[i * 4 + 2] = 0.0;
				colorsType[i * 4 + 3] = 1.0;
			}
			if(types[i] == AIR)
			{
				colorsType[i * 4 + 0] = 0.0;
				colorsType[i * 4 + 1] = 1.0; 
				colorsType[i * 4 + 2] = 0.0;
				colorsType[i * 4 + 3] = 1.0;
			}
		}
   		glBindBuffer(GL_ARRAY_BUFFER, objectPressures->colorsVboId);
    	glBufferData(GL_ARRAY_BUFFER, objectPressures->nbVertices*4*sizeof(GLfloat), colorsType, GL_DYNAMIC_DRAW);
    }
    // Updates the corresponding data if stored on GPU
    if (objectParticles!=NULL)
    {
        glBindBuffer(GL_ARRAY_BUFFER, objectParticles->vboId);
        glBufferData(GL_ARRAY_BUFFER, objectParticles->nbVertices*4*sizeof(GLfloat), particles, GL_DYNAMIC_DRAW);
        
        
        // Color advect DESACTIVATED
        //glBindBuffer(GL_ARRAY_BUFFER, objectParticles->colorsVboId);
        //glBufferData(GL_ARRAY_BUFFER, objectParticles->nbVertices*4*sizeof(GLfloat), particleColors, GL_DYNAMIC_DRAW);
    }
    if (objectParticleVelocities!=NULL)
    {
        glBindBuffer(GL_ARRAY_BUFFER, objectParticleVelocities->vboId);
        glBufferData(GL_ARRAY_BUFFER, objectParticleVelocities->nbVertices*4*sizeof(GLfloat), particleVelocities, GL_DYNAMIC_DRAW);
        
        // Color advect DESACTIVATED
        //glBindBuffer(GL_ARRAY_BUFFER, objectParticleVelocities->colorsVboId);
        //glBufferData(GL_ARRAY_BUFFER, objectParticleVelocities->nbVertices*4*sizeof(GLfloat), velocityColors, GL_DYNAMIC_DRAW);
    }
}


// Modifies the color samples so that color derivative is null (particles moving through the grid keep their color)
// 3D not implemented
// Semi-Lagrangian advection : 1 - interpolates velocity at center, 
//                             2 - finds the previous position for material at center
//                             3 - copies the value found at this position to center
void Simulation::advectColors(GLfloat dt)
{	
	for (GLuint iY = 0; iY < nbSamplesY; ++iY)
	{
		for (GLuint iX = 0; iX < nbSamplesX; ++iX)
		{
			GLuint iSample = iY * nbSamplesX + iX;
			GLuint indexSampleXLeft = iY * (nbSamplesX +1) + iX;
			GLuint indexSampleXRight = iY * (nbSamplesX +1) + iX +1;
			GLuint indexSampleYBottom = iY * nbSamplesX + iX;
			GLuint indexSampleYTop = (iY+1) * nbSamplesX + iX;

			GLfloat knownVelocity[]={0.0, 0.0, 0.0, 0.0};

			knownVelocity[0] = (velocitiesX[indexSampleXLeft] + velocitiesX[indexSampleXRight])/2.0f;
			knownVelocity[1] = (velocitiesY[indexSampleYBottom] + velocitiesY[indexSampleYTop])/2.0f;

			GLfloat oldPosition[]={0.0f, 0.0f, 0.0f, 0.0f};
			// What was the position before dt ?
			integrate(-1.0,                            	// direction
					&(samples[iSample*4+0]),     		// startValue
					  knownVelocity,                  	// startVariation 
					  dt,                             	// dt
           	 		  oldPosition);						// resultValue
           	
           	GLuint typeSample = type(int(oldPosition[0]), int(oldPosition[1]), 0);
           	if( typeSample == FLUID )
           	{
				GLfloat color[]={0.0, 0.0, 0.0, 1.0};
				interpolateFromCenters(colors, oldPosition, color);
				colors[iSample*4+0] = color[0];
				colors[iSample*4+1] = color[1];
				colors[iSample*4+2] = color[2];
				colors[iSample*4+3] = color[3];
			}
		}
	}
	// Updates the corresponding data if stored on GPU
	if (objectSamples != NULL)
	{
   		glBindBuffer(GL_ARRAY_BUFFER, objectSamples->colorsVboId);
    	glBufferData(GL_ARRAY_BUFFER, objectSamples->nbVertices*4*sizeof(GLfloat), colors, GL_DYNAMIC_DRAW);
    }
    
    // Réduire la taille des cellules ralentit l'effet "mélange des couleurs"
}


/// Capitalization : update velocities from centers from the centeredVel.
void Simulation::updateVelocitiesFromCenters(GLfloat * centeredVel)
{
	// Resets velocities components to 0
	for (GLuint iBordersX=0 ; iBordersX<nbBordersX ; iBordersX++)
		velocitiesX[iBordersX]=0.0;
	for (GLuint iBordersY=0 ; iBordersY<nbBordersY ; iBordersY++)
		velocitiesY[iBordersY]=0.0;
		
	// Uses the new centers velocity to update the separated velocity components
	for (GLuint iY=0 ; iY<nbSamplesY ; iY++)
	{
		for (GLuint iX=0 ; iX<nbSamplesX ; iX++)
		{
			GLuint iSamples=iY*nbSamplesX+iX;
			GLuint indexVelocitiesXLeft =iY *(nbSamplesX+1)+ iX;
			GLuint indexVelocitiesXRight =iY *(nbSamplesX+1)+(iX+1);
			GLuint indexVelocitiesYBottom=iY * nbSamplesX + iX;
			GLuint indexVelocitiesYTop =(iY+1) * nbSamplesX + iX;
			GLfloat coefLeft=0.5;
			GLfloat coefRight=0.5;
			GLfloat coefBottom=0.5;
			GLfloat coefTop=0.5;
			if (iX==0) coefLeft =1.0;
			if (iX==(nbSamplesX-1)) coefRight =1.0;
			if (iY==0) coefBottom=1.0;
			if (iY==(nbSamplesY-1)) coefTop =1.0;

			velocitiesX[indexVelocitiesXLeft] +=coefLeft *centeredVel[iSamples*4+0];
			velocitiesX[indexVelocitiesXRight] +=coefRight *centeredVel[iSamples*4+0];
			velocitiesY[indexVelocitiesYBottom]+=coefBottom*centeredVel[iSamples*4+1];
			velocitiesY[indexVelocitiesYTop] +=coefTop *centeredVel[iSamples*4+1];
		}
	}
	
	enforceVelocitiesBorders();
	
	// Updates the corresponding data if stored on GPU
	if (objectVelocitiesCenters!=NULL)
	{
		GLfloat vertices[nbSamples*4*2];
		for (GLuint iSamples=0 ; iSamples<nbSamples ; iSamples++)
		{
			for (GLuint iCoord=0 ; iCoord<4 ; iCoord++)
			{
			vertices[(iSamples*2+0)*4+iCoord]=samples[iSamples*4+iCoord];
			vertices[(iSamples*2+1)*4+iCoord]=samples[iSamples*4+iCoord]-centeredVel[iSamples*4+iCoord]*vectorScale;
			}
		}
		glBindBuffer(GL_ARRAY_BUFFER, objectVelocitiesCenters->vboId);
		glBufferData(GL_ARRAY_BUFFER, objectVelocitiesCenters->nbVertices*4*sizeof(GLfloat), vertices, GL_DYNAMIC_DRAW);
	}
	if (objectVelocitiesBorders!=NULL)
	{
		GLfloat colors[objectVelocitiesBorders->nbVertices*4];
		GLuint iBorders=0;
		for (GLuint iY=0 ; iY<nbSamplesY ; iY++)
		{
			for (GLuint iX=0 ; iX<(nbSamplesX+1) ; iX++)
			{
			GLfloat normVelocityX=velocitiesX[iY*(nbSamplesX+1)+iX];
			// if (normVelocityX<0.0) normVelocityX=-normVelocityX;
			colors[iBorders*4+0]=normVelocityX;
			colors[iBorders*4+1]=0.0;
			colors[iBorders*4+2]=-normVelocityX;
			colors[iBorders*4+3]=1.0;
			iBorders++;
			}
		}
		for (GLuint iY=0 ; iY<(nbSamplesY+1) ; iY++)
		{
			for (GLuint iX=0 ; iX<nbSamplesX ; iX++)
			{
			GLfloat normVelocityY=velocitiesY[iY*nbSamplesX+iX];
			//if (normVelocityY<0.0) normVelocityY=-normVelocityY;
			colors[iBorders*4+0]=0.0;
			colors[iBorders*4+1]=normVelocityY;
			colors[iBorders*4+2]=-normVelocityY;
			colors[iBorders*4+3]=1.0;
			iBorders++;
			}
		}
		glBindBuffer(GL_ARRAY_BUFFER, objectVelocitiesBorders->colorsVboId);
		glBufferData(GL_ARRAY_BUFFER, objectVelocitiesBorders->nbVertices*4*sizeof(GLfloat), colors, GL_DYNAMIC_DRAW);
	}
}


void Simulation::updateVelocitiesFromBorders()
{
	if (objectVelocitiesCenters!=NULL)
	{
		GLfloat vertices[nbSamples*4*2];
		for (GLuint iY=0 ; iY<nbSamplesY ; iY++)
		{
			for (GLuint iX=0 ; iX<nbSamplesX ; iX++)
			{
				GLuint iSamples=iY*nbSamplesX+iX;
				GLfloat sampleVelocity[]={0.0, 0.0, 0.0, 0.0};
				// Interpolates the velocity at cell center (sampleVelocity)
				GLuint indexVelocitiesXLeft =iY *(nbSamplesX+1)+ iX;
				GLuint indexVelocitiesXRight =iY *(nbSamplesX+1)+(iX+1);
				GLuint indexVelocitiesYBottom=iY * nbSamplesX + iX;
				GLuint indexVelocitiesYTop =(iY+1) * nbSamplesX + iX;
				sampleVelocity[0]=(velocitiesX[indexVelocitiesXLeft] + velocitiesX[indexVelocitiesXRight])/2.0;
				sampleVelocity[1]=(velocitiesY[indexVelocitiesYBottom] + velocitiesY[indexVelocitiesYTop])/2.0;
				for (GLuint iCoord=0 ; iCoord<4 ; iCoord++)
				{
					vertices[(iSamples*2+0)*4+iCoord]=samples[iSamples*4+iCoord];
					vertices[(iSamples*2+1)*4+iCoord]=samples[iSamples*4+iCoord]-sampleVelocity[iCoord]*vectorScale;
				}
			}
		}
		glBindBuffer(GL_ARRAY_BUFFER, objectVelocitiesCenters->vboId);
		glBufferData(GL_ARRAY_BUFFER, objectVelocitiesCenters->nbVertices*4*sizeof(GLfloat), vertices, GL_DYNAMIC_DRAW);
	}
	if (objectVelocitiesBorders!=NULL)
	{
		GLfloat colors[objectVelocitiesBorders->nbVertices*4];
		GLuint iBorders=0;
		for (GLuint iY=0 ; iY<nbSamplesY ; iY++)
		{
			for (GLuint iX=0 ; iX<(nbSamplesX+1) ; iX++)
			{
			GLfloat normVelocityX=velocitiesX[iY*(nbSamplesX+1)+iX];
			if (normVelocityX<0.0) normVelocityX=-normVelocityX;
			colors[iBorders*4+0]=normVelocityX*10.0;//new 10
			colors[iBorders*4+1]=0.0;
			colors[iBorders*4+2]=0.0;
			colors[iBorders*4+3]=1.0;
			iBorders++;
			}
		}
		for (GLuint iY=0 ; iY<(nbSamplesY+1) ; iY++)
		{
			for (GLuint iX=0 ; iX<nbSamplesX ; iX++)
			{
				GLfloat normVelocityY=velocitiesY[iY*nbSamplesX+iX];
				if (normVelocityY<0.0) normVelocityY=-normVelocityY;
				colors[iBorders*4+0]=0.0;
				colors[iBorders*4+1]=normVelocityY*100.0;//new 10
				colors[iBorders*4+2]=0.0;
				colors[iBorders*4+3]=1.0;
				iBorders++;
				std::cout<<"vélocity : " << normVelocityY <<std::endl;
			}
		}
		glBindBuffer(GL_ARRAY_BUFFER, objectVelocitiesBorders->colorsVboId);
		glBufferData(GL_ARRAY_BUFFER, objectVelocitiesBorders->nbVertices*4*sizeof(GLfloat), colors, GL_DYNAMIC_DRAW);
	}
}



// Modifies the velocity samples so that velocity derivative is null (particles moving through the grid keep their velocity)
// 3D not implemented
// Semi-Lagrangian advection : 1 - interpolates velocity at center, 
//                             2 - finds the previous position for material at center
//                             3 - copies the value found at this position to center
void Simulation::advectVelocities(GLfloat dt)
{
	GLfloat * newVelocities = new GLfloat[ nbSamples * 4 ];
	GLfloat sampleVelocity[]={0.0f, 0.0f, 0.0f, 0.0f};
	GLfloat oldPosition[]={0.0f, 0.0f, 0.0f, 1.0f};

	for (GLuint iY=0 ; iY<nbSamplesY ; ++iY)
	{
		for (GLuint iX=0 ; iX<nbSamplesX ; ++iX)
		{
			GLuint iSamples = iY * nbSamplesX + iX;
			//Interpolates the velocity at cell center (sampleVelocity)
			GLuint indexVelocitiesXLeft = iY * (nbSamplesX +1) + iX;
			GLuint indexVelocitiesXRight = iY * (nbSamplesX +1) + iX +1;
			GLuint indexVelocitiesYBottom = iY * nbSamplesX + iX;
			GLuint indexVelocitiesYTop = (iY+1) * nbSamplesX + iX;

			sampleVelocity[0]=(velocitiesX[indexVelocitiesXLeft]+velocitiesX[indexVelocitiesXRight])/2.0;
			sampleVelocity[1]=(velocitiesY[indexVelocitiesYBottom]+velocitiesY[indexVelocitiesYTop])/2.0;

			// at t-dt where was the material now arriving at sample location ?
			integrate(-1.0, &(samples[iSamples*4]), sampleVelocity, dt, oldPosition);

			// What is the velocity of the material now arriving on sample location ?
			newVelocities[iSamples*4+0]=sampleVelocity[0];
			newVelocities[iSamples*4+1]=sampleVelocity[1];
			newVelocities[iSamples*4+2]=sampleVelocity[2];
			newVelocities[iSamples*4+3]=sampleVelocity[3];
			int iXvirtual, iYvirtual, iZvirtual;
			getCell(oldPosition, &iXvirtual, &iYvirtual, &iZvirtual);
			if (type(iXvirtual, iYvirtual, 0)==FLUID)
				interpolateFromBorders(velocitiesX, velocitiesY, velocitiesZ, oldPosition, &(newVelocities[iSamples*4+0]));
		}
	}

	updateVelocitiesFromCenters(newVelocities);
}

void Simulation::applyForces(GLfloat dt)
{
	GLfloat * newVelocities = new GLfloat[nbSamples * 4];
	for( GLuint iY = 0; iY < nbSamplesY ; ++iY)
	{
		for(GLuint iX=0 ; iX < nbSamplesX ; ++iX)
		{
			
			GLuint iSamples = iY * nbSamplesX + iX;
			GLfloat samplesVelocity[] = {0.0, 0.0, 0.0, 0.0};
		
			GLuint indexVelocityXLeft = iY * (nbSamplesX + 1) + iX;
			GLuint indexVelocityXRight = iY * (nbSamplesX + 1) + iX + 1;
			GLuint indexVelocityYBottom = iY * nbSamplesX + iX;
			GLuint indexVelocityYTop = (iY +1 ) * nbSamplesX + iX;
		
			samplesVelocity[0] = (velocitiesX[indexVelocityXLeft] + velocitiesX[indexVelocityXRight] )/2.0f;
			samplesVelocity[1] = (velocitiesY[indexVelocityYTop] + velocitiesY[indexVelocityYBottom] )/2.0f;
		
			for( GLuint iCoord = 0; iCoord < 4; iCoord++)
			{
				newVelocities[iSamples *4 + iCoord] = samplesVelocity[iCoord];
			}
			
			GLuint typeSample = type(int(iX), int(iY), 0);
           	if( typeSample == FLUID )
			{
				for( GLuint iCoord = 0; iCoord < 4; iCoord++)
           		{
           			newVelocities[iSamples *4 + iCoord] += dt * forces[iSamples *4 + iCoord];
           		}
			}
			
		 }
	}
	
	updateVelocitiesFromCenters(newVelocities);

	//enforceVelocitiesBorders();
}

void Simulation::buildPressures(Object * object)
{
	std::cout<<"Building on GPU : samples positions"<<std::endl;
	object->nbVertices=nbSamples;
	// No indices are necessary since we draw GL_POINTS
	GLuint * indices=NULL;
	GLfloat * colors = new GLfloat[object->nbVertices * 4];
	for (GLuint iSamples=0 ; iSamples<nbSamples ; iSamples++)
	{
		colors[iSamples * 4 + 0] = pressures[iSamples];
		colors[iSamples * 4 + 0] = pressures[iSamples];
		colors[iSamples * 4 + 0] = pressures[iSamples];
		colors[iSamples * 4 + 0] = 0.0;
	}
	// Sends the data into buffers on the GPU
	object->sendPrimitives(samples, indices);
	object->sendColors(colors);
	// The colors buffer will be updated each frame so is stored in dynamic memory
	glBindBuffer(GL_ARRAY_BUFFER, objectPressures->colorsVboId);
	glBufferData(GL_ARRAY_BUFFER, objectPressures->nbVertices * 4 * sizeof(GLfloat), colors, GL_DYNAMIC_DRAW);
}

// Creates, builds and add to draw list an object for pressures visualization
void Simulation::drawPressures()
{
	this->objectPressures=new Object(GL_POINTS);
	GLuint storedObjectPressures=scene->storeObject(objectPressures);
	this->buildPressures(objectPressures);
	GLuint pressuresID=scene->addObjectToDraw(storedObjectPressures);
	scene->setDrawnObjectShaderID(pressuresID, defaultShaderID);
}

void Simulation::setDivergences(double * divergences)
{
	GLfloat scale = 1.f/h;
	for (GLuint iY=0 ; iY< nbSamplesY ; iY++)
	{
		for (GLuint iX=0 ; iX< nbSamplesX ; iX++)
		{
			GLuint iSamples = iY * nbSamplesX + iX;

			GLuint indexVelocityXLeft = iY * (nbSamplesX + 1) + iX;
			GLuint indexVelocityXRight = iY * (nbSamplesX + 1) + iX + 1;
			GLuint indexVelocityYBottom = iY * nbSamplesX + iX;
			GLuint indexVelocityYTop = (iY +1 ) * nbSamplesX + iX;

			divergences[iSamples] = scale*(velocitiesX[indexVelocityXRight] - velocitiesX[indexVelocityXLeft] + velocitiesY[indexVelocityYTop] - velocitiesY[indexVelocityYBottom]);
		}
	}
}

// Permet d'éviter les regroupements de particules grâce à pressures
void Simulation::project(GLfloat dt)
{
	double scaleRhs = 1.0 / h;
	double rhs[nbSamples]; //Right Hand Side
	
	double scaleA = dt / (density * h *h );
	double Aright[nbSamples];
	double Adiag[nbSamples];
	double Atop[nbSamples];
	
	//Initialisation
	for(GLuint iSamples = 0; iSamples < nbSamples; ++iSamples)
	{
		rhs[iSamples] = 0.0;
		Adiag[iSamples] = 0.0;
		Atop[iSamples] = 0.0;
		Aright[iSamples] = 0.0;
	}
	
	double velocitySolids = 0.0;
	setDivergences(rhs);
	
	for (GLuint iY = 0; iY < nbSamplesY; ++iY)
	{
		for (GLuint iX = 0; iX < nbSamplesX; ++iX)
		{
			GLuint iSamples=iY*nbSamplesX+iX;
			
			GLuint indexVelocityXLeft = iY * (nbSamplesX + 1) + iX;
			GLuint indexVelocityXRight = iY * (nbSamplesX + 1) + iX + 1;
			GLuint indexVelocityYBottom = iY * nbSamplesX + iX;
			GLuint indexVelocityYTop = (iY +1 ) * nbSamplesX + iX;
			
			GLuint typeSample = type(int(iX), int(iY), 0);
			GLuint typeNeighbourLeft = type(int(iX)-1, int(iY), 0);
			GLuint typeNeighbourRight = type(int(iX)+1, int(iY), 0);
			GLuint typeNeighbourBottom = type(int(iX), int(iY)-1, 0);
			GLuint typeNeighbourTop = type(int(iX), int(iY)+1, 0);
			
			//Sets Right Hand Side						
			if( typeSample == FLUID ) 
			{
				rhs[iSamples] = -rhs[iSamples];
				
				if( typeNeighbourLeft == SOLID ) // Solid on the left
					rhs[iSamples] -= scaleRhs * (velocitiesX[indexVelocityXLeft] - velocitySolids);
		
				if( typeNeighbourRight == SOLID ) // Solid on the right
					rhs[iSamples] += scaleRhs * (velocitiesX[indexVelocityXRight] - velocitySolids);
			
				if( typeNeighbourBottom == SOLID ) // Solid on the right
					rhs[iSamples] -= scaleRhs * (velocitiesY[indexVelocityYBottom] - velocitySolids);
			
				if( typeNeighbourTop == SOLID ) // Solid on the right
					rhs[iSamples] += scaleRhs * (velocitiesY[indexVelocityYTop] - velocitySolids);
				
				if( typeNeighbourLeft != SOLID ) Adiag[iSamples] += scaleA;		
				if( typeNeighbourRight != SOLID ) Adiag[iSamples] += scaleA;		
				if( typeNeighbourBottom != SOLID ) Adiag[iSamples] += scaleA;
				if( typeNeighbourTop != SOLID ) Adiag[iSamples] += scaleA;
			
				if( typeNeighbourTop == FLUID ) Atop[iSamples] = -scaleA;
				if( typeNeighbourRight == FLUID ) Aright[iSamples] = -scaleA;
			}
			else
			{
				rhs[iSamples] = 0.f;
			}
		}	
	}
	conjugateGradient(Adiag, Aright, Atop, rhs);
	
	GLfloat scale = dt / (density * h );
	for (GLuint iY = 0; iY < nbSamplesY; ++iY)
	{
		for (GLuint iX = 0; iX < nbSamplesX; ++iX)
		{
			GLuint iSamples=iY*nbSamplesX+iX;
			
			GLuint indexVelocityXLeft = iY * (nbSamplesX + 1) + iX;
			GLuint indexVelocityXRight = iY * (nbSamplesX + 1) + iX + 1;
			GLuint indexVelocityYBottom = iY * nbSamplesX + iX;
			GLuint indexVelocityYTop = (iY +1 ) * nbSamplesX + iX;

			if( types[iSamples] == FLUID ) 
			{
				velocitiesX[indexVelocityXLeft] -= scale * pressures[iSamples]; 
				velocitiesX[indexVelocityXRight] += scale * pressures[iSamples]; 
				velocitiesY[indexVelocityYBottom] -= scale * pressures[iSamples]; 
				velocitiesY[indexVelocityYTop] += scale * pressures[iSamples];
			}

			/*if( solidWalls )
			{
				if(iX == 0) velocitiesX[indexVelocityXLeft] = 0.0;
				if(iX == nbSamplesX - 1) velocitiesX[indexVelocityXRight] = 0.0;
				if(iY == 0) velocitiesY[indexVelocityYBottom] = 0.0;
				if(iY == nbSamplesY - 1) velocitiesY[indexVelocityYTop] = 0.0;
			}*/
		}
	}

	enforceVelocitiesBorders();
	updateVelocitiesFromBorders();
	
	/*if(objectVelocitiesBorders != NULL)
	{
		GLfloat * colors = new GLfloat[objectVelocitiesBorders->nbVertices * 4];
		GLuint iBorders = 0;
		for (GLuint iY = 0; iY < nbSamplesY; ++iY)
		{
			for (GLuint iX = 0; iX < nbSamplesX + 1; ++iX)
			{
				GLfloat normVelocityX = velocitiesX[iY * (nbSamplesX + 1) + iX];
				if(normVelocityX < 0.0)
				{
					normVelocityX = -normVelocityX;
				}
				colors[(iBorders*4 + 0)] = normVelocityX;
				colors[(iBorders*4 + 1)] = 0.0;
				colors[(iBorders*4 + 2)] = 0.0;
				colors[(iBorders*4 + 3)] = 1.0;
				iBorders++;
			}
		}
		
		for (GLuint iY = 0; iY < nbSamplesY +1 ; ++iY)
		{
			for (GLuint iX = 0; iX < nbSamplesX; ++iX)
			{
				GLfloat normVelocityY = velocitiesY[iY * (nbSamplesX) + iX];
				if(normVelocityY < 0.0)
				{
					normVelocityY = -normVelocityY;
				}
				colors[(iBorders*4 + 0)] = 0.0;
				colors[(iBorders*4 + 1)] = normVelocityY;
				colors[(iBorders*4 + 2)] = 0.0;
				colors[(iBorders*4 + 3)] = 1.0;
				iBorders++;
			}
		}
		glBindBuffer(GL_ARRAY_BUFFER, objectVelocitiesBorders->colorsVboId);
    	glBufferData(GL_ARRAY_BUFFER, objectVelocitiesBorders->nbVertices*4*sizeof(GLfloat), colors, GL_DYNAMIC_DRAW);
	}*/
}

void Simulation::multiplySparseMatrix(double * Adiag, double * Aright, double * Atop, Eigen::VectorXd v, Eigen::VectorXd * result)
{
	for (GLuint iY = 0; iY < nbSamplesY ; ++iY)
	{
		for (GLuint iX = 0; iX < nbSamplesX; ++iX)
		{
			GLuint iSamples=iY*nbSamplesX+iX;
			
			(*result)[iSamples] = Adiag[iSamples] * v[iSamples];
			//Right
			if(iX < (nbSamplesX - 1) )
				(*result)[iSamples] += Aright[iSamples] * v[iSamples+1];
			//Top
			if(iY < (nbSamplesY - 1) )
				(*result)[iSamples] += Atop[iSamples] * v[iSamples + nbSamplesX] ;
			// Left
			if(iX > 0 )
				(*result)[iSamples] += Aright[iSamples-1] * v[iSamples-1];
			//Bottom
			if(iY > 0 )
				(*result)[iSamples] += Atop[iSamples-nbSamplesX] * v[iSamples - nbSamplesX];
		}
	}
}

// Fills pressures from sparse matrix and rhs=b
void Simulation::conjugateGradient(double * Adiag, double * Aright, double * Atop, double * b)
{
    GLuint maxIterations=100;
    double tol=1e-6;
    // Initializes pressures to null
    Eigen::VectorXd p=Eigen::VectorXd::Zero(nbSamples);
    // Copy r from b and return if r null
    Eigen::VectorXd r(nbSamples);
    for (GLuint iSamples=0 ; iSamples<nbSamples ; iSamples++)
         r[iSamples]=b[iSamples];

	// rInfinityNorm is the maximum of absolute values of r
    double rInfinityNorm=r.lpNorm<Eigen::Infinity>();
    if (rInfinityNorm<=tol)
    {
         for (GLuint iSamples=0 ; iSamples<nbSamples ; iSamples++)
              pressures[iSamples]=(GLfloat)p[iSamples];
         std::cout<<"Stability reached."<<std::endl;
         return;
    }
    
    // Builds Preconditionner
	Eigen::VectorXd precon=Eigen::VectorXd::Zero(nbSamples);
	MICPreconditioner(Adiag, Aright, Atop, &precon);
	// Apply preconditionner from r to z
	Eigen::VectorXd z=Eigen::VectorXd::Zero(nbSamples);
	applyPreconditioner(Aright, Atop, precon, r, &z);
	// Copy z to s
	Eigen::VectorXd s(z);
	double sigma=z.dot(r);
	GLuint cpt = 0;;
	for (GLuint i=0 ; i<maxIterations ; i++)
	{
		// Multiply matrix A to s to get z;
		multiplySparseMatrix(Adiag, Aright, Atop, s, &z);
		double alpha=0;
		double div=z.dot(s);
		if (div!=0.0) alpha=sigma/div;
		p+=alpha*s;
		r-=alpha*z;
		rInfinityNorm=r.lpNorm<Eigen::Infinity>();
		if (rInfinityNorm<=tol)
		{
		     i=maxIterations;
		}
		else
		{
		     applyPreconditioner(Aright, Atop, precon, r, &z);
		     double beta=0;
		     double newSigma=z.dot(r);
		     if (sigma!=0.0) beta=newSigma/sigma;
		     s=z+beta*s;
		     sigma=newSigma;
		}
		++cpt;
	}
	std::cout<<"Nombre iterations : "<<cpt<<std::endl;
	for (GLuint iSamples=0 ; iSamples<nbSamples ; iSamples++)
		pressures[iSamples]=(GLfloat)p[iSamples];
}

void Simulation::MICPreconditioner(double * Adiag, double * Aright, double * Atop, Eigen::VectorXd* precon)
{
    double tau=0.97;
    double safetyConstant=0.25;
    for (int iY=0 ; iY<(int)nbSamplesY ; iY++)
    {
        for (int iX=0 ; iX<(int)nbSamplesX ; iX++)
        {
            int iS=iY*nbSamplesX+iX;
            int iSLeft=iY*nbSamplesX+(iX-1);
            int iSBottom=(iY-1)*nbSamplesX+iX;
            double e=Adiag[iS];
            if (iX>0)
            {
                 e-=pow(Aright[iSLeft]*(*precon)[iSLeft], 2)
                   +tau*(Aright[iSLeft] *Atop[iSLeft] *pow((*precon)[iSLeft] , 2));
            }
            if (iY>0)
            {
                 e-=pow(Atop[iSBottom]*(*precon)[iSBottom], 2)
                   +tau*(Atop [iSBottom]*Aright[iSBottom]*pow((*precon)[iSBottom], 2));
            }
            if (e<(safetyConstant*Adiag[iS])) e=Adiag[iS];
            if (e!=0.0) (*precon)[iS]=1.0/sqrt(e);
        }
    }
}
// Multiplies preconditionneur to r and stores it in z
void Simulation::applyPreconditioner(double * Aright, double * Atop, Eigen::VectorXd precon, Eigen::VectorXd r, Eigen::VectorXd * z)
{
    double t=0.0;
	Eigen::VectorXd q=Eigen::VectorXd::Zero(nbSamples);
	for (int iY=0 ; iY<(int)nbSamplesY ; iY++)
	{
	  for (int iX=0 ; iX<(int)nbSamplesX ; iX++)
	  {
		  int iS=iY*nbSamplesX+iX;
		  int iSLeft=iY*nbSamplesX+(iX-1);
		  int iSBottom=(iY-1)*nbSamplesX+iX;
		  t=r[iS];
		  if (iX>0) t-=Aright[iSLeft]*precon[iSLeft] *q[iSLeft];
		  if (iY>0) t-=Atop[iSBottom]*precon[iSBottom]*q[iSBottom];
		  q[iS]=t*precon[iS];
	  }
	}
	for (int iY=((int)nbSamplesY-1) ; iY>=0; iY--)
	{
	  for (int iX=((int)nbSamplesX-1) ; iX>=0 ; iX--)
	  {
		  int iS=iY*nbSamplesX+iX;
		  int iSRight=iY*nbSamplesX+(iX+1);
		  int iSTop=(iY+1)*nbSamplesX+iX;
		  t=q[iS];
		  if (iX<((int)nbSamplesX-1)) t-=Aright[iS]*precon[iS]*(*z)[iSRight];
		  if (iY<((int)nbSamplesY-1)) t-=Atop[iS] *precon[iS]*(*z)[iSTop];
		  (*z)[iS]=t*precon[iS];
	  }
	}
}


void Simulation::buildForces(Object * object)
{
	std::cout<<"Building on GPU : cell forces"<<std::endl;
	object->nbVertices=nbSamples*2;
	object->nbIndices=object->nbVertices;
	GLfloat * vertices = new GLfloat[object->nbVertices*4];
	GLfloat * colors = new GLfloat[object->nbVertices*4];
	GLuint * indices = new GLuint[object->nbIndices];
	for (GLuint iY=0 ; iY<nbSamplesY ; iY++)
	{
		for (GLuint iX=0 ; iX<nbSamplesX ; iX++)
		{
		    GLuint iSamples=iY*nbSamplesX+iX;
		    GLuint index=iSamples*2+0;
		    // Position of vectors beginning
		    vertices[index*4+0]=samples[iSamples*4+0];
		    vertices[index*4+1]=samples[iSamples*4+1];
		    vertices[index*4+2]=samples[iSamples*4+2];
		    vertices[index*4+3]=1.0;
		    colors[index*4+0]=1.0;
		    colors[index*4+1]=1.0;
		    colors[index*4+2]=1.0;
		    colors[index*4+3]=1.0;
		    indices[index]=index;
		    index=iSamples*2+1;
		    // The vector is printed from white to red
		    vertices[index*4+0]=samples[iSamples*4+0]+forces[iSamples*4+0]*vectorScale;
		    vertices[index*4+1]=samples[iSamples*4+1]+forces[iSamples*4+1]*vectorScale;
		    vertices[index*4+2]=samples[iSamples*4+2]+forces[iSamples*4+2]*vectorScale;
		    vertices[index*4+3]=1.0;
		    colors[index*4+0]=1.0;
		    colors[index*4+1]=0.0;
		    colors[index*4+2]=0.0;
		    colors[index*4+3]=1.0;
		    indices[index]=index;
		}
	}
	// Sends the data into buffers on the GPU
	object->sendPrimitives(vertices, indices);
	object->sendColors(colors);
}

void Simulation::drawForces()
{
	Object * objectForces = new Object(GL_LINES);
	GLuint storedObjectForces=scene->storeObject(objectForces);
	this->buildForces(objectForces);
	GLuint forcesID=scene->addObjectToDraw(storedObjectForces);
	scene->setDrawnObjectShaderID(forcesID, defaultShaderID);

}

// Builds an object to visualize the types
// 3D working
void Simulation::buildTypes(Object* object)
{
	std::cout<<"Building on GPU : samples types"<<std::endl;
	object->nbVertices=nbSamples;
	// No indices are necessary since we draw GL POINTS
	GLuint * indices = NULL;
	GLfloat* colors = new GLfloat[object->nbVertices * 4];
	for (GLuint iSamples=0 ; iSamples<nbSamples ; iSamples++)
	{
		colors[iSamples * 4 + 0] = 0.0;
		colors[iSamples * 4 + 1] = 0.0;
		colors[iSamples * 4 + 2] = 0.0;
		colors[iSamples * 4 + 3] = 1.0;
		colors[iSamples * 4 + types[iSamples]] = 1.0;
	}
	// Sends the data into buffers on the GPU
	object->sendPrimitives(samples, indices);
	object->sendColors(colors);
}

// Creates, builds and add to draw list an object for types visualization
void Simulation::drawTypes()
{
	Object* objectTypes = new Object(GL_POINTS);
	GLuint storedObjectTypes=scene->storeObject(objectTypes);
	this->buildTypes(objectTypes);
	GLuint typesID=scene->addObjectToDraw(storedObjectTypes);
	scene->setDrawnObjectShaderID(typesID, defaultShaderID);
}

// Returns type flag for neighbours with provided coordinates (works for boundaries as well)
// 3D not implemented
GLuint Simulation::type(int iX, int iY, int iZ)
{
	if ( (iX<0) || (iX>=(int)nbSamplesX) || (iY<0) || (iY>=(int)nbSamplesY) )
	{
		if (solidWalls) 
			return 2;
		else
			return 1;
	}
	return types[iY * nbSamplesX + iX];
}

// Cancels normal velocities on solid boundaries
// 3D not implemented
void Simulation::enforceVelocitiesBorders()
{
	for (GLuint iY=0 ; iY<nbSamplesY ; iY++)
	{
		for (GLuint iX=0 ; iX<nbSamplesX ; iX++)
		{
			GLuint iSamples = iY * nbSamplesX + iX;
			GLuint indexVelocitiesXLeft =iY * (nbSamplesX+1)+ iX;
			GLuint indexVelocitiesXRight =iY * (nbSamplesX+1)+(iX+1);
			GLuint indexVelocitiesYBottom=iY * nbSamplesX + iX;
			GLuint indexVelocitiesYTop = (iY+1) * nbSamplesX + iX;
			// If solid cell
			if (types[iSamples]==SOLID)
			{
				// normal velocity component must be null
				velocitiesX[indexVelocitiesXLeft] =0.0;
				velocitiesX[indexVelocitiesXRight] =0.0;
				velocitiesY[indexVelocitiesYBottom]=0.0;
				velocitiesY[indexVelocitiesYTop] =0.0;
			}
			if (solidWalls)
			{
				if (iX==0) velocitiesX[indexVelocitiesXLeft] =0.0;
				if (iX==(nbSamplesX-1)) velocitiesX[indexVelocitiesXRight] =0.0;
				if (iY==0) velocitiesY[indexVelocitiesYBottom]=0.0;
				if (iY==(nbSamplesY-1)) velocitiesY[indexVelocitiesYTop] =0.0;
			}
		}
	}
}

// Simulation update
void Simulation::update()
{
	float dt = 1.0f/30.0f;
    // Colors semi-lagrangian advection (and potential visualization update)
    this->advectColors(dt);
    
    // Velocities semi-lagrangian advection (and potential visualization update)
    this->advectVelocities(dt);
    
	// Extern forces are applied
	this->applyForces(dt);
	// Velocities are uptated to satisfy incompressiblity and limit conditons
	this->project(dt);


}


// Particle visualization for the simulation
void Simulation::render()
{
    // Moves particles according to current velocity field 
    // Updates potential visualization
    this->advectParticles(1.0f/30.0f);
}

