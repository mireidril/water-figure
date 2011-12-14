// Application.cpp
// Template Fluids Simulation
// N. Dommanget dommange@univ-mlv.fr

#include "Application.hpp"
#include "Scene.hpp"
#include "Tools.hpp"
#include "Camera.hpp"
#include "Simulation.hpp"

#include <iostream>
#define _USE_MATH_DEFINES
#include <math.h>
#include <sstream>


// Default constructor
Application::Application()
{
	bStoreFrame = false;
    this->init();
}

// Cleans before the application can be closed
Application::~Application()
{
	//SDL_RemoveTimer(this->renderTimer);
	SDL_RemoveTimer(this->moveTimer);
    SDL_RemoveTimer(this->animateTimer);
    SDL_Quit();
    if (scene!=NULL) delete scene;
    if (simulation!=NULL) delete simulation;
}


// Sets the application parameters and does all the initialisation
void Application::init()
{
    this->scene=NULL;
    this->simulation = NULL;
    
    // False as long as we don't want to quit the application
    this->done=false;
	
	this->moveFlags[0]=0.0;
	this->moveFlags[1]=0.0;
	this->moveFlags[2]=0.0;
	
	// Times measurements initialisation
	this->cntFrame=0;
	this->lastStartTime=0;
	this->frameDuration=0; 
	
    // Move counter
    this->cntMove=0;
    eventCnt=0;
    // Mouse is seeable at first
    this->mouse=true;
    
    // Window size initialization (for windowed mode)
	this->windowedWidth=800;
	this->windowedHeight=600;
	
	// Mouse position, pressed position and scroll data initilaization
    // Positions : floats, origin in center, sides at +/-1 or more
	this->xMousePosition=0.0;
	this->yMousePosition=0.0;
	this->xMouseLeftDownPosition=0.0; 
	this->yMouseLeftDownPosition=0.0;

	this->scroll=0;
	
	// Initialisation of SDL and creation of OpenGL context
    initSDLOpenGL();
    
    // Customize a few OpenGL and SDL states (after context creation)
    customizeStates();
}


// Inits SDL and OpenGL context, sets a few states
void Application::initSDLOpenGL()
{
    // Initialization of SDL
    
        // Initialize timer, audio, video, CD_ROM, and joystick.
    int sdlError=SDL_Init(SDL_INIT_EVERYTHING);
    if (sdlError<0) 
		std::cout<<"Unable to init SDL : "<<SDL_GetError()<<std::endl;

    	// Sets openGL parameters before opening the draw context
    SDL_GL_SetAttribute(SDL_GL_DOUBLEBUFFER, 1);    // Double buffering
    SDL_GL_SetAttribute(SDL_GL_DEPTH_SIZE, 16);		// Depth buffer size of 16-bit
    
    SDL_GL_SetAttribute(SDL_GL_RED_SIZE, 8);    // Color components size of 8-bit each
    SDL_GL_SetAttribute(SDL_GL_GREEN_SIZE, 8);
    SDL_GL_SetAttribute(SDL_GL_BLUE_SIZE, 8);
    SDL_GL_SetAttribute(SDL_GL_ALPHA_SIZE, 8);

	// Operating systems uniform keyboard handling 
	// (overhead, disable for performance : param->0)
	SDL_EnableUNICODE(1);	

    // Capture of screen size (for fullscreen mode)
    
    const SDL_VideoInfo* videoInfo=SDL_GetVideoInfo();
	this->fullScreenWidth=videoInfo->current_w;
	this->fullScreenHeight=videoInfo->current_h;

	// Creation of the openGL draw context
	
		// Window size
	this->width=this->windowedWidth;
	this->height=this->windowedHeight;	
	    // Options about the OpenGL window for SDL_SetVideoMode(...)  
	this->videoModeFlags = SDL_OPENGL | SDL_RESIZABLE; // Resizable OpenGL display
        // Specifies the size and other options about the window and OpenGL context
     this->drawContext = SDL_SetVideoMode(this->width, this->height, 0, this->videoModeFlags);
}


// Customize a few OPenGL states to fit the application's needs
void Application::customizeStates()
{
   // Puts the window top-left corner at the following coordinates 
    // Resizing doesn't impact so it is not very usefull
	//putenv("SDL_VIDEO_WINDOW_POS=100,0"); 
	 
    	// Glew initialisation : to register all available extentions
    GLenum glewError=glewInit();
    if (glewError!=GLEW_OK)
		std::cout<<"GLEW Error : "<<glewGetErrorString(glewError)<<std::endl;
	
	// Mouse motion will not generate events
	// Instead we will check the position ourselves when we need them
	SDL_EventState(SDL_MOUSEMOTION, SDL_IGNORE);
		

	// Initialization of the mouse position in the middle of the window	

	    // WarpMouse changes the mouse position and 
        // generates motion events which we need to ignore.
	SDL_EventState(SDL_MOUSEMOTION, SDL_IGNORE);
    SDL_WarpMouse(this->width/2, this->height/2);
        // After, we can reactivate the mouse motion events
        // But we instead choose to check directly the position 
        // ourselves when we need it (in the camera update)
        // It is better then to disable the unused events
    //SDL_EventState(SDL_MOUSEMOTION, SDL_ENABLE);


    	// Depth test
    //glEnable(GL_DEPTH_TEST);
    glDisable(GL_DEPTH_TEST);
    
        // Decides the background color used after this call
    GLfloat black[]={0.0, 0.0, 0.0, 1.0};
    setBackgroundColor(black);
    
        // Sets the with of the lines
    glLineWidth(2); 
    glPointSize(5);
    
    glEnable(GL_POINT_SMOOTH);
    glEnable(GL_POINT_SPRITE);
    glEnable(GL_PROGRAM_POINT_SIZE);

    // Turn on additive blending
    glEnable(GL_BLEND);
    glBlendFunc(GL_ONE, GL_ONE);
    //glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

        // Disables culling
    //glDisable(GL_CULL_FACE);
}


// Sets one scene for the application to draw
void Application::setScene(Scene * scene)
{
	this->scene=scene;
	resize(this->width, this->height);
}


// Initializes timers
void Application::initTimers()
{
    //Each timer event has an arbitrairy int ID
    
    //this->renderEventID=1; // Identifies our render loop event
    this->moveEventID=2; // Identifies our render loop event
    this->animateEventID=3; // Identifies our move loop event

    // Creates a timer to send a render event every 20 ms
    //this->renderTimer=SDL_AddTimer(20, genericTimer, (void*)(&(this->renderEventID)));
    
    // Creates a timer to send a move event every 20 ms
    this->moveTimer=SDL_AddTimer(20, genericTimer, (void*)(&(this->moveEventID)));    

    // Creates a timer to send a move event every 60 ms    
    this->animateTimer=SDL_AddTimer(60, genericTimer, (void*)(&(this->animateEventID)));
}


// Create an user event, give it the id (int) passed in "param",
// and registers it : it should now be send every "interval" set of time
Uint32 genericTimer(Uint32 interval, void* param)
{
    SDL_Event event;
    event.type = SDL_USEREVENT;
    event.user.code = *(int *)param; 
    event.user.data1 = 0;
    event.user.data2 = 0;
    SDL_PushEvent(&event);
    return interval;
}


// Turns ON/OFF wireframe mode
void Application::switchWireframe()
{
    std::cout<<"Wireframe switch"<<std::endl;
	GLint wireframe[2];
	glGetIntegerv(GL_POLYGON_MODE, wireframe);
	if (wireframe[0]==GL_FILL)
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	else
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
}


// Adapts the drawing to the new size of the window
// resize doesn't work on Mac os (Windows ?)
void Application::resize(GLuint w, GLuint h)
{
    //std::cout<<"Window resize  : ["<<w<<","<<h<<"]"<<std::endl;
    
    this->width=w;
    this->height=h;

	//SDL_VideoMode update (restart the OpenGL context on windows, does not work on mac os...)
	//#ifdef __APPLE__ | WIN32 // was not tested
        //this->drawContext = SDL_SetVideoMode(this->width, this->height, 0, this->videoModeFlags);
        //std::cout<<SDL_GetError()<<std::endl;
        //this->customizeStates();
    //#else
        this->drawContext = SDL_SetVideoMode(this->width, this->height, 0, this->videoModeFlags);
    //#endif
	
    
    // Viewport transformation update to fit initial window size
    glViewport(0, 0, this->width, this->height);

    // Projection transformation update
        // Keeping the angle constant
    GLfloat fovy=(float)M_PI/3.0;
    GLfloat aspectRatio=this->width/GLfloat(this->height);
    this->scene->camera->setPerspectiveFromAngle(fovy, aspectRatio);
}

// Switch to or from fullScreen mode
// Does not work on Mac os : full screen apps must be started as so.
void Application::switchFullScreen()
{
    if (this->videoModeFlags==(SDL_OPENGL | SDL_RESIZABLE))
    {
        this->videoModeFlags = SDL_OPENGL | SDL_FULLSCREEN;
        this->width=this->fullScreenWidth;
        this->height=this->fullScreenHeight;
    }
    else
    {
        this->videoModeFlags = SDL_OPENGL | SDL_RESIZABLE;
        this->width=this->windowedWidth;
        this->height=this->windowedHeight;
    }
    resize(this->width,  this->height);
}

// Switch from fps mouse mode (not limited nor seeable) to the regular mode
// Does not work on Mac os.
void Application::switchMouse()
{
    if (this->mouse)
    {
        // Mouse is not going to leave the window
	    SDL_WM_GrabInput(SDL_GRAB_ON);
        // Mouse won't be seeable
	    SDL_ShowCursor(SDL_DISABLE);
	    this->mouse=false;
    }
    else
    {
	    SDL_WM_GrabInput(SDL_GRAB_OFF);
	    SDL_ShowCursor(SDL_ENABLE);
	    this->mouse=true;
    }
}


// Sets the background color
void Application::setBackgroundColor(GLfloat * color)
{
    glClearColor(color[0], color[1], color[2], color[3]);
}


// Prints the frame duration and framerate (both averaged over 10 frames)
void Application::printFPS()
{
    float milliseconds=((float)(this->frameDuration)/1000.0f);
    float FPS=1000.0f/milliseconds; // Unit : [s^-1]= 1/seconds
    std::cout<<"frame : "<<milliseconds<<" ms       FPS : "<<FPS<<std::endl;
    
    // To print fps in the title bar
    //std::ostringstream title; title.precision(3); title.width(5);
    //title<<"frame : "<<milliseconds<<" ms       FPS : "<<FPS;
    //SDL_WM_SetCaption(title.str().data(), 0);
}


// Distributes task for the "user" kind of events 
// For example : rendrFrame action occurs if a timer event is passed
void Application::handleUserEvent(SDL_Event& event)
{
    int eventID=event.user.code;
    
    // In case a rendering event occured
    /*if (eventID==renderEventID) 
    {
        renderFrame();
    }
    else
    {*/
        // In case a moving event occured
        if (eventID==moveEventID)
        {
            //moveFPS();
		    moveSpherical();
        }
        else
        {
            // In case an animate event occured      
            if (eventID==animateEventID)
                animate();
        }
    //}
}


// Distributes task for the "key" kind of events 
// For example : std::cout when b key is pressed
// down is true when the key is pressed, false when released
void Application::handleKeyEvent(SDL_keysym& keysym, bool down)
{
    // unit is at 1 when the key is pressed and at -1 when released
    float unit=-1.0;
    if (down) unit=1.0;
    
    if (down)
    {
        switch(keysym.sym)
        {
        
            std::cout<<"down"<<std::endl;
        	case SDLK_ESCAPE:
          		this->done=true;
          	break;
          	case SDLK_SPACE:
          	  /*if(this->simulation != NULL)
				{
					this->simulation->update();
					this->simulation->render();
				}*/
			break;
          	
        	case SDLK_c :
          		std::cout<<"Key \"c\" was pressed."<<std::endl;
                this->scene->camera->switchCameraProjection();
          	break;
          	
        	case SDLK_w :
                std::cout<<"Key \"w\" was pressed."<<std::endl;
          		switchWireframe();
          	break;
          	
          	case SDLK_f :
          	    std::cout<<"Key \"f\" was pressed."<<std::endl;
          	    printFPS();
          	break;
          	
         	case SDLK_F5 :
                std::cout<<"Key \"F5\" was pressed."<<std::endl;
          		switchFullScreen();
          	break; 
          	    	
          	case SDLK_m :
                std::cout<<"Key \"m\" was pressed."<<std::endl;
          		switchMouse();
          	break;
          	
          	case SDLK_s :
                std::cout<<"Key \"s\" was pressed."<<std::endl;
          		bStoreFrame = !bStoreFrame;
          	break;
          	
          	case SDLK_o :
          		std::cout<<"Key \"o\" was pressed."<<std::endl;
          		//the real story is on the first 5 ppm
          		//the 6 and the 7 are using for the rendering of our video
		      	if(simulation->ppmPhase < 7){
		      		++ simulation->ppmPhase;
					GLuint width;
					GLuint height;
					
					//Get files string (complicated for nothing ...)
					std::string image_path;
					std::string xml_string_path;
					const char * xml_path;
					std::stringstream out;
					std::stringstream out2;
					out<<"../ppm/image"<<simulation->ppmPhase<<".ppm";
					image_path = out.str();
					out2<<"../xml/forces"<<simulation->ppmPhase<<".xml";
					
					xml_string_path = out2.str();
					xml_path = xml_string_path.data();
					//Load PPM
					unsigned char * image_ppm = loadPPM(image_path, &width, &height);
					std::cout<<"Charging "<<image_path<<std::endl;
					this->simulation->threeColorImageHandler(image_ppm);
					
					//Load forces from XML
					this->simulation->loadForcesFrom(xml_path);
					
					//Testing
					/*if(simulation->ppmPhase == 4)
					{
						simulation->drawForces();
						simulation->drawTypes();
					}*/
				}
			break;
          	
          	default:
      	    break;
      	}
    }

    switch(keysym.sym)
    {
        // Camera controls
      	// Keys are associated with x, y, z direction for move.
      	// moveFlags holds a 3D vector of direction
      	// each coordinates of moveFlags is in the end : -1, 0 or 1.   	
        
        case SDLK_UP :
        case SDLK_z :
            //std::cout<<"Forward."<<std::endl;
      		this->moveFlags[2]-=unit;
      	break;
      	
    	case SDLK_DOWN :
    	case SDLK_s :
            //std::cout<<"Backward."<<std::endl;
      		this->moveFlags[2]+=unit;
      	break;
      	
      	case SDLK_LEFT :
      	case SDLK_q :
            //std::cout<<"Left."<<std::endl;
      		this->moveFlags[0]-=unit;
      	break;
      	
      	case SDLK_RIGHT :
      	case SDLK_d :
            //std::cout<<"Right."<<std::endl;
      		this->moveFlags[0]+=unit;
      	break;
      	
      	case SDLK_PAGEUP :
            //std::cout<<"Up."<<std::endl;
      		this->moveFlags[1]+=unit;
      	break;
      	
      	case SDLK_PAGEDOWN :
            //std::cout<<"Down."<<std::endl;
      		this->moveFlags[1]-=unit;
      	break;  

    	default:
      	break;
  	}
}


// Handles mouse motion
// This function is only called in case of mouse motion events
// If mouse position is only used by the camera update, we chose
// to disable the events and to get the position directly in the 
// camera update instead, the function then useless.
void Application::handleMouseMotion(SDL_MouseMotionEvent& motion)
{
    // Stores coordinates of the mouse between -1 and 1 on x and y axis (origin : center)
    // If mouse is graped it goes beyond -1.0 and 1.0.
    // SDL window origin is top-left so negating the y axis is necessary...
    
    // When the mouse is not grabed, absolute positions can be used
	//this->xMousePosition=2.0*(motion.x/(GLfloat)this->windowWidth)-1.0;
	//this->yMousePosition=-2.0*(motion.y/(GLfloat)this->windowHeight)+1.0;
	
    // When the mouse is grabed, relatives position are the only option
	this->xMousePosition+=2.0f*motion.xrel/(GLfloat)this->width;
	this->yMousePosition+=-2.0f*motion.yrel/(GLfloat)this->height;
	
	if (motion.state&SDL_BUTTON(1))
    {	
	    this->xMouseLeftDownPosition+=2.0f*motion.xrel/(GLfloat)this->width;
	    this->yMouseLeftDownPosition+=-2.0f*motion.yrel/(GLfloat)this->height;
	}
}


// Listens to events during the whole time of the application
// and distributes corresponding tasks
void Application::handleEvent(SDL_Event& event)
{
    switch(event.type) 
    {
        // User events
        case SDL_USEREVENT:
            handleUserEvent(event);
            break;
                
        // Key presses
        case SDL_KEYDOWN:
            handleKeyEvent(event.key.keysym, true);
            break;
                
        // Key releases    
        case SDL_KEYUP:
            handleKeyEvent(event.key.keysym, false);
            break;                
                
        // Mouse motion
        // Will only ever occur if mouse motion events enabled               
        case  SDL_MOUSEMOTION: 
            handleMouseMotion(event.motion);
            break;
        
        // Scrollwheel
        case  SDL_MOUSEBUTTONDOWN: 
            if (event.button.button==SDL_BUTTON_WHEELUP)
                this->scroll++;
            if (event.button.button==SDL_BUTTON_WHEELDOWN)
                this->scroll--;            
            break;            
                
        // Window resize                  
        case  SDL_VIDEORESIZE: 
            this->windowedWidth=event.resize.w;
            this->windowedHeight=event.resize.h;
            resize(this->windowedWidth, this->windowedHeight);
            break;                      
            
        // Quit event (for example sent when the window is closed)
        case SDL_QUIT:
            this->done=true;
            break;
                
        default:
            break;
    }
}

void Application::storeFrame()
{
	std::string fileName="../images/render";
	char iChar[]={0, 0, 0, 0};

	sprintf(iChar, "%d", this->cntFrame);
	if (this->cntFrame<10)
		fileName+="0";
	if (this->cntFrame<100)
		fileName+="0";
	if (this->cntFrame<1000)
		fileName+="0";
	fileName+=iChar;
	fileName+=".ppm";
	saveFrameBufferPPM(fileName, this->width, this->height);
	//std::cout<<fileName<<" stored"<<std::endl;
}

// Clears the current frame buffer (the image on screen) 
// draws the scene in the other available frame buffer (double buffering)
// and prints the new filled frame buffer on screen
void Application::renderFrame()
{
	// Clears the window with current clearing color, clears also the depth buffer
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	
	// Draws this->scene
	if (this->scene!=NULL)
    	this->scene->drawObjectsOfScene();
    
    if(bStoreFrame)	
		this->storeFrame();
    
    // Performs the buffer swap between the current shown buffer, 
    // and the one we just worked on
    SDL_GL_SwapBuffers();
    
    // Gets the average frame duration over 10 frames
    this->cntFrame++;
    if (this->cntFrame%20==0)
    {
        uint64_t time=getTime();
        this->frameDuration=(time-this->lastStartTime)/20LL;
        this->lastStartTime=time;
    }
    
    // Reports any possible glError
    //std::cout<<"renderFrame error"<<std::endl;
    printGlErrors();
}


// Listens to events during the whole time of the application
// and distributes corresponding tasks
void Application::loop()
{
    SDL_Event event;
    
    // While the application is running (this->done==false), 
    // it waits for events and catch them in "event" when appearing
    // With SDL_wait event, you must render in a special timer event
    // Framarate is then reguled
    //while ((!this->done) && (SDL_WaitEvent(&event)))
    //    handleEvent(event);
    
    // If you expect to render as often as possible, 
    // call SDL_pollEvent between each rendered frame
    while (!this->done)
    {
        while (SDL_PollEvent(&event))
            handleEvent(event);
           
        if(this->simulation != NULL)
        {
		    this->simulation->update();
			this->simulation->render();
		}
		renderFrame();
    }    
        
}


// Moves with a FPS camera
void Application::moveFPS()
{
    // Motion motion information retrieval
    // Comment if mouse motion events are enabled
    SDL_PumpEvents();
    int mouseRelX, mouseRelY;
	#ifdef __APPLE__
        int mystery=0;
        SDL_GetRelativeMouseState(mystery, &mouseRelX, &mouseRelY);
    #else
        SDL_GetRelativeMouseState(&mouseRelX, &mouseRelY);
    #endif
    
	this->xMousePosition+=2.0f*mouseRelX/(GLfloat)this->width;
	this->yMousePosition+=-2.0f*mouseRelY/(GLfloat)this->height;



	GLfloat moveStep=0.02f;
	GLfloat cameraNewPos[3];
	GLfloat moveOnX=this->moveFlags[0]*moveStep;
	GLfloat moveOnY=this->moveFlags[1]*moveStep;
	GLfloat moveOnZ=this->moveFlags[2]*moveStep;
	for (GLuint iCoord=0 ; iCoord<3 ; iCoord++)
	{
		cameraNewPos[iCoord]=this->scene->camera->c[iCoord]
							+this->scene->camera->x[iCoord]*moveOnX
							+this->scene->camera->y[iCoord]*moveOnY
							+this->scene->camera->z[iCoord]*moveOnZ;
	}

	GLfloat angleForWindowWidth=(float)M_PI;
	GLfloat angleForWindowHeight=(float)M_PI/2.0;
	GLfloat angleLong=this->xMousePosition*angleForWindowWidth;
	GLfloat angleLat=this->yMousePosition*angleForWindowHeight;


	// Method with rotates

	GLfloat xAxis[]={1.0, 0.0, 0.0}; GLfloat yAxis[]={0.0, 1.0, 0.0};
	GLfloat rotateAroundX[16]; setToRotate(rotateAroundX, -angleLat, xAxis);
	GLfloat rotateAroundY[16]; setToRotate(rotateAroundY, angleLong, yAxis);
	GLfloat t[]={-cameraNewPos[0], -cameraNewPos[1], -cameraNewPos[2]};
	GLfloat translate[16]; setToTranslate(translate, t);

	setToIdentity(this->scene->camera->view);
	multMatrixBtoMatrixA(this->scene->camera->view, rotateAroundX);
	multMatrixBtoMatrixA(this->scene->camera->view, rotateAroundY);
	multMatrixBtoMatrixA(this->scene->camera->view, translate);

	for (GLuint iCoord=0 ; iCoord<3 ; iCoord++)
	{
		// Updates the axis with values in view
		this->scene->camera->x[iCoord]=this->scene->camera->view[iCoord*4+0];
		this->scene->camera->y[iCoord]=this->scene->camera->view[iCoord*4+1];
		this->scene->camera->z[iCoord]=this->scene->camera->view[iCoord*4+2];
		// Updates the position of the camera c
		this->scene->camera->c[iCoord]=cameraNewPos[iCoord];
	}
}


// Moves with a spherical camera
void Application::moveSpherical()
{
    // Motion motion information retrieval
    // Comment if mouse motion events are enabled
    SDL_PumpEvents();
    int mouseX, mouseY, mouseRelX, mouseRelY;
	#ifdef __APPLE__
        int mystery=0;
        SDL_GetRelativeMouseState(mystery, &mouseRelX, &mouseRelY);
    #else
        SDL_GetRelativeMouseState(&mouseRelX, &mouseRelY);
    #endif
	if (SDL_GetMouseState(&mouseX, &mouseY)&SDL_BUTTON(1))
    {	
	    this->xMouseLeftDownPosition+=2.0f*mouseRelX/(GLfloat)this->width;
	    this->yMouseLeftDownPosition+=-2.0f*mouseRelY/(GLfloat)this->height;
	}



    GLfloat moveStep=0.02f;
    GLfloat cameraNewPos[3];
    GLfloat moveOnX=this->moveFlags[0]*moveStep;
    GLfloat moveOnY=this->moveFlags[1]*moveStep;
    GLfloat moveOnZ=this->moveFlags[2]*moveStep;

    for (GLuint iCoord=0 ; iCoord<3 ; iCoord++)
	{
		cameraNewPos[iCoord]=this->scene->camera->c[iCoord]
							+this->scene->camera->x[iCoord]*moveOnX
							+this->scene->camera->y[iCoord]*moveOnY
							+this->scene->camera->z[iCoord]*moveOnZ;
	}
    
                             
    GLfloat angleForWindowWidth=(GLfloat)M_PI;
    GLfloat angleForWindowHeight=(GLfloat)M_PI/2.0;
    GLfloat angleLong=this->xMouseLeftDownPosition*angleForWindowWidth;
	GLfloat angleLat=this->yMouseLeftDownPosition*angleForWindowHeight;
	GLfloat radius=1.0f+scroll*0.5f; 
	

    // Method with rotates

    GLfloat xAxis[]={1.0f, 0.0f, 0.0f}; GLfloat yAxis[]={0.0f, 1.0f, 0.0f};
    GLfloat rotateAroundX[16]; setToRotate(rotateAroundX, -angleLat, xAxis);
    GLfloat rotateAroundY[16]; setToRotate(rotateAroundY, angleLong, yAxis);
    GLfloat tC[]={-cameraNewPos[0], -cameraNewPos[1], -cameraNewPos[2]};
    GLfloat translateC[16]; setToTranslate(translateC, tC);
    GLfloat tRadius[]={0.0f, 0.0f, -radius};
    GLfloat translateRadius[16]; setToTranslate(translateRadius, tRadius);    
    
    setToIdentity(this->scene->camera->view);
    multMatrixBtoMatrixA(this->scene->camera->view, translateRadius);
    multMatrixBtoMatrixA(this->scene->camera->view, rotateAroundX);
    multMatrixBtoMatrixA(this->scene->camera->view, rotateAroundY);
    multMatrixBtoMatrixA(this->scene->camera->view, translateC);

    for (GLuint iCoord=0 ; iCoord<3 ; iCoord++)
    {
        // Updates the axis with values in view
        this->scene->camera->x[iCoord]=this->scene->camera->view[iCoord*4+0];
        this->scene->camera->y[iCoord]=this->scene->camera->view[iCoord*4+1];
        this->scene->camera->z[iCoord]=this->scene->camera->view[iCoord*4+2];
        // Updates the position of the camera c
        this->scene->camera->c[iCoord]=cameraNewPos[iCoord];
    }
}


// Animates an object
void Application::animate()
{
    this->cntMove++;

	/*GLfloat r=2.0; GLfloat s[]={r, r, r};

    GLfloat S[16]; setToScale(S, s);

    GLfloat d=10.0; GLfloat t[]={0.0, d, 0.0};
    GLfloat T[16]; setToTranslate(T, t);
    
    GLfloat tFrame=1.0;
    GLfloat dt=M_PI/(50.0*tFrame);
    GLfloat a=this->cntMove*dt;
    GLfloat axis[]={0.0, 0.0, 1.0};
    GLfloat R[16]; setToRotate(R, a, axis); 

    GLfloat modelSun[16]; 
    setToIdentity(modelSun);
    multMatrixBtoMatrixA(modelSun, R);
    multMatrixBtoMatrixA(modelSun, T);
    multMatrixBtoMatrixA(modelSun, S);

    GLuint sunID=2;
    this->scene->setDrawnObjectModel(sunID, modelSun);
    GLfloat lightPosition[]={d*modelSun[4], d*modelSun[5], d*modelSun[6], 1.0};
    this->scene->setLight(lightPosition, 1.0);*/    
}



