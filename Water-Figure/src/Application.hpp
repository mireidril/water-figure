// Application.hpp
// Template Fluids Simulation
// N. Dommanget dommange@univ-mlv.fr


#ifndef __APPLICATION_HPP__
#define __APPLICATION_HPP__

#include "GLHeaders.hpp"
// Windowing system SDL
#include <SDL.h>

// Forward declaration
struct Scene;
struct Simulation;

// All the initialisation of states and events for SDL and OpenGL
class Application
{
    public:
    
        //int renderEventID;
        int moveEventID;
        int animateEventID;
	    SDL_TimerID renderTimer;                // Timer for the rendering
        SDL_TimerID moveTimer;                  // Timer for the moves
        SDL_TimerID animateTimer;               // Timer for the animation     
        
	    SDL_Surface* drawContext;

	    unsigned int windowedWidth;               // Window dimentions when not fullscreen - horizontal axis
	    unsigned int windowedHeight;              // Window dimentions when not fullscreen - vertical axis
	    unsigned int fullScreenWidth;               // Screen dimentions - horizontal axis
	    unsigned int fullScreenHeight;              // Screen dimentions - vertical axis	
	    unsigned int width;                     // Window actual dimentions - horizontal axis
	    unsigned int height;                    // Window actual dimentions - vertical axis
	    
        Uint32 videoModeFlags;
        
	    GLfloat xMousePosition;                 // Mouse position - horizontal axis (-1 to 1)
	    GLfloat yMousePosition;                 // mouse position - vertical axis (-1 to 1)
	    GLfloat xMouseLeftDownPosition;         // Mouse position - updated only when left button down
	    GLfloat yMouseLeftDownPosition;         // mouse position - updated only when left button down
	    int scroll;                             // scroll value (up : ++, down : --)
	    
	    GLfloat moveFlags[4];
	    bool mouse;                             // True if mouse is seeable
	    
	    GLuint cntFrame;                        // Frame counter
	    uint64_t lastStartTime;                 // Time updated every 10 frames
	    uint64_t frameDuration;                 // Frame time duration
	    
	    GLuint cntMove;                         // Move counter
	    bool done;                              // True when the window is closed to end the application
        Scene * scene;                          // Scene to draw
        Simulation * simulation;
        GLuint eventCnt;

        Application();
        ~Application();
        
        void init();
        void initSDLOpenGL();
        void customizeStates();
             
        void setScene(Scene * scene);
        
        void initTimers();

        void switchWireframe();
        void resize(GLuint w, GLuint h);
        void switchFullScreen();
        void switchMouse();
        void setBackgroundColor(GLfloat * color);
        void printFPS();
        
        void handleUserEvent(SDL_Event& event);
        void handleKeyEvent(SDL_keysym& keysym, bool down);
        void handleMouseMotion(SDL_MouseMotionEvent& motion);
        void handleEvent(SDL_Event& event);
        void storeFrame();
        
        void renderFrame();
        
        void loop();

        
        void moveFPS();
        void moveSpherical();
        void animate();
};

Uint32 genericTimer(Uint32 interval, void* param);


#endif //__APPLICATION_HPP__ 
