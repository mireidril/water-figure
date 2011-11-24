// Tools.hpp
// Template Fluids Simulation
// N. Dommanget dommange@univ-mlv.fr


#ifndef __TOOLS_HPP__
#define __TOOLS_HPP__

// To include GL
#include "GLHeaders.hpp"

// C++
#include <string>
#include <vector>

#include <time.h>
#include <stdint.h>

// Mersenne Twister random generation
#include "../api/random/MersenneTwister.h"

GLfloat getNorm (GLfloat * a);
void normalize (GLfloat * a);
void vectorProduct (GLfloat * a, GLfloat * b, GLfloat * result);
void normalFace(GLfloat * O, GLfloat * A, GLfloat * B, GLfloat * normal, bool toNormalize);

void setToIdentity(GLfloat * matrix);
void setToTranslate(GLfloat * matrix, GLfloat * t);
void setToScale(GLfloat * matrix, GLfloat * s);
void setToRotate(GLfloat * matrix, GLfloat angle, GLfloat * axis);
void setPerspective(GLfloat * mat, GLfloat l, GLfloat r, GLfloat b, GLfloat t, GLfloat n, GLfloat f);

void multMatrixBtoMatrixA(GLfloat * A, GLfloat * B);
double getRand(MTRand & mtrand, bool uniform, double min, double max);
void linearInterpolation(int nbCoords, GLfloat * x0, GLfloat * x1, GLfloat toX0, GLfloat * result);
void biLinearInterpolation(int nbCoords, GLfloat * x0y0, GLfloat * x1y0, GLfloat *x0y1, GLfloat *x1y1, GLfloat toX0, GLfloat toY0, GLfloat * result);
void triLinearInterpolation(int nbCoords, GLfloat * x0y0z0, GLfloat * x1y0z0, GLfloat *x0y1z0, GLfloat *x1y1z0, GLfloat * x0y0z1, GLfloat * x1y0z1, GLfloat *x0y1z1, GLfloat *x1y1z1, GLfloat toX0, GLfloat toY0, GLfloat toZ0, GLfloat * result);

void printVec2(GLfloat * vect);
void printVec3(GLfloat * vect);
void printVec4(GLfloat * vect);
void printMat16(GLfloat * mat);
void printGlErrors();

GLuint loadTexture(const std::string &fileName);

std::string * loadFile(const std::string & fileName);
void printShaderLog(GLint shaderId);
GLuint loadProgram(const std::vector<std::string> & files);
void setMatricesInShader(GLuint shaderId, GLfloat * model, GLfloat * view, GLfloat * c, GLfloat * projection);
void setLightInShader(GLuint shaderID, GLfloat * position, GLfloat power);
void setMaterialInShader(GLuint shaderID, GLfloat * ambient, GLfloat * diffuse, GLfloat * specular, GLfloat ka, GLfloat kd, GLfloat ks, GLfloat shininess);
void changeMaterialColorInShader(GLuint shaderID, GLfloat * color);
void setFilledDataInShader(GLuint shaderID, GLboolean positions, GLboolean normals, GLboolean uvs, GLboolean colors);
void setTextureUnitsInShader(GLuint shaderID);

unsigned char * loadPPM(const std::string & filename, GLuint * width, GLuint * height);

uint64_t getTime();

void saveFrameBufferPPM(const std::string & fileName, unsigned int _width, unsigned int _height);

#endif //  __TOOLS_HPP__
