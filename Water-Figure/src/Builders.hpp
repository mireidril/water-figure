// Builders.hpp
// Template Fluids Simulation
// N. Dommanget dommange@univ-mlv.fr


#ifndef __BUILDERS_HPP__
#define __BUILDERS_HPP__

#include "GLHeaders.hpp"

#include <vector>
#include <string>

//Forward declaration
class Object;

// Building functions
void buildAxis(Object * object);
void buildEarthPlane(Object * object, GLfloat width, GLfloat length);
void buildSphere_TrSmoothNonRed(Object * object, GLfloat radius, GLint discLat, GLint discLong);
void buildWindRose(Object * object);
void buildCube(Object * object);

// Mesh conformation functions
void setNormalsFlatTr(Object * object, GLfloat * vertices, GLuint * indices, GLfloat * normals);
void centerAndNormalizeMesh(Object * object, GLfloat * vertices);


// Obj file reading fuctions
void split(std::string * string, std::vector<std::string> * tokens, std::string& delim);
void readVec3(std::istringstream& line, std::vector<GLfloat> * vertices);
void readVec2(std::istringstream& line, std::vector<GLfloat> * vertices);
void readFace(std::istringstream * line, std::vector<GLuint> * indices,  std::vector<GLuint> * uvIndices, std::vector<GLuint> * normalsIndices);
void reorderUvsAndNormalsIfSmooth(std::vector<GLfloat> * vertices, std::vector<GLfloat> * uvs, std::vector<GLfloat> * normals,
                          std::vector<GLuint> * indices, std::vector<GLuint> * uvIndices, std::vector<GLuint> * normalsIndices);
void reorderUvsAndNormalsIfNonSmooth(std::vector<GLfloat> * vertices, std::vector<GLfloat> * uvs, std::vector<GLfloat> * normals,
                                     std::vector<GLuint> * indices, std::vector<GLuint> * uvIndices, std::vector<GLuint> * normalsIndices);
void addHomogeneousToVertices(std::vector<GLfloat> * vertices);
bool buildObjectGeometryFromOBJ(Object * object, const std::string& fileName, bool smoothObject);


#endif
