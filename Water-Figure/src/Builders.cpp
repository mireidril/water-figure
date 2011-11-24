// Builders.cpp
// Template Fluids Simulation
// N. Dommanget dommange@univ-mlv.fr


#include "Builders.hpp"
#include "Object.hpp"
#include "Tools.hpp"

#include <iostream>
#define _USE_MATH_DEFINES
#include <math.h>
#include <fstream>
#include <sstream>
#include <cstdlib>


//______________________________________________________________________________
// Building functions


// Build axis
void buildAxis(Object * object)
{
    std::cout<<"       - axis"<<std::endl;

    object->nbVertices=4;
    object->nbIndices=6;
    
    // The vertices of the axis
	GLfloat vertices[]={0.0, 0.0, 0.0, 1.0,  // O
	                    1.0, 0.0, 0.0, 1.0,  // x
	                    0.0, 1.0, 0.0, 1.0,  // y
	                    0.0, 0.0, 1.0, 1.0}; // z
	                    	
	// The indices of the vertices of the 3 lines          
	GLuint indices[]={0, 1,
                      0, 2,
                      0, 3};
    GLfloat colors[]={1.0, 1.0, 1.0, 1.0,  // white
	                  1.0, 0.0, 0.0, 1.0,  // red
	                  0.0, 1.0, 0.0, 1.0,  // green
	                  0.0, 0.0, 1.0, 1.0}; // blue
	
    // Sends the data into buffers on the GPU
    object->sendPrimitives(vertices, indices);
    object->sendColors(colors);                
}


// Building one square with Offset
void buildEarthPlane(Object * object, GLfloat width, GLfloat length)
{
    std::cout<<"       - earth plane"<<std::endl;

    object->nbVertices=4;
    object->nbIndices=6;

    // The 3 vertices of a triangle
    GLfloat w=width/2.0f;
    GLfloat l=length/2.0f;
	GLfloat vertices[]={-w, 0.0f, l, 1.0f,
	                     w, 0.0f, l, 1.0f,
	                    -w, 0.0f, -l, 1.0f,
	                     w, 0.0f, -l, 1.0f};

    // The 3 indices of the vertices to draw the face
    GLuint indices[]={0, 2, 3, // First triangle : top left
    				  0, 1, 3}; // Second triangle : bottom right

    // One simple normal
	GLfloat normals[]={0.0f, 1.0f, 0.0f,
	                   0.0f, 1.0f, 0.0f,
	                   0.0f, 1.0f, 0.0f,
	                   0.0f, 1.0f, 0.0f};

    // Sends the data into buffers on the GPU
    object->sendPrimitives(vertices, indices);
    object->sendNormals(normals);
}





// Builds sphere with triangles, smooth normals (method 1) with no redondancy
void buildSphere_TrSmoothNonRed(Object * object, GLfloat radius, GLint discLat, GLint discLong)
{
    std::cout<<"       - sphere with triangles, smooth and non redondant"<<std::endl;

    object->nbVertices=discLong*(discLat+1);
    object->nbIndices=(discLong*discLat)*6;
    
    GLfloat deltaLong= (float) (2.0f*M_PI)/GLfloat(discLong);
    GLfloat deltaLat=(float)M_PI/GLfloat(discLat);
    
    GLfloat localRadius=0.0f;
    GLfloat bottom=-1.0f;
    GLfloat localHeight=bottom;
    
    GLfloat * vertices = new GLfloat[object->nbVertices*4];
    GLfloat * normals = new GLfloat[object->nbVertices*3];
    GLuint * indices = new GLuint[object->nbIndices];

    
    // Index of the first vertex of a pair of triangles : the BL
    int indexBL=0;
    // Index of this vertex (BL) in the indices table, to draw two triangle from it
    int indexBLIndices=0;
    for (int iLat=0 ; iLat<(discLat+1) ; iLat++)
    {
        int indexBLFirstOfRow=indexBL;
        for (int iLong=0 ; iLong<discLong ; iLong++)
        {
            GLfloat latB=(iLat*deltaLat)-((float)M_PI/2.0f);  // Latitude Bottom
            GLfloat longL=iLong*deltaLong;  // Longitude Left
            
            localRadius=cos(latB); 
            localHeight=sin(latB);
            GLfloat BL[4]={localRadius*cos(longL), localHeight, -localRadius*sin(longL), 1.0}; // Vertex Bottom Left
            
            for (int iCoord=0 ; iCoord<4 ; iCoord++)            
            {
                GLfloat coef=radius;
                if (iCoord==3)
                    coef=1.0;
                vertices[(indexBL*4)+iCoord]=coef*BL[iCoord];
                
                if (iCoord<3)
                    normals[(indexBL*3)+iCoord]=BL[iCoord];
            }
            
            if (iLat<discLat)
            {
                indices[indexBLIndices+0]=indexBL;                 // BL
                indices[indexBLIndices+1]=indexBL + 1;             // BR
                indices[indexBLIndices+2]=indexBL + 1 + discLong;  // TR 
                indices[indexBLIndices+3]=indexBL;                 // BL
                indices[indexBLIndices+4]=indexBL + 1 + discLong;  // TR
                indices[indexBLIndices+5]=indexBL     + discLong;  // TL
            }
            
            if (iLong==(discLong-1))
            {
                indices[indexBLIndices+0]=indexBL;                 // BL
                indices[indexBLIndices+1]=indexBLFirstOfRow;             // BR
                indices[indexBLIndices+2]=indexBLFirstOfRow + discLong;  // TR 
                indices[indexBLIndices+3]=indexBL;                 // BL
                indices[indexBLIndices+4]=indexBLFirstOfRow + discLong;  // TR
                indices[indexBLIndices+5]=indexBL     + discLong;  // TL            
            }
            
            indexBL++;
            indexBLIndices+=6;
        }
    }

    // Sends the data into buffers on the GPU
    object->sendPrimitives(vertices, indices);
    //object->sendNormals(normals);
}



// Builds Wind Rose
void buildWindRose(Object * object)
{
    std::cout<<"       - wind rose"<<std::endl;

    object->nbVertices=24*2;
    object->nbIndices=24*2;
    
    GLfloat z=0.1f;
    GLfloat l=0.1f;
    GLfloat L=0.5f;
    
    GLfloat Af[]= {0.0f,  0.0f,   z,  1.0f};
    GLfloat Ab[]= {0.0f,  0.0f,  -z,  1.0f};
    GLfloat B[]= {  L,  0.0f,  0.0f,  1.0f};
    GLfloat bc[]={  l,    l,  0.0f,  1.0f};
    GLfloat C[]= {0.0f,    L,  0.0f,  1.0f};
    GLfloat cd[]={- l,    l,  0.0f,  1.0f};
    GLfloat D[]= {- L,  0.0,  0.0f,  1.0f};
    GLfloat de[]={- l,  - l,  0.0f,  1.0f};
    GLfloat E[]= {0.0f,  - L,  0.0f,  1.0f};
    GLfloat eb[]={  l,  - l,  0.0f,  1.0f};
    
    // The 3 vertices of a triangle
    GLfloat vertices[]={Af[0], Af[1], Af[2], Af[3], //Af
                        C[0], C[1], C[2], C[3], //C
                        cd[0], cd[1], cd[2], cd[3], //cd
                        
                        Af[0], Af[1], Af[2], Af[3], //Af
                        cd[0], cd[1], cd[2], cd[3], //cd
                        D[0], D[1], D[2], D[3], //D
                        
                        Af[0], Af[1], Af[2], Af[3], //Af
                        D[0], D[1], D[2], D[3], //D
                        de[0], de[1], de[2], de[3], //de
                        
                        Af[0], Af[1], Af[2], Af[3], //Af
                        de[0], de[1], de[2], de[3], //de
                        E[0], E[1], E[2], E[3], //E
                         
                        Af[0], Af[1], Af[2], Af[3], //Af
                        E[0], E[1], E[2], E[3], //E
                        eb[0], eb[1], eb[2], eb[3], //eb  
                               
                        Af[0], Af[1], Af[2], Af[3], //Af 
                        eb[0], eb[1], eb[2], eb[3], //eb
                        B[0], B[1], B[2], B[3], //B
                           
                        Af[0], Af[1], Af[2], Af[3], //Af
                        B[0], B[1], B[2], B[3], //B
                        bc[0], bc[1], bc[2], bc[3], //bc
                           
                        Af[0], Af[1], Af[2], Af[3], //Af
                        bc[0], bc[1], bc[2], bc[3], //bc
                        C[0], C[1], C[2], C[3], //C   
                         
                         
                         
                        Ab[0], Ab[1], Ab[2], Ab[3], //Ab
                        cd[0], cd[1], cd[2], cd[3], //cd
                        C[0], C[1], C[2], C[3], //C
                        
                        Ab[0], Ab[1], Ab[2], Ab[3], //Ab
                        D[0], D[1], D[2], D[3], //D
                        cd[0], cd[1], cd[2], cd[3], //cd
                        
                        Ab[0], Ab[1], Ab[2], Ab[3], //Ab
                        de[0], de[1], de[2], de[3], //de
                        D[0], D[1], D[2], D[3], //D
                        
                        Ab[0], Ab[1], Ab[2], Ab[3], //Ab
                        E[0], E[1], E[2], E[3], //E
                        de[0], de[1], de[2], de[3], //de
                         
                        Ab[0], Ab[1], Ab[2], Ab[3], //Ab
                        eb[0], eb[1], eb[2], eb[3], //eb
                        E[0], E[1], E[2], E[3], //E
                               
                        Ab[0], Ab[1], Ab[2], Ab[3], //Ab 
                        B[0], B[1], B[2], B[3], //B
                        eb[0], eb[1], eb[2], eb[3], //eb
                           
                        Ab[0], Ab[1], Ab[2], Ab[3], //Ab
                        bc[0], bc[1], bc[2], bc[3], //bc
                        B[0], B[1], B[2], B[3], //B
                           
                        Ab[0], Ab[1], Ab[2], Ab[3], //Ab
                        C[0], C[1], C[2], C[3], //C
                        bc[0], bc[1], bc[2], bc[3]}; //bc
				
    // The 3 indices of the vertices to draw the face
    GLuint * indices = new GLuint[object->nbIndices];
    for (GLuint iIndices=0 ; iIndices<object->nbIndices ; iIndices++)
        indices[iIndices]=iIndices;

	GLfloat normals[3*(24*2)];
    setNormalsFlatTr(object, vertices, indices,  normals);


    GLfloat * uvs = new GLfloat[object->nbVertices*2];
    for (GLuint iVertices=0 ; iVertices<object->nbVertices ; iVertices++)
    {
        uvs[iVertices*2+0]=vertices[iVertices*4+0]+L;
        uvs[iVertices*2+1]=vertices[iVertices*4+1]+L;
    }

	GLfloat * colors = new GLfloat[object->nbVertices*4];
    for (GLuint iVertices=0 ; iVertices<object->nbVertices ; iVertices++)
    {
        for (GLuint iCoord=0 ; iCoord<3 ; iCoord++)
            colors[iVertices*4+iCoord]=0.0;
        colors[iVertices*4+3]=1.0;
    }   
    colors[1*4+0]=1.0;
    colors[23*4+0]=1.0;
    colors[26*4+0]=1.0;
    colors[46*4+0]=1.0;

    // Sends the data into buffers on the GPU
    object->sendPrimitives(vertices, indices);
    object->sendNormals(normals);
    object->sendUvs(uvs);
    object->sendColors(colors);
}


// Builds a cube
void buildCube(Object * object)
{
    std::cout<<"       - cube"<<std::endl;

    object->nbVertices=4*6;
    object->nbIndices=6*6;

    // The 3 vertices of a triangle
    GLfloat L=1.0;
    GLuint i=0;
    GLfloat A[]={ L,  L,  L, 1.0}; normalize(A);
    GLfloat B[]={ L,  L, -L, 1.0}; normalize(B);
    GLfloat C[]={-L,  L, -L, 1.0}; normalize(C);
    GLfloat D[]={-L,  L,  L, 1.0}; normalize(D);
    GLfloat E[]={ L, -L,  L, 1.0}; normalize(E);
    GLfloat F[]={ L, -L, -L, 1.0}; normalize(F);
    GLfloat G[]={-L, -L, -L, 1.0}; normalize(G);
    GLfloat H[]={-L, -L,  L, 1.0}; normalize(H);
    
    GLfloat  X[]={1.0, 0.0, 0.0};
    GLfloat mX[]={-1.0, 0.0, 0.0};
    GLfloat  Y[]={0.0, 1.0, 0.0};
    GLfloat mY[]={0.0, -1.0, 0.0};
    GLfloat  Z[]={0.0, 0.0, 1.0};
    GLfloat mZ[]={0.0, 0.0, -1.0};
         
    GLuint a0=i++; GLuint a1=i++; GLuint a2=i++;
    GLuint b0=i++; GLuint b1=i++; GLuint b2=i++;
    GLuint c0=i++; GLuint c1=i++; GLuint c2=i++;
    GLuint d0=i++; GLuint d1=i++; GLuint d2=i++;
    GLuint e0=i++; GLuint e1=i++; GLuint e2=i++;
    GLuint f0=i++; GLuint f1=i++; GLuint f2=i++;  
    GLuint g0=i++; GLuint g1=i++; GLuint g2=i++;     
    GLuint h0=i++; GLuint h1=i++; GLuint h2=i++;   
    
    
	GLfloat vertices[]=
	{A[0], A[1], A[2], A[3],   A[0], A[1], A[2],   A[3],A[0], A[1], A[2], A[3],
	 B[0], B[1], B[2], B[3],   B[0], B[1], B[2],   B[3],B[0], B[1], B[2], B[3],
	 C[0], C[1], C[2], C[3],   C[0], C[1], C[2],   C[3],C[0], C[1], C[2], C[3],
	 D[0], D[1], D[2], D[3],   D[0], D[1], D[2],   D[3],D[0], D[1], D[2], D[3],
	 E[0], E[1], E[2], E[3],   E[0], E[1], E[2],   E[3],E[0], E[1], E[2], E[3],
	 F[0], F[1], F[2], F[3],   F[0], F[1], F[2],   F[3],F[0], F[1], F[2], F[3],
	 G[0], G[1], G[2], G[3],   G[0], G[1], G[2],   G[3],G[0], G[1], G[2], G[3],
	 H[0], H[1], H[2], H[3],   H[0], H[1], H[2],   H[3],H[0], H[1], H[2], H[3]};
			    
	GLfloat normals[]=
	{Y[0], Y[1], Y[2],      X[0], X[1], X[2],      Z[0], Z[1], Z[2],     //A
	 Y[0], Y[1], Y[2],      X[0], X[1], X[2],      mZ[0], mZ[1], mZ[2],  //B
	 Y[0], Y[1], Y[2],      mX[0], mX[1], mX[2],   mZ[0], mZ[1], mZ[2],  //C
	 Y[0], Y[1], Y[2],      mX[0], mX[1], mX[2],   Z[0], Z[1], Z[2],     //D
	 mY[0], mY[1], mY[2],   X[0], X[1], X[2],      Z[0], Z[1], Z[2],     //E
	 mY[0], mY[1], mY[2],   X[0], X[1], X[2],      mZ[0], mZ[1], mZ[2],  //F
	 mY[0], mY[1], mY[2],   mX[0], mX[1], mX[2],   mZ[0], mZ[1], mZ[2],  //G
	 mY[0], mY[1], mY[2],   mX[0], mX[1], mX[2],   Z[0], Z[1], Z[2] };   //H		    
			    
			    
	GLfloat uvs3[]={A[0], A[1], A[2],   A[0], A[1], A[2],   A[0], A[1], A[2],
			        B[0], B[1], B[2],   B[0], B[1], B[2],   B[0], B[1], B[2],
			        C[0], C[1], C[2],   C[0], C[1], C[2],   C[0], C[1], C[2],
			        D[0], D[1], D[2],   D[0], D[1], D[2],   D[0], D[1], D[2],
			        E[0], E[1], E[2],   E[0], E[1], E[2],   E[0], E[1], E[2],
			        F[0], F[1], F[2],   F[0], F[1], F[2],   F[0], F[1], F[2],
			        G[0], G[1], G[2],   G[0], G[1], G[2],   G[0], G[1], G[2],
			        H[0], H[1], H[2],   H[0], H[1], H[2],   H[0], H[1], H[2]};
			    
	GLuint indices[]={ a0, c0, b0,
				       a0, d0, c0,
				       a1, b1, f1,
				       a1, f1, e1,
				       a2, e2, h2,
				       a2, h2, d2,
				       b2, c2, g2,
				       b2, g2, f2,
				       c1, d1, g1,
				       d1, h1, g1,
				       h0, e0, f0,
				       h0, f0, g0};

    // Sends the data into buffers on the GPU
    object->sendPrimitives(vertices, indices);
    object->sendNormals(normals);
    object->sendUvs(uvs3);
}

//______________________________________________________________________________
// Mesh conformation functions


// To fill the normals for vertices in flat triangles
void setNormalsFlatTr(Object * object, GLfloat * vertices, GLuint * indices, GLfloat * normals)
{
	for (GLuint iTriangle=0 ; iTriangle<(object->nbIndices) ; iTriangle+=3)
	{
		GLuint iIndices0=indices[iTriangle+0];
		GLuint iIndices1=indices[iTriangle+1];
		GLuint iIndices2=indices[iTriangle+2];
		GLuint iVertices0=iIndices0*4;
		GLuint iVertices1=iIndices1*4;
		GLuint iVertices2=iIndices2*4;
		GLuint iNormals0=iIndices0*3;
		GLuint iNormals1=iIndices1*3;
		GLuint iNormals2=iIndices2*3;
		
		normalFace(&(vertices[iVertices0]), &(vertices[iVertices1]), &(vertices[iVertices2]), &(normals[iNormals0]), true);
		normals[iNormals1+0]=normals[iNormals0+0]; 
		normals[iNormals1+1]=normals[iNormals0+1]; 
		normals[iNormals1+2]=normals[iNormals0+2];
		normals[iNormals2+0]=normals[iNormals0+0]; 
		normals[iNormals2+1]=normals[iNormals0+1]; 
		normals[iNormals2+2]=normals[iNormals0+2];
	}
}


// Centers and normalize any mesh
void centerAndNormalizeMesh(Object * object, GLfloat * vertices)
{
    GLfloat inf=100000.0;
    // Normalizes to 1.0 the size of the mesh
    GLfloat min[]={inf, inf, inf};
    GLfloat max[]={-inf, -inf, -inf};

    for(GLuint iVerticesValues=0 ; iVerticesValues<object->nbVertices ; iVerticesValues++)
    {
        for(GLuint iCoord=0 ; iCoord<3 ; iCoord++)
        {
            GLfloat value=vertices[iVerticesValues*4+iCoord];
            if (value<min[iCoord])
                min[iCoord]=value;
            if (value>max[iCoord])
                max[iCoord]=value;
        }
    }
    GLfloat maxSize[]={max[0]-min[0], max[1]-min[1], max[2]-min[2]};
    GLuint maxCoord=0;
    if ((maxSize[1]>=maxSize[0]) && (maxSize[1]>=maxSize[2])) maxCoord=1;    
    if ((maxSize[2]>=maxSize[0]) && (maxSize[2]>=maxSize[1])) maxCoord=2;

    for(GLuint iVerticesValues=0 ; iVerticesValues<object->nbVertices ; iVerticesValues++)
    {
        for(GLuint iCoord=0 ; iCoord<3 ; iCoord++)
        {
            vertices[iVerticesValues*4+iCoord]-=min[iCoord]+(maxSize[iCoord]/2.0f);
            vertices[iVerticesValues*4+iCoord]/=maxSize[maxCoord];
        }
    } 
}






//______________________________________________________________________________
// Obj file reading fuctions


// Fills a std::vector with strings from cuts at every "delim" in "string"
void split(const std::string& string, std::vector<std::string>& tokens, const std::string& delim)
{
        // Clears the vector where the results are put
        tokens.clear();
        // Goes through the string and extract the tokens
        for (std::string::size_type p1=0, p2=0; p1!=std::string::npos; )
        {
                p1=string.find_first_not_of(delim, p1);
                if (p1!=std::string::npos)
                {
                        p2=string.find_first_of(delim, p1);
                        tokens.push_back(string.substr(p1, p2-p1));
                        p1=p2;
                }
        }
}

// Fills a std::vector with three GLfloat values in a string
void readVec3(std::istringstream& line, std::vector<GLfloat> * vertices)
{
    GLfloat vec[3]={0.0, 0.0, 0.0};
    line >> vec[0] >> vec[1] >> vec[2];
    vertices->push_back(vec[0]);
    vertices->push_back(vec[1]);
    vertices->push_back(vec[2]);
}

// Fills a std::vector with two GLfloat values in a string
void readVec2(std::istringstream& line, std::vector<GLfloat> * vertices)
{
    GLfloat vec[2]={0.0, 0.0};
    line >> vec[0] >> vec[1];
    vertices->push_back(vec[0]);
    vertices->push_back(1.0f-vec[1]);
}

// Fills std::vectors with indices of three type for the vertices in the face
void readFace(std::istringstream& line, std::vector<GLuint> * indices,  std::vector<GLuint> * uvIndices, std::vector<GLuint> * normalsIndices)
{
    std::string face;
    while (!line.eof())
	{
        line >> face;
            std::vector<std::string> sequence;
            split(face, sequence, "/");
            size_t index=strtoul(sequence[0].c_str(), NULL, 10)-1;
            indices->push_back(index);

            size_t uvIndex=strtoul(sequence[1].c_str(), NULL, 10)-1;
            uvIndices->push_back(uvIndex);

            size_t normalIndex=strtoul(sequence[2].c_str(), NULL, 10)-1;
            normalsIndices->push_back(normalIndex);
    }
}

// Regroup in four arrays instead of six by assigning a uv coordinates and a normal to the corresponding vertex. 
// Does not create new vertices : works only if one normal is sufficient per location in space : smoothed objects
void reorderUvsAndNormalsIfSmooth(std::vector<GLfloat> * vertices, std::vector<GLfloat> * uvs, std::vector<GLfloat> * normals,
                          std::vector<GLuint> * indices, std::vector<GLuint> * uvIndices, std::vector<GLuint> * normalIndices)
{
    std::vector<GLfloat> newUvs(vertices->size(), 0.0);
    std::vector<GLfloat> newNormals(vertices->size(), 0.0);
    GLuint uvIndex;
    GLuint normalIndex;
    for(GLuint iIndices=0 ; iIndices<indices->size() ; iIndices++)
    {
        uvIndex=(*uvIndices)[iIndices];
        normalIndex=(*normalIndices)[iIndices];
        for (GLuint iCoord=0 ; iCoord<3 ; iCoord++)
        {
            newNormals[((*indices)[iIndices]*3)+iCoord]=(*normals)[(normalIndex*3)+iCoord];
            if (iCoord<2)
                newUvs[((*indices)[iIndices]*2)+iCoord]=(*uvs)[(uvIndex*2)+iCoord];
        }
    }
    
    uvs->clear();
    normals->clear();
    for(GLuint iVertices=0 ; iVertices<vertices->size() ; iVertices++)
    {
        normals->push_back(newNormals[iVertices]);
    }
    for(GLuint iUvs=0 ; iUvs<newUvs.size() ; iUvs++)
    {      
        uvs->push_back(newUvs[iUvs]);
    }
}

// Regroup in four arrays instead of six by assigning a uv coordinates and a normal to corresponding vertices. 
// Creates new vertices to be able to give several normals per location in space : non smoothed objects
void reorderUvsAndNormalsIfNonSmooth(std::vector<GLfloat> * vertices, std::vector<GLfloat> * uvs, std::vector<GLfloat> * normals,
                                     std::vector<GLuint> * indices, std::vector<GLuint> * uvIndices, std::vector<GLuint> * normalIndices)
{
    std::vector<GLfloat> newVertices;
    std::vector<GLfloat> newUvs;
    std::vector<GLfloat> newNormals;
    GLuint vertexIndex;
    GLuint uvIndex;
    GLuint normalIndex;
    for(GLuint iIndices=0 ; iIndices<indices->size() ; iIndices++)
    {
        vertexIndex=(*indices)[iIndices];
        uvIndex=(*uvIndices)[iIndices];
        normalIndex=(*normalIndices)[iIndices];
        for (GLuint iCoord=0 ; iCoord<3 ; iCoord++)
        {
            newVertices.push_back((*vertices)[(vertexIndex*3)+iCoord]);
            newNormals.push_back((*normals)[(normalIndex*3)+iCoord]);
            if (iCoord<2)
                newUvs.push_back((*uvs)[(uvIndex*2)+iCoord]);
        }
        (*indices)[iIndices]=iIndices;
    }
    
    vertices->clear();
    uvs->clear();
    normals->clear();
    for(GLuint iVertices=0 ; iVertices<newVertices.size() ; iVertices++)
    {
        vertices->push_back(newVertices[iVertices]);
        normals->push_back(newNormals[iVertices]);
    }
    for(GLuint iUvs=0 ; iUvs<newUvs.size() ; iUvs++)
    {      
        uvs->push_back(newUvs[iUvs]);
    }
}

// Adds w to every vertices
void addHomogeneousToVertices(std::vector<GLfloat> * vertices)
{
    std::vector<GLfloat> newVertices;
    for(GLuint iVertices=0 ; iVertices<(vertices->size())/3 ; iVertices++)
    {
        for (GLuint iCoord=0 ; iCoord<3 ; iCoord++)
            newVertices.push_back((*vertices)[(iVertices*3)+iCoord]);        
        newVertices.push_back(1.0); 
    }
    
    vertices->clear();
    for(GLuint iVertices=0 ; iVertices<newVertices.size() ; iVertices++)
        vertices->push_back(newVertices[iVertices]);   
}


// Changes the array to conform our Object's vbos formats
void conformToObject(std::vector<GLfloat> * vertices, std::vector<GLfloat> * uvs, std::vector<GLfloat> * normals)
{
    addHomogeneousToVertices(vertices);
    // Normalizes normals
    for(GLuint iNormals=0 ; iNormals<(normals->size())/3 ; iNormals++)
    {
        normalize(&((*normals)[iNormals*3]));
    }
}


// Builds an object made from an OBJ file, taking only geometry into account (not materials)
bool buildObjectGeometryFromOBJ(Object * object, const std::string& fileName, bool smoothObject)
{
    std::ifstream file(fileName.c_str(), std::ios_base::in);
    if(!file)
    {
        std::cout<<"       Error while loading object from .obj file : "<<fileName<<"."<<std::endl;
        return false;
    }
    std::cout<<"       Loading object from .obj file : "<<fileName<<"."<<std::endl; 

    bool hasVt=false;
    bool hasVn=false;
	std::vector<GLfloat> vertices;
    std::vector<GLfloat> uvs;
    std::vector<GLfloat> normals;
    std::vector<GLuint> indices;
    std::vector<GLuint> uvIndices;
    std::vector<GLuint> normalIndices;
    std::string buf, key, name, MTLFileName;

    while (getline(file, buf))
	{
            std::istringstream line(buf);
            line >> key;    
            if (key=="o")
      		    line >> name;
            else if(key=="v")
      		     readVec3(line, &vertices);
      		else if(key == "vt")
            {    readVec2(line, &uvs); hasVt=true; }
            else if(key=="vn")
      		{    readVec3(line, &normals); hasVn=true; }
            else if(key=="f")
      		    readFace(line, &indices,  &uvIndices, &normalIndices);
            else if (key=="mltlib")
      		    line >> MTLFileName;
    }

    std::cout<<"       Obj mesh "<<name<<" loading..."<<std::endl;
    std::cout<<"       Obj meshes should only be made of triangles (for this loader), make sure this is correct in file."<<std::endl;
    if (!hasVt) 
        std::cout<<"       Obj file "<<name<<" has no texture coordinates, add some in modeler."<<std::endl;

    if (!hasVn) 
    {
        std::cout<<"       Obj file "<<name<<" has no normals, add some in modeler."<<std::endl;
    }
    else 
    {
        if (normals.size()*3==(indices.size()*3)) 
        {
            std::cout<<"       Obj file "<<name<<" was not smoothed in modeler."<<std::endl;
            if (smoothObject==true)
                std::cout<<"       Normals will be wrong with smoothObject parameter as true, change it to false."<<std::endl;
        }
    }
    bool onenormalPerTriangle=(normals.size()*3==(indices.size()*3));

    if (smoothObject)
        reorderUvsAndNormalsIfSmooth(&vertices, &uvs, &normals, &indices, &uvIndices, &normalIndices);
    else
        reorderUvsAndNormalsIfNonSmooth(&vertices, &uvs, &normals, &indices, &uvIndices, &normalIndices);
    conformToObject(&vertices, &uvs, &normals);

    object->nbVertices=vertices.size()/4;
    object->nbIndices=indices.size();
    
    // Normalizes to 1.0 the size of the mesh and centers it
    centerAndNormalizeMesh(object, vertices.data());


    object->sendPrimitives(vertices.data(), indices.data());
    if (!hasVt) 
        std::cout<<"       WARNING : Obj file "<<name<<" has no texture coordinates, add some in modeler."<<std::endl;
    else
        object->sendUvs(uvs.data());

    if (!hasVn) 
    {
        std::cout<<"       WARNING : Obj file "<<name<<" has no normals, add some in modeler."<<std::endl;
    }
    else 
    {
        if (onenormalPerTriangle) 
        {
            std::cout<<"       Obj file "<<name<<" was not smoothed in modeler."<<std::endl;
            if (smoothObject)
                std::cout<<"       WARNING : smoothObject==true. Normals will be wrong : change it to false."<<std::endl;
        }
        object->sendNormals(normals.data());
    }
    std::cout<<"       Material files are not taken into account by this loader."<<std::endl;
    
    return true;
}
