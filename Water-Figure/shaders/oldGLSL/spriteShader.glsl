// lightingTexturingShader.glsl
// Template Fluids Simulation
// N. Dommanget dommange@univ-mlv.fr


#ifdef _VERTEX_

// Attributes : per vertex data
attribute vec4 vertexPosition;
attribute vec3 vertexNormal;
attribute vec2 vertexUvs;
attribute vec4 vertexColor;

// Varyings : data to transmit to fragments
varying vec4 position;
varying vec4 normal;
varying vec2 uvs;
varying vec4 localColor;

void main()
{
    if (filledData[0]) position = model * vertexPosition;
    if (filledData[1]) normal = model * vec4(vertexNormal, 0.0);
    if (filledData[2]) uvs = vertexUvs;
    if (filledData[3]) localColor = vertexColor;
    
    gl_PointSize=15;
    
    gl_Position = projection * view * model * vertexPosition;
}

#endif




#ifdef _FRAGMENT_


// Varyings : data receved and interpolated from the vertex shaders
varying vec4 position;
varying vec4 normal;
varying vec2 uvs;
varying vec4 localColor;

void main()
{
    vec4 text=vec4(texture2D(textureUnitDiffuse, gl_PointCoord).rgb, 1.0);
    float energy=((text.r+text.g+text.b)/3.0)-0.1;
    
    fragColor=clamp(energy, 0.0, 1.0)*localColor;
}

#endif
