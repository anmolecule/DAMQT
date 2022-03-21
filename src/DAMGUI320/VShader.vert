 #version 330 core
// #ifdef GL_ES
// // Set default precision to medium
// precision mediump int;
// precision mediump float;
// #endif

uniform mat4 m_matrix;
uniform mat4 v_matrix;
uniform mat4 mv_matrix;
uniform mat4 mvp_matrix;

uniform vec3 LightPosition_worldspace;
uniform vec3 Ambient_Color;
uniform vec3 Light_Color;
uniform vec3 Specular_Color;
uniform float Light_Power;
uniform bool Linear_Attenuation;
uniform float Specular_Index;

layout(location=0) in vec4 a_position;
layout(location=1) in vec3 a_normals;
layout(location=2) in vec4 a_Colors;

out vec4 colorvertex;
out vec3 Position_worldspace;
out vec3 Normal_cameraspace;
out vec3 EyeDirection_cameraspace;
out vec3 LightDirection_cameraspace;
out vec3 ambientcolor;
out vec3 lightcolor;
out vec3 specularcolor;
out float lightpower;
out float linearattenuation;
out float specularindex;

//varying vec2 v_texcoord;

void main()
{

    // Output position of the vertex, in clip space : MVP * position
    gl_Position = mvp_matrix * a_position;

//    Old version (sealed): Position of the vertex, in worldspace : M * position
//    Position_worldspace = (m_matrix * vec4(a_position.xyz,1)).xyz;
//    New version (working): Modified to get the object always illuminated from the observer's eye:
    vec4 v1 = mvp_matrix * a_position;
    Position_worldspace = v1.xyz;

//    Vector that goes from the vertex to the camera, in camera space.
//    In camera space, the camera is at the origin (0,0,0).
    v1 = mv_matrix * a_position;
    vec3 vertexPosition_cameraspace = v1.xyz;
    EyeDirection_cameraspace = vec3(0., 0., 0.) - vertexPosition_cameraspace;

//    Old version (sealed): Vector that goes from the vertex to the light, in camera space.
//    vec3 LightPosition_cameraspace = ( mv_matrix * vec4(LightPosition_worldspace,1)).xyz;
//    New version (working): Modified to get the object always illuminated from the observer's eye:
    vec3 LightPosition_cameraspace = LightPosition_worldspace;
    LightDirection_cameraspace = LightPosition_cameraspace + EyeDirection_cameraspace;

//    Normal of the the vertex, in camera space
//    Only correct if ModelMatrix does not scale the model ! Use its inverse transpose if not.
    v1 = mv_matrix * vec4(a_normals,0);
    Normal_cameraspace = v1.xyz;

    colorvertex = a_Colors;

//  Pass color of ambientlight to fragment shader
    ambientcolor = Ambient_Color;
    lightcolor = Light_Color;
    specularcolor = Specular_Color;
    specularindex = Specular_Index;

    lightpower = Light_Power;
    if (Linear_Attenuation)
        linearattenuation = 1.0;
    else
        linearattenuation = -1.0;
//  Pass texture coordinate to fragment shader
//  Value will be automatically interpolated to fragments inside polygon faces
//    v_texcoord = a_texcoord;
}
