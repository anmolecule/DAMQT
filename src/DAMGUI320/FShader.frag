// #ifdef GL_ES
// // Set default precision to medium
// precision mediump int;
// precision mediump float;
// #endif

#version 330 core

//uniform sampler2D texture;

//varying vec2 v_texcoord;

in vec4 colorvertex;
in vec3 Position_worldspace;
in vec3 Normal_cameraspace;
in vec3 EyeDirection_cameraspace;
in vec3 LightDirection_cameraspace;
in vec3 ambientcolor;
in vec3 lightcolor;
in vec3 specularcolor;

in float lightpower;
in float linearattenuation;

in float specularindex;

out vec4 color;
//varying out vec4 color;

uniform mat4 MV;
uniform vec3 LightPosition_worldspace;

void main()
{
    // Light emission properties
    // You probably want to put them as uniforms
    vec3 LightColor = lightcolor;
    float LightPower = lightpower;

    // Material properties
    vec3 MaterialDiffuseColor = colorvertex.xyz;
    vec3 MaterialAmbientColor = ambientcolor * MaterialDiffuseColor;
    vec3 MaterialSpecularColor = specularcolor;

    // Distance to the light
    float distance = length( LightPosition_worldspace - Position_worldspace );

    // Type of light attenuation
    float decay;
    if (linearattenuation > 0.)
        decay = 1. / distance;
    else
        decay = 1. / (distance*distance);

    // Normal of the computed fragment, in camera space
    vec3 n = normalize( Normal_cameraspace );
    // Direction of the light (from the fragment to the light)
    vec3 l = normalize( LightDirection_cameraspace );
    // Cosine of the angle between the normal and the light direction,
    // clamped above 0
    //  - light is at the vertical of the triangle -> 1
    //  - light is perpendicular to the triangle -> 0
    //  - light is behind the triangle -> 0
    float cosTheta = clamp( dot( n,l ), 0., 1. );

    // Eye vector (towards the camera)
    vec3 E = normalize(EyeDirection_cameraspace);
    // Direction in which the triangle reflects the light
    vec3 R = reflect(-l,n);
    // Cosine of the angle between the Eye vector and the Reflect vector,
    // clamped to 0
    //  - Looking into the reflection -> 1
    //  - Looking elsewhere -> < 1
    float cosAlpha = clamp( dot( E,R ), 0., 1. );

    float cospow = specularindex;
//    cospow = 5;
    vec3 rawcolor =
        // Ambient : simulates indirect lighting
        MaterialAmbientColor +
        // Diffuse : "color" of the object
        MaterialDiffuseColor * LightColor * LightPower * cosTheta * decay +
        // Specular : reflective highlight, like a mirror
        MaterialSpecularColor * LightColor * LightPower * pow(cosAlpha,cospow) * decay;

    // Set fragment color from texture
    float opacity = colorvertex.w;
//    gl_FragColor = vec4(rawcolor, opacity);
    color = vec4(rawcolor, opacity);
}
