varying vec4 v_color;
varying vec2 v_texCoord;

varying vec3 normal;
varying vec3 vertex;

void main( void )
{

    normal = normalize(gl_NormalMatrix * -gl_Normal);
    vertex = normalize(gl_ModelViewMatrix * gl_Vertex);

    gl_TexCoord[0] = gl_MultiTexCoord0*gl_TextureMatrix[0];
    gl_Position = gl_ModelViewMatrix * (gl_Vertex);
    gl_Position = gl_Position;
   // v_texCoord = gl_MultiTexCoord0.xy;




}
