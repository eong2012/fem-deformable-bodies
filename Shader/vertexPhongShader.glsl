
varying vec2 v_texCoord;

varying vec3 normal;


void main( void )
{

    normal = normalize(gl_NormalMatrix * gl_Normal);
  

    gl_TexCoord[0] = gl_MultiTexCoord0*gl_TextureMatrix[0];
    gl_Position = gl_ProjectionMatrix * gl_ModelViewMatrix * (gl_Vertex);
    gl_Position = gl_Position;
   // v_texCoord = gl_MultiTexCoord0.xy;




}
