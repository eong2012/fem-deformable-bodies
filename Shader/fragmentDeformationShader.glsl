uniform sampler2D positionTex;
void main (void)
{
   vec4 position = texture2D(positionTex, gl_TexCoord[0].st);
   position = position/2.0;
   gl_FragData[0] = vec4(position);

}

