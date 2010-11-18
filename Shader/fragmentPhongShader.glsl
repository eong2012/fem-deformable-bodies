varying vec4 v_color;
varying vec2 v_texCoord;

varying vec3 normal;

void main (void)
{
	vec3 n,lightDir;
	vec4 diffuse, ambient;
	float nDotLight;

	n = normalize(normal);
    n=n;
    lightDir = normalize(vec3(-1.0,1.0,1.0));
	nDotLight = max(dot(n,lightDir),0.0);


	ambient = vec4(0.2,0.2,0.2,1.0);
	diffuse = vec4(0.8,0.8,0.8,1.0)*nDotLight;
	diffuse = clamp(diffuse,0.0,1.0);

	gl_FragColor = ambient + diffuse;

}

