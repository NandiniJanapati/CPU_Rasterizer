#include "Triangle.h"
#include <GL/glew.h>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtx/normal.hpp>
#include <algorithm>
#define WINDOW_WIDTH 1024
#define WINDOW_HEIGHT 1024

// A function clamping the input values to the lower and higher bounds
#define CLAMP(in, low, high) ((in) < (low) ? (low) : ((in) > (high) ? (high) : in))

bool isInteger(float N)
{

	// Convert float value
	// of N to integer
	int X = N;

	float temp2 = N - X;

	// If N is not equivalent
	// to any integer
	if (temp2 > 0) {
		return false;
	}
	return true;
}

Triangle::Triangle()
{
	v[0] = glm::vec3(0.0f, 0.0f, 0.0f);
	v[1] = glm::vec3(0.0f, 0.0f, 0.0f);
	v[2] = glm::vec3(0.0f, 0.0f, 0.0f);

	c[0] = glm::vec3(0.0f, 0.0f, 0.0f);
	c[1] = glm::vec3(0.0f, 0.0f, 0.0f);
	c[2] = glm::vec3(0.0f, 0.0f, 0.0f);

	t[0] = glm::vec2(0.0f, 0.0f);
	t[1] = glm::vec2(0.0f, 0.0f);
	t[2] = glm::vec2(0.0f, 0.0f);
}

Triangle::Triangle(glm::vec3& v0, glm::vec3& v1, glm::vec3& v2)
{
	v[0] = v0;
	v[1] = v1;
	v[2] = v2;

	c[0] = glm::vec3(1.0f, 1.0f, 1.0f);
	c[1] = glm::vec3(1.0f, 1.0f, 1.0f);
	c[2] = glm::vec3(1.0f, 1.0f, 1.0f);

	t[0] = glm::vec2(0.0f, 0.0f);
	t[1] = glm::vec2(0.0f, 0.0f);
	t[2] = glm::vec2(0.0f, 0.0f);

};

Triangle::Triangle(glm::vec3& v0, glm::vec3& v1, glm::vec3& v2, glm::vec2& t0, glm::vec2& t1, glm::vec2& t2)
{
	v[0] = v0;
	v[1] = v1;
	v[2] = v2;

	t[0] = t0;
	t[1] = t1;
	t[2] = t2;

	c[0] = glm::vec3(1.0f, 1.0f, 1.0f);
	c[1] = glm::vec3(1.0f, 1.0f, 1.0f);
	c[2] = glm::vec3(1.0f, 1.0f, 1.0f);

};


Triangle::Triangle(glm::vec3& v0, glm::vec3& v1, glm::vec3& v2, glm::vec3& color) {
	v[0] = v0;
	v[1] = v1;
	v[2] = v2;

	c[0] = color;
	c[1] = color;
	c[2] = color;

	t[0] = glm::vec2(0.0f, 0.0f);
	t[1] = glm::vec2(0.0f, 0.0f);
	t[2] = glm::vec2(0.0f, 0.0f);
}

Triangle::Triangle(glm::vec3& v0, glm::vec3& v1, glm::vec3& v2, glm::vec2& t0, glm::vec2& t1, glm::vec2& t2, glm::vec3& color) {
	v[0] = v0;
	v[1] = v1;
	v[2] = v2;

	t[0] = t0;
	t[1] = t1;
	t[2] = t2;

	c[0] = color;
	c[1] = color;
	c[2] = color;
}

Triangle::Triangle(glm::vec3& v0, glm::vec3& v1, glm::vec3& v2, glm::vec3& c0, glm::vec3& c1, glm::vec3& c2) {
	v[0] = v0;
	v[1] = v1;
	v[2] = v2;

	c[0] = c0;
	c[1] = c1;
	c[2] = c2;

	t[0] = glm::vec2(0.0f, 0.0f);
	t[1] = glm::vec2(0.0f, 0.0f);
	t[2] = glm::vec2(0.0f, 0.0f);
}

// Rendering the triangle using OpenGL
void Triangle::RenderOpenGL(glm::mat4 &modelViewMatrix, glm::mat4 &projectionMatrix, bool isTextured)
{

	glMatrixMode(GL_MODELVIEW);
	glLoadMatrixf(glm::value_ptr(modelViewMatrix));

	glMatrixMode(GL_PROJECTION);
	glLoadMatrixf(glm::value_ptr(projectionMatrix));
	
	// For textured object
	if (isTextured)
	{
		glEnable(GL_TEXTURE_2D);

		// Avoid modulating the texture by vertex color
		glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);

		glBegin(GL_TRIANGLES);

			glTexCoord2f(t[0].x, t[0].y);
			glVertex3f(v[0].x, v[0].y, v[0].z);

			glTexCoord2f(t[1].x, t[1].y);
			glVertex3f(v[1].x, v[1].y, v[1].z);

			glTexCoord2f(t[2].x, t[2].y);
			glVertex3f(v[2].x, v[2].y, v[2].z);

		glEnd();

		glDisable(GL_TEXTURE_2D);


	}
	// For object with only vertex color
	else
	{
		glBegin(GL_TRIANGLES);

			glColor3f(c[0].x, c[0].y, c[0].z);
			glVertex3f(v[0].x, v[0].y, v[0].z);

			glColor3f(c[1].x, c[1].y, c[1].z);
			glVertex3f(v[1].x, v[1].y, v[1].z);

			glColor3f(c[2].x, c[2].y, c[2].z);
			glVertex3f(v[2].x, v[2].y, v[2].z);
		
		glEnd();
	}

}

void wrapUV(float& u, float& v) {
	while (u >= 1024) {
		u = u - 1024;
	}
	while (u < 0.0) {
		u = 1024 + u;
	}
	while (v >= 1024) {
		v = v - 1024;
	}
	while (v < 0.0) {
		v = 1024 + v;
	}
	return;
}

void wrapUV(int& u, int& v) {
	while (u >= 1024) {
		u = u - 1204;
	}
	while (u < 0.0) {
		u = 1024 + u;
	}
	while (v >= 1024) {
		v = v - 1024;
	}
	while (v < 0.0) {
		v = 1024 + v;
	}
	return;
}

bool Triangle::inside(float x, float y, glm::vec3 vertexs[3], float& alpha, float& beta, float& gamma) {
	float xA = vertexs[0].x; float yA = vertexs[0].y; 
	float xB = vertexs[1].x; float yB = vertexs[1].y; 
	float xC = vertexs[2].x; float yC = vertexs[2].y;

	alpha = ( -1*(x - xB) * (yC - yB) + (y - yB) * (xC - xB) ) / ( -1*(xA - xB) * (yC - yB) + (yA - yB) * (xC - xB) );
	if (alpha < 0) return false;
	beta = (-1 * (x - xC) * (yA - yC) + (y - yC) * (xA - xC)) / (-1 * (xB - xC) * (yA - yC) + (yB - yC) * (xA - xC));
	if (beta < 0) return false;
	gamma = 1 - alpha - beta;
	if (gamma < 0) return false;


	/*
	//we will draw to screen, so save depth in depth array
	int y_int = y;
	int x_int = x;
	//alpha*vertex[0].color + beta*vertex[1].color + gamma*vertex[2].color
	//alpha*vertex[0].z + beta*vertex[1].z + gamma*vertex[2].z

	float zdepth = alpha * vertexs[0].z + beta * vertexs[1].z + gamma * vertexs[2].z;

	if ((*depth)[y_int][x_int] > zdepth) {
		(*color)[y_int][x_int][0] = CLAMP(alpha * c[0].r + beta * c[1].r + gamma * c[2].r, 0.0, 1.0);
		(*color)[y_int][x_int][1] = CLAMP(alpha * c[0].g + beta * c[1].g + gamma * c[2].g, 0.0, 1.0);
		(*color)[y_int][x_int][2] = CLAMP(alpha * c[0].b + beta * c[1].b + gamma * c[2].b, 0.0, 1.0);
		
		(*depth)[y_int][x_int] = zdepth;
	}
	*/
	
	
	return true;

}

void Triangle::ColorPixel(float alpha, float beta, float gamma, float x, float y, glm::vec3 vertexs[3], float(*color)[WINDOW_HEIGHT][WINDOW_WIDTH][3], float(*depth)[WINDOW_HEIGHT][WINDOW_WIDTH]) {
	//we will draw to screen, so save depth in depth array
	int y_int = y;
	int x_int = x;
	//alpha*vertex[0].color + beta*vertex[1].color + gamma*vertex[2].color
	//alpha*vertex[0].z + beta*vertex[1].z + gamma*vertex[2].z

	float zdepth = alpha * vertexs[0].z + beta * vertexs[1].z + gamma * vertexs[2].z;

	if ((*depth)[y_int][x_int] > zdepth) {
		(*color)[y_int][x_int][0] = CLAMP(alpha * c[0].r + beta * c[1].r + gamma * c[2].r, 0.0, 1.0);
		(*color)[y_int][x_int][1] = CLAMP(alpha * c[0].g + beta * c[1].g + gamma * c[2].g, 0.0, 1.0);
		(*color)[y_int][x_int][2] = CLAMP(alpha * c[0].b + beta * c[1].b + gamma * c[2].b, 0.0, 1.0);

		(*depth)[y_int][x_int] = zdepth;
	}
}

void BaryCentricCoords(float& alpha, float& beta, float& gamma, float x, float y, glm::vec3 vertexs[3]) {
	float xA = vertexs[0].x; float yA = vertexs[0].y;
	float xB = vertexs[1].x; float yB = vertexs[1].y;
	float xC = vertexs[2].x; float yC = vertexs[2].y;
	alpha = (-1 * (x - xB) * (yC - yB) + (y - yB) * (xC - xB)) / (-1 * (xA - xB) * (yC - yB) + (yA - yB) * (xC - xB));
	beta = (-1 * (x - xC) * (yA - yC) + (y - yC) * (xA - xC)) / (-1 * (xB - xC) * (yA - yC) + (yB - yC) * (xA - xC));
	gamma = 1 - alpha - beta;
	return;
}

void Triangle::FindTexCoord(float alpha, float beta, float gamma, float x, float y, glm::vec3 vertexs[3], glm::vec4 worldSpaceCoords[3], float& u_final, float& v_final, int textureWidth, int textureHeight) {
	//alpha, beta, and gamma are calculated from screen space coordinates' vertices
	
	
	float Zinv0 = 1 / worldSpaceCoords[0][2];
	float Zinv1 = 1 / worldSpaceCoords[1][2];
	float Zinv2 = 1 / worldSpaceCoords[2][2];

	glm::vec2 Qsca0 = t[0] / worldSpaceCoords[0][2];
	glm::vec2 Qsca1 = t[1] / worldSpaceCoords[1][2];
	glm::vec2 Qsca2 = t[2] / worldSpaceCoords[2][2];

	float Zinv_for_pixel = alpha * Zinv0 + beta * Zinv1 + gamma * Zinv2;
	glm::vec2 Qsca_for_pixel = alpha * Qsca0 + beta * Qsca1 + gamma * Qsca2;
	

	float u = Qsca_for_pixel[0] / Zinv_for_pixel;
	float v = Qsca_for_pixel[1] / Zinv_for_pixel;

	if (u >= 1.0) {
		u = u - 1;
	}
	else if (u < 0.0) {
		u = 1 + u;
	}
	if (v >= 1.0) {
		v = v - 1;
	}
	else if (v < 0.0) {
		v = 1 + v;
	}

	//scaling
	u = u * 1024;
	v = v * 1024;

	u_final = u;
	v_final = v;

	return;
}

void Triangle::FindUnwrappedTexCoord(float alpha, float beta, float gamma, float x, float y, glm::vec3 vertexs[3], glm::vec4 worldSpaceCoords[3], float& u_final, float& v_final) {
	float Zinv0 = 1 / worldSpaceCoords[0][2];
	float Zinv1 = 1 / worldSpaceCoords[1][2];
	float Zinv2 = 1 / worldSpaceCoords[2][2];

	glm::vec2 Qsca0 = t[0] / worldSpaceCoords[0][2];
	glm::vec2 Qsca1 = t[1] / worldSpaceCoords[1][2];
	glm::vec2 Qsca2 = t[2] / worldSpaceCoords[2][2];

	float Zinv_for_pixel = alpha * Zinv0 + beta * Zinv1 + gamma * Zinv2;
	glm::vec2 Qsca_for_pixel = alpha * Qsca0 + beta * Qsca1 + gamma * Qsca2;


	float u = Qsca_for_pixel[0] / Zinv_for_pixel;
	float v = Qsca_for_pixel[1] / Zinv_for_pixel;

	u = u * 1024;
	v = v * 1024;

	u_final = u;
	v_final = v;
	return;
}

void ColorNearestNeighbor(float x, float y, float u, float v, float alpha, float beta, float gamma, int texWidth, glm::vec3 vertexs[3], float(*color)[WINDOW_HEIGHT][WINDOW_WIDTH][3], float(*depth)[WINDOW_HEIGHT][WINDOW_WIDTH], std::vector<float*>& texture) {
	int y_int = y;
	int x_int = x;
	

	float zdepth = alpha * vertexs[0].z + beta * vertexs[1].z + gamma * vertexs[2].z;
	
	if ((*depth)[y_int][x_int] > zdepth) {

		int v_int = floor(v);
		int u_int = floor(u);
		float r = texture[0][v_int * texWidth * 3 + u_int * 3 + 0];
		float g = texture[0][v_int * texWidth * 3 + u_int * 3 + 1];
		float b = texture[0][v_int * texWidth * 3 + u_int * 3 + 2];
		(*color)[y_int][x_int][0] = r;
		(*color)[y_int][x_int][1] = g;
		(*color)[y_int][x_int][2] = b;
		(*depth)[y_int][x_int] = zdepth;
	}

}

float linterpolate(float a, float b, float t) {
	float tmp = a * (1 - t) + b*t;
	return tmp;
}

void bilinterpolation(float u, float v, float &r, float& g, float& b, int D, std::vector<float*>& texture) {

	float ushift, vshift;
	ushift = u - 0.5;
	vshift = v - 0.5;

	float flooru, floorv, ceilu, ceilv;

	ceilu = ceil(ushift);
	ceilv = ceil(vshift);
	flooru = floor(ushift);
	floorv = floor(vshift);

	float S = ushift - flooru;
	float T = vshift - floorv;

	wrapUV(ushift, vshift);
	wrapUV(ceilu, ceilv);
	wrapUV(flooru, floorv);

	int P1[2] = { flooru, ceilv };
	float P1rgb[3] = { texture[D][P1[1] * 1024 * 3 + P1[0] * 3 + 0], texture[D][P1[1] * 1024 * 3 + P1[0] * 3 + 1], texture[D][P1[1] * 1024 * 3 + P1[0] * 3 + 2] };

	int P2[2] = { ceilu, ceilv };
	float P2rgb[3] = { texture[D][P2[1] * 1024 * 3 + P2[0] * 3 + 0], texture[D][P2[1] * 1024 * 3 + P2[0] * 3 + 1], texture[D][P2[1] * 1024 * 3 + P2[0] * 3 + 2] };

	int P3[2] = { flooru, floorv };
	float P3rgb[3] = { texture[D][P3[1] * 1024 * 3 + P3[0] * 3 + 0], texture[D][P3[1] * 1024 * 3 + P3[0] * 3 + 1], texture[D][P3[1] * 1024 * 3 + P3[0] * 3 + 2] };

	int P4[2] = { ceilu, floorv };
	float P4rgb[3] = { texture[D][P4[1] * 1024 * 3 + P4[0] * 3 + 0], texture[D][P4[1] * 1024 * 3 + P4[0] * 3 + 1], texture[D][P4[1] * 1024 * 3 + P4[0] * 3 + 2] };

	float Argb[3];
	float Brgb[3];
	for (int i = 0; i < 3; i++) {
		Argb[i] = linterpolate(P1rgb[i], P2rgb[i], S);
		Brgb[i] = linterpolate(P3rgb[i], P4rgb[i], S);
	}

	r = linterpolate(Brgb[0], Argb[0], T);
	g = linterpolate(Brgb[1], Argb[1], T);
	b = linterpolate(Brgb[2], Argb[2], T);

	return;
}

void ColorBilinear(float x, float y, float u, float v, float alpha, float beta, float gamma, glm::vec3 vertexs[3], float(*color)[WINDOW_HEIGHT][WINDOW_WIDTH][3], float(*depth)[WINDOW_HEIGHT][WINDOW_WIDTH], std::vector<float*>& texture) {
	int y_int = y;
	int x_int = x;

	float zdepth = alpha * vertexs[0].z + beta * vertexs[1].z + gamma * vertexs[2].z;

	if ((*depth)[y_int][x_int] > zdepth) {

		float r, g, b;
		bilinterpolation(u, v, r, g, b, 0, texture);

		/*float ushift, vshift;
		ushift = u - 0.5;
		vshift = v - 0.5;

		float flooru, floorv, ceilu, ceilv;

		ceilu = ceil(ushift);
		ceilv = ceil(vshift);
		flooru = floor(ushift);
		floorv = floor(vshift);

		float S = ushift - flooru;
		float T = vshift - floorv;

		wrapUV(ushift, vshift);
		wrapUV(ceilu, ceilv);
		wrapUV(flooru, floorv);


		int P1[2] = { flooru, ceilv };
		float P1rgb[3] = { texture[0][P1[1] * 1024 * 3 + P1[0] * 3 + 0], texture[0][P1[1] * 1024 * 3 + P1[0] * 3 + 1], texture[0][P1[1] * 1024 * 3 + P1[0] * 3 + 2] };

		int P2[2] = { ceilu, ceilv };
		float P2rgb[3] = { texture[0][P2[1] * 1024 * 3 + P2[0] * 3 + 0], texture[0][P2[1] * 1024 * 3 + P2[0] * 3 + 1], texture[0][P2[1] * 1024 * 3 + P2[0] * 3 + 2] };

		int P3[2] = { flooru, floorv };
		float P3rgb[3] = { texture[0][P3[1] * 1024 * 3 + P3[0] * 3 + 0], texture[0][P3[1] * 1024 * 3 + P3[0] * 3 + 1], texture[0][P3[1] * 1024 * 3 + P3[0] * 3 + 2] };

		int P4[2] = { ceilu, floorv };
		float P4rgb[3] = { texture[0][P4[1] * 1024 * 3 + P4[0] * 3 + 0], texture[0][P4[1] * 1024 * 3 + P4[0] * 3 + 1], texture[0][P4[1] * 1024 * 3 + P4[0] * 3 + 2] };
		
		float Argb[3];
		float Brgb[3];
		for (int i = 0; i < 3; i++) {
			Argb[i] = linterpolate(P1rgb[i], P2rgb[i], S);
			Brgb[i] = linterpolate(P3rgb[i], P4rgb[i], S);
		}

		float r = linterpolate(Brgb[0], Argb[0], T);
		float g = linterpolate(Brgb[1], Argb[1], T);
		float b = linterpolate(Brgb[2], Argb[2], T);*/
		
		//bilinterpolation(u, v, r, g, b, 0, texture);

		(*color)[y_int][x_int][0] = r;
		(*color)[y_int][x_int][1] = g;
		(*color)[y_int][x_int][2] = b;
		(*depth)[y_int][x_int] = zdepth;
	}

}

void Triangle::ColorMipMap(float x, float y, float u, float v, float alpha, float beta, float gamma, glm::vec4 worldSpaceCoords[3], glm::vec3 vertexs[3], float(*color)[WINDOW_HEIGHT][WINDOW_WIDTH][3], float(*depth)[WINDOW_HEIGHT][WINDOW_WIDTH], std::vector<float*>& texture) {
	int x_int = x;
	int y_int = y;

	float zdepth = alpha * vertexs[0].z + beta * vertexs[1].z + gamma * vertexs[2].z;

	if ((*depth)[y_int][x_int] > zdepth) {

		float rightPoint[2] = { x + 1, y };
		float upPoint[2] = { x, y + 1 };

		float Texel[2]; 
		FindUnwrappedTexCoord(alpha, beta, gamma, x, y, vertexs, worldSpaceCoords, Texel[0], Texel[1]); //texel is scaled but unwrapped

		float TexelRight[2]; //u,v point
		//calculate texel coords for right point
		float tmpalpha, tmpbeta, tmpgamma;
		BaryCentricCoords(tmpalpha, tmpbeta, tmpgamma, rightPoint[0], rightPoint[1], vertexs);
		FindUnwrappedTexCoord(tmpalpha, tmpbeta, tmpgamma, rightPoint[0], rightPoint[1], vertexs, worldSpaceCoords, TexelRight[0], TexelRight[1]); //scaled but unwrapped


		float TexelUp[2]; //u,v point
		//calculate texel coords for up point
		BaryCentricCoords(tmpalpha, tmpbeta, tmpgamma, upPoint[0], upPoint[1], vertexs);
		FindUnwrappedTexCoord(tmpalpha, tmpbeta, tmpgamma, upPoint[0], upPoint[1], vertexs, worldSpaceCoords, TexelUp[0], TexelUp[1]);


		float vectorA[2] = { TexelRight[0] - Texel[0], TexelRight[1] - Texel[1]}; //vector from initial texel to pixel to the right's texel; ( du/dx, dy/dx )
		float vectorB[2] = { TexelUp[0] - Texel[0], TexelUp[1] - Texel[1]}; //vector from initial texel to pixel above's texel; ( du/dy, dv/dy )

		float lengthA = sqrt(pow(vectorA[0], 2) + pow(vectorA[1], 2));
		float lengthB = sqrt(pow(vectorB[0], 2) + pow(vectorB[1], 2));

		float L = std::max(lengthA, lengthB);
	
		if (L < 1) {
			L = 1;
		}
		else if (L > 1024) {
			L = 1024;
		}
		float D = log2(L);

		int Dmin, Dmax;
		Dmin = floor(D);
		Dmax = ceil(D);

		float RDmin, GDmin, BDmin;
		float RDmax, GDmax, BDmax;

		bilinterpolation(Texel[0], Texel[1], RDmin, GDmin, BDmin, Dmin, texture);
		bilinterpolation(Texel[0], Texel[1], RDmax, GDmax, BDmax, Dmax, texture);

		float rFinal = linterpolate(RDmin, RDmax, D - Dmin);
		float gFinal = linterpolate(GDmin, GDmax, D - Dmin);
		float bFinal = linterpolate(BDmin, BDmax, D - Dmin);

		(*color)[y_int][x_int][0] = rFinal;
		(*color)[y_int][x_int][1] = gFinal;
		(*color)[y_int][x_int][2] = bFinal;
		(*depth)[y_int][x_int] = zdepth;

	}

}

// Render the triangle on CPU
void Triangle::RenderCPU(glm::mat4& modelViewMatrix, glm::mat4& projectionMatrix, float (*color)[WINDOW_HEIGHT][WINDOW_WIDTH][3], float (*depth)[WINDOW_HEIGHT][WINDOW_WIDTH])
{
	//values for determining the bounding box of triangle
	float minX = INT_MAX;
	float minY = INT_MAX;
	float maxX = INT_MIN;
	float maxY = INT_MIN;

	std::vector<glm::vec3> inputs; //triangle's v[0-2]
	inputs.push_back(v[0]); inputs.push_back(v[1]); inputs.push_back(v[2]);
	glm::vec3 vertexs[3]; //for storing viewport versions of triangle verts

	for (int i = 0; i < 3; i++) { //for each vertex in the triangle and the eye direction vector
		glm::vec4 vertex;
		vertex[0] = inputs[i].x; vertex[1] = inputs[i].y; vertex[2] = inputs[i].z; 
		
		vertex[3] = 1;
		//if (i < 3) vertex[3] = 1;
		//else vertex[3] = 0;
		vertex = modelViewMatrix * vertex;
		vertex = projectionMatrix * vertex;

		//divide by w
		vertex[0] = vertex[0] / vertex[3];
		vertex[1] = vertex[1] / vertex[3];
		vertex[2] = vertex[2] / vertex[3];
		vertex[3] = vertex[3] / vertex[3];
		

		float temp = static_cast<float>(1024 / 2);
		glm::mat4 viewport( temp, 0.0f, 0.0f, 0.0f, 0.0f, temp, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f, temp, temp, 0.0f, 1.0f ); //column major
		vertex = viewport * vertex;

		if (vertex[0] > maxX) { //finding the bounding box for the triangle
			maxX = vertex[0];
		}
		if (vertex[1] > maxY) {
			maxY = vertex[1];
		}
		if (vertex[0] < minX) {
			minX = vertex[0];
		}
		if (vertex[1] < minY) {
			minY = vertex[1];
		}
		vertexs[i] = vertex;
	}


	//check if pixel is inside triangle
	int ceilXmax = ceil(maxX) + 0.4; //ceil returns a float so ceil(maxX) when maxX = 2.3 could be 2.9999999999
	int ceilYmax = ceil(maxY) + 0.4; //add a small decimal to so that it is for sure above int we want so int will floor to correct value
	for (int i = floor(minX) + 0.4; i < ceilXmax; i++) { //test the pixels inside the bounding box.
		for (int j = floor(minY) + 0.4; j < ceilYmax; j++) { 
			float alpha, beta, gamma;
			bool isInside = inside(i+0.5, j+0.5, vertexs, alpha, beta, gamma);
			if (isInside) {
				ColorPixel(alpha, beta, gamma, i+0.5, j+0.5, vertexs, color, depth);
			}
		}
	}
	
}

void Triangle::RenderCPUforTexture(glm::mat4& modelViewMatrix, glm::mat4& projectionMatrix, bool isTextured, float(*color)[WINDOW_HEIGHT][WINDOW_WIDTH][3], float(*depth)[WINDOW_HEIGHT][WINDOW_WIDTH], std::vector<float*>& texture, int textureWidth, int textureHeight, int textureMode)
{
	//values for determining the bounding box of triangle
	float minX = INT_MAX;
	float minY = INT_MAX;
	float maxX = INT_MIN;
	float maxY = INT_MIN;

	std::vector<glm::vec3> inputs; //triangle's v[0-2]
	inputs.push_back(v[0]); inputs.push_back(v[1]); inputs.push_back(v[2]);
	glm::vec3 vertexs[3]; //for storing viewport versions of triangle verts

	glm::vec4 worldSpaceCoords[3]; //for storing worldspace versions of triangle verts (used for texture coord interpolation)

	for (int i = 0; i < 3; i++) { //for each vertex in the triangle and the eye direction vector
		glm::vec4 vertex;
		vertex[0] = inputs[i].x; vertex[1] = inputs[i].y; vertex[2] = inputs[i].z;

		vertex[3] = 1;
		//if (i < 3) vertex[3] = 1;
		//else vertex[3] = 0;
		vertex = modelViewMatrix * vertex;

		worldSpaceCoords[i][0] = vertex[0];
		worldSpaceCoords[i][1] = vertex[1];
		worldSpaceCoords[i][2] = vertex[2];
		worldSpaceCoords[i][3] = vertex[3];
		vertex = projectionMatrix * vertex;

		//divide by w
		vertex[0] = vertex[0] / vertex[3];
		vertex[1] = vertex[1] / vertex[3];
		vertex[2] = vertex[2] / vertex[3];
		vertex[3] = vertex[3] / vertex[3];


		float temp = static_cast<float>(1024 / 2);
		glm::mat4 viewport(temp, 0.0f, 0.0f, 0.0f, 0.0f, temp, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f, temp, temp, 0.0f, 1.0f); //column major
		vertex = viewport * vertex;

		if (vertex[0] > maxX) { //finding the bounding box for the triangle
			maxX = vertex[0];
		}
		if (vertex[1] > maxY) {
			maxY = vertex[1];
		}
		if (vertex[0] < minX) {
			minX = vertex[0];
		}
		if (vertex[1] < minY) {
			minY = vertex[1];
		}
		vertexs[i] = vertex;
	}

	//check if pixel is inside triangle
	int ceilXmax = ceil(maxX) + 0.4; //ceil returns a float so ceil(maxX) when maxX = 2.3 could be 2.9999999999
	int ceilYmax = ceil(maxY) + 0.4; //add a small decimal to it so that it is for sure above int we want so int will floor to correct value
	if (floor(minX) < 0 || floor(minY) < 0) {
		int g = 0;
	}
	for (int i = floor(minX) + 0.4; i < ceilXmax; i++) { //test the pixels inside the bounding box.
		for (int j = floor(minY) + 0.4; j < ceilYmax; j++) {
			float alpha, beta, gamma;
			bool isInside = inside(i + 0.5, j + 0.5, vertexs, alpha, beta, gamma);
			if (isInside) {
				float u, v;
				FindTexCoord(alpha, beta, gamma, i + 0.5, j + 0.5, vertexs, worldSpaceCoords, u, v, textureWidth, textureHeight);
				if (textureMode == 0) { //nearest neighbors
					ColorNearestNeighbor(i + 0.5, j + 0.5, u, v, alpha, beta, gamma, textureWidth, vertexs, color, depth, texture);
				}
				else if (textureMode == 1) { //bilinear interpolation
					FindUnwrappedTexCoord(alpha, beta, gamma, i + 0.5, j + 0.5, vertexs, worldSpaceCoords, u, v);
					ColorBilinear(i + 0.5, j + 0.5, u, v, alpha, beta, gamma, vertexs, color, depth, texture);
				}
				else {
					//FindUnwrappedTexCoord(alpha, beta, gamma, i + 0.5, j + 0.5, vertexs, worldSpaceCoords, u, v);
					ColorMipMap(i + 0.5, j + 0.5, u, v, alpha, beta, gamma, worldSpaceCoords, vertexs, color, depth, texture);
				}
			}
		}
	}
}


void Triangle::SetSolidTriangle(glm::vec3& color) {
	c[0] = color;
	c[1] = color;
	c[2] = color;
}

void Triangle::SetColorfulTriangle(glm::vec3& c0, glm::vec3& c1, glm::vec3& c2) {
	c[0] = c0;//c0 for vertex 0
	c[1] = c1;//c1 for vertex 1
	c[2] = c2;//c2 for vertex 2
}

glm::vec3 Triangle::getZcoords() {
	return glm::vec3(v[0][2], v[1][2], v[2][2]);
}