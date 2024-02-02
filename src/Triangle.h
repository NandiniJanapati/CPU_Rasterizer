#pragma once

#include <stdlib.h>
#include <math.h>
#include <vector>

#include <glm/glm.hpp>

#define WINDOW_WIDTH 1024
#define WINDOW_HEIGHT 1024


class Triangle {
	private:
		glm::vec3 v[3];		// Triangle vertices
		glm::vec3 c[3];		// Vertex color
		glm::vec2 t[3];		// Texture coordinates

	public:

		// Default constructor
		Triangle();

		// Constructor without texture coordinates
		Triangle(glm::vec3& v0, glm::vec3& v1, glm::vec3& v2);

		// Constructor with texture coordinates
		Triangle(glm::vec3& v0, glm::vec3& v1, glm::vec3& v2, glm::vec2& t0, glm::vec2& t1, glm::vec2& t2);

		//MY constructor: solid color triangle constructor (no texture coords)
		Triangle(glm::vec3& v0, glm::vec3& v1, glm::vec3& v2, glm::vec3& color);

		//MY constructor: each vertex different color triangle constructor (no tex coords)
		Triangle(glm::vec3& v0, glm::vec3& v1, glm::vec3& v2, glm::vec3& c0, glm::vec3& c1, glm::vec3& c2);

		//MY constructor: solid color triangle constructor (with texture coords)
		Triangle(glm::vec3& v0, glm::vec3& v1, glm::vec3& v2, glm::vec2& t0, glm::vec2& t1, glm::vec2& t2, glm::vec3& color);

		// Rendering the triangle using OpenGL
		void RenderOpenGL(glm::mat4 &modelViewMatrix, glm::mat4 &projectionMatrix, bool textureMode);

		// Rendering the triangle using CPU
		void RenderCPU(glm::mat4& modelViewMatrix, glm::mat4& projectionMatrix, float(*color)[WINDOW_HEIGHT][WINDOW_WIDTH][3], float(*depth)[WINDOW_HEIGHT][WINDOW_WIDTH]);

		void RenderCPUforTexture(glm::mat4& modelViewMatrix, glm::mat4& projectionMatrix, bool isTextured, float(*color)[WINDOW_HEIGHT][WINDOW_WIDTH][3], float(*depth)[WINDOW_HEIGHT][WINDOW_WIDTH], std::vector<float*>& texture, int textureWidth, int textureHeight, int textureMode);
		
		bool inside(float x, float y, glm::vec3 vertexs[3], float& alpha, float& beta, float& gamma);

		void FindTexCoord(float alpha, float beta, float gamma, float x, float y, glm::vec3 vertexs[3], glm::vec4 worldSpaceCoords[3], float& u_final, float& v_final, int textureWidth, int textureHeight);
		
		void FindUnwrappedTexCoord(float alpha, float beta, float gamma, float x, float y, glm::vec3 vertexs[3], glm::vec4 worldSpaceCoords[3], float& u_final, float& v_final);

		void ColorPixel(float alpha, float beta, float gamma, float x, float y, glm::vec3 vertexs[3], float(*color)[WINDOW_HEIGHT][WINDOW_WIDTH][3], float(*depth)[WINDOW_HEIGHT][WINDOW_WIDTH]);

		void ColorMipMap(float x, float y, float u, float v, float alpha, float beta, float gamma, glm::vec4 worldSpaceCoords[3], glm::vec3 vertexs[3], float(*color)[WINDOW_HEIGHT][WINDOW_WIDTH][3], float(*depth)[WINDOW_HEIGHT][WINDOW_WIDTH], std::vector<float*>& texture);
		
		void SetSolidTriangle(glm::vec3& color);

		void SetColorfulTriangle(glm::vec3& c0, glm::vec3& c1, glm::vec3& c2);

		glm::vec3 getZcoords();

};
