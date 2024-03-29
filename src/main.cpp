#define GLEW_STATIC
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#define STB_IMAGE_RESIZE_IMPLEMENTATION
#include "stb_image_resize.h"
#define TINYOBJLOADER_IMPLEMENTATION
#include "tiny_obj_loader.h"
#include "Triangle.h"


#define WINDOW_WIDTH 1024
#define WINDOW_HEIGHT 1024

GLFWwindow *window;
bool lButtonPressed;
bool rButtonPressed;


float color[WINDOW_HEIGHT][WINDOW_WIDTH][3];
float depth[WINDOW_HEIGHT][WINDOW_WIDTH];


std::vector<Triangle> triangleVector;
float maxZcoord = INT_MIN;
float minZcoord = INT_MAX;
std::string loadOBJ = "sphere.obj";
std::vector<float*> texture;

bool isOpenGL = true; //TODO: change back to true
bool isTextured = false;
float eyeDistance = 5.0f;
int textureMode = 0;
int colorMode = 0;
float angle = 0;

std::string mainName = "Assignment3 - <Nandini Janapati>";

int texWidth, texHeight;

GLuint texID;

void ClearFrameBuffer()
{
	memset(&color[0][0][0], 0.0f, sizeof(float) * WINDOW_WIDTH * WINDOW_HEIGHT * 3);
}

void Display()
{	
	glm::mat4 projectionMatrix = glm::perspective(glm::radians(60.0f), float(WINDOW_WIDTH) / float(WINDOW_HEIGHT), 0.1f, 100.0f); //projection matrix
	glm::mat4 modelViewMatrix = glm::lookAt(eyeDistance * glm::vec3(0.0f, 0.0f, 1.0f), glm::vec3(0.0f, 0.0f, 0.0f), glm::vec3(0.0f, 1.0f, 0.0f)) * glm::rotate(glm::mat4(1.0f), angle, glm::vec3(0.0f, 1.0f, 0.0f)); //view matrix

	if (isOpenGL)
	{
		if (isTextured)
			glBindTexture(GL_TEXTURE_2D, texID);
		
		for (int i = 0; i < triangleVector.size(); i++)
			triangleVector[i].RenderOpenGL(modelViewMatrix, projectionMatrix, isTextured);
		
		if (isTextured)
			glBindTexture(GL_TEXTURE_2D, 0);
	}
	else
	{
		if (isTextured) {
			for (int i = 0; i < triangleVector.size(); i++) {
				triangleVector[i].RenderCPUforTexture(modelViewMatrix, projectionMatrix, isTextured, &color, &depth, texture, texWidth, texHeight, textureMode);

			}
		}
		else {
			for (int i = 0; i < triangleVector.size(); i++) {
				triangleVector[i].RenderCPU(modelViewMatrix, projectionMatrix, &color, &depth);

			}
		}
		
		glDrawPixels(WINDOW_WIDTH, WINDOW_HEIGHT, GL_RGB, GL_FLOAT, &color[0][0][0]);
		ClearFrameBuffer();
	}

	glFlush();
}


void ChangeColorMode() {
	float Lx = minZcoord;
	float Mx = maxZcoord;
	float m = 1 / (-1 * Lx + Mx);
	float b = -1 * Lx * m;

	for (int i = 0; i < triangleVector.size(); i++) {
		if (colorMode == 0) {
			float rand_red = static_cast<float>(rand()) / static_cast<float>(RAND_MAX);
			float rand_green = static_cast<float>(rand()) / RAND_MAX;
			float rand_blue = rand() / static_cast<float>(RAND_MAX);
			glm::vec3 rand_color(rand_red, rand_green, rand_blue);
			triangleVector[i].SetSolidTriangle(rand_color);
			
		}
		else if(colorMode == 1) {
			float rand_red = static_cast<float>(rand()) / RAND_MAX;
			float rand_green = static_cast<float>(rand()) / RAND_MAX;
			float rand_blue = static_cast<float>(rand()) / RAND_MAX;
			glm::vec3 rand_color0(rand_red, rand_green, rand_blue);
			rand_red = static_cast<float>(rand()) / RAND_MAX;
			rand_green = static_cast<float>(rand()) / RAND_MAX;
			rand_blue = static_cast<float>(rand()) / RAND_MAX;
			glm::vec3 rand_color1(rand_red, rand_green, rand_blue);
			rand_red = static_cast<float>(rand()) / RAND_MAX;
			rand_green = static_cast<float>(rand()) / RAND_MAX;
			rand_blue = static_cast<float>(rand()) / RAND_MAX;
			glm::vec3 rand_color2(rand_red, rand_green, rand_blue);
			triangleVector[i].SetColorfulTriangle(rand_color0, rand_color1, rand_color2);
		}
		else {
			glm::vec3 coords = triangleVector[i].getZcoords();
			coords[0] = m * coords[0] + b;
			coords[1] = m * coords[1] + b;
			coords[2] = m * coords[2] + b;
			glm::vec3 vert0(0, coords[0], 0);
			glm::vec3 vert1(0, coords[1], 0);
			glm::vec3 vert2(0, coords[2], 0);
			triangleVector[i].SetColorfulTriangle(vert0, vert1, vert2);

		}

	}

}

// Keyboard character callback function
void CharacterCallback(GLFWwindow* lWindow, unsigned int key)
{
	switch (key) 
	{
	case '0':
		colorMode = 0;
		ChangeColorMode();
		break;
	case '1':
		colorMode = 1;
		ChangeColorMode();
		break;
	case '2':
		colorMode = 2;
		ChangeColorMode();
		break;
	case 'w':
		eyeDistance *= (1 - 0.05);
		break;
	case 's':
		eyeDistance *= (1 + 0.05);
		break;
	case 'a':
		angle -= 0.01;
		break;
	case 'd':
		angle += 0.01;
		break;
	case ' ':
		isOpenGL = !isOpenGL;
		break;
	case 't':
	{
		if (!texture.empty())
			isTextured = !isTextured;
		break;
	}
		
	case 'n':
		glBindTexture(GL_TEXTURE_2D, texID);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		glBindTexture(GL_TEXTURE_2D, 0);
		textureMode = 0;
		break;
	case 'l':
		glBindTexture(GL_TEXTURE_2D, texID);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
		glBindTexture(GL_TEXTURE_2D, 0);
		textureMode = 1;
		break;
	case 'm':
		glBindTexture(GL_TEXTURE_2D, texID);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
		glBindTexture(GL_TEXTURE_2D, 0);
		textureMode = 2;
		break;
	case 'q':
		glfwSetWindowShouldClose(window, GLFW_TRUE);
		break;
	default:
		break;
	}


}

// Create a vector of triangles. Considers the texture coordinates if they are available.
void CreateTriangleVector(std::vector<glm::vec3> &vertices, std::vector<glm::vec2>& texCoords)
{
	for (int i = 0; i < vertices.size() / 3; i++)
	{
		Triangle myTriangle;

		if (texCoords.empty()) {
			float rand_red = static_cast<float>(rand()) / static_cast<float>(RAND_MAX);
			float rand_green = static_cast<float>(rand()) / RAND_MAX;
			float rand_blue = rand() / static_cast<float>(RAND_MAX);
			glm::vec3 rand_color(rand_red, rand_green, rand_blue);
			myTriangle = Triangle(vertices[i * 3 + 0], vertices[i * 3 + 1], vertices[i * 3 + 2], rand_color); 
			if (vertices[i * 3 + 2][2] > maxZcoord) {
				maxZcoord = vertices[i * 3 + 2][2];
			}
			if (vertices[i * 3 + 2][2] < minZcoord) {
				minZcoord = vertices[i * 3 + 2][2];
			}
		}	
		else {
			float rand_red = static_cast<float>(rand()) / static_cast<float>(RAND_MAX);
			float rand_green = static_cast<float>(rand()) / RAND_MAX;
			float rand_blue = rand() / static_cast<float>(RAND_MAX);
			glm::vec3 rand_color(rand_red, rand_green, rand_blue);
			myTriangle = Triangle(vertices[i * 3 + 0], vertices[i * 3 + 1], vertices[i * 3 + 2],
				texCoords[i * 3 + 0], texCoords[i * 3 + 1], texCoords[i * 3 + 2], rand_color);
			if (vertices[i * 3 + 2][2] > maxZcoord) {
				maxZcoord = vertices[i * 3 + 2][2];
			}
			if (vertices[i * 3 + 2][2] < minZcoord) {
				minZcoord = vertices[i * 3 + 2][2];
			}
		}
			

		triangleVector.push_back(myTriangle);
	}
}

// Load the geometry and texture coordinates if available
void LoadModel(char* name, std::vector<glm::vec3> &vertices, std::vector<glm::vec2>& texCoords)
{
	// Taken from Shinjiro Sueda with slight modification
	std::string meshName(name);
	tinyobj::attrib_t attrib;
	std::vector<tinyobj::shape_t> shapes;
	std::vector<tinyobj::material_t> materials;
	std::string errStr;
	bool rc = tinyobj::LoadObj(&attrib, &shapes, &materials, &errStr, meshName.c_str());
	if (!rc) {
		std::cerr << errStr << std::endl;
	}
	else {
		// Some OBJ files have different indices for vertex positions, normals,
		// and texture coordinates. For example, a cube corner vertex may have
		// three different normals. Here, we are going to duplicate all such
		// vertices.
		// Loop over shapes
		for (size_t s = 0; s < shapes.size(); s++) {
			// Loop over faces (polygons)
			size_t index_offset = 0;
			for (size_t f = 0; f < shapes[s].mesh.num_face_vertices.size(); f++) {
				size_t fv = shapes[s].mesh.num_face_vertices[f];
				// Loop over vertices in the face.
				for (size_t v = 0; v < fv; v++) {
					// access to vertex
					tinyobj::index_t idx = shapes[s].mesh.indices[index_offset + v];
					vertices.push_back(glm::vec3(attrib.vertices[3 * idx.vertex_index + 0],
												 attrib.vertices[3 * idx.vertex_index + 1],
												 attrib.vertices[3 * idx.vertex_index + 2]));
					if (!attrib.texcoords.empty()) {
						texCoords.push_back(glm::vec2(attrib.texcoords[2 * idx.texcoord_index + 0],
							attrib.texcoords[2 * idx.texcoord_index + 1]));
					}
				}
				index_offset += fv;
			}
		}
	}
}

// Load texture and create downsampled versions of it for mipmapping
void LoadTexture(char* name)
{
	std::string texName(name);
	int c;
	stbi_set_flip_vertically_on_load(true);
	stbi_hdr_to_ldr_gamma(1.0f);
	float* image = stbi_loadf(texName.c_str(), &texWidth, &texHeight, &c, 0);
	
	if (!image)
		std::cerr << texName << " not found" << std::endl;
	else if (c != 3)
		std::cerr << texName << " must have 3 channels (RGB)" << std::endl;
	else if ((texWidth % 2) != 0 || (texHeight % 2) != 0)
		std::cerr << " must be a power of 2" << std::endl;
	else
		texture.push_back(image);

	int length = std::min(texWidth, texHeight);
	int numLevels = log2(length);
	
	float** downImages = new float* [numLevels];
	for (int i = 0; i < numLevels; i++)
		downImages[i] = new float[texWidth * texHeight * c];

	for (int i = 0; i < numLevels; i++)
	{
		int curWidth = texWidth / pow(2, i + 1);
		int curHeight = texHeight / pow(2, i + 1);
		float* temp = new float[curWidth * curHeight * c];
		stbir_resize_float(image, texWidth, texHeight, 0, temp, curWidth, curHeight, 0, c);
		stbir_resize_float(temp, texWidth / pow(2, i + 1), texHeight / pow(2, i + 1), 0, downImages[i], texWidth, texHeight, 0, c);
		texture.push_back(downImages[i]);
		stbi_image_free(temp);
	}


	if (!texture.empty())
	{
		glGenTextures(1, &texID);
		glBindTexture(GL_TEXTURE_2D, texID);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, texWidth, texHeight, 0, GL_RGB, GL_FLOAT, texture[0]);
		glGenerateMipmap(GL_TEXTURE_2D);
		glBindTexture(GL_TEXTURE_2D, 0);
	}
		
}

std::string WindowTitle(std::string mainName)
{
	std::string hardwareName;
	if (isOpenGL)
		hardwareName = " - GPU";
	else
		hardwareName = " - CPU";

	std::string textureMethod;
	if (textureMode == 0)
		textureMethod = " - Nearest";
	else if (textureMode == 1)
		textureMethod = " - Bilinear";
	else if (textureMode == 2)
		textureMethod = " - Mipmap";

	std::string colorMethod;
	if (textureMode == 0)
		colorMethod = " - Mode 0";
	else if (textureMode == 1)
		colorMethod = " - Mode 1";
	else if (textureMode == 2)
		colorMethod = " - Mode 2";

	if (isTextured)
		return (mainName + hardwareName + std::string(" - Textured") + textureMethod);
	else
		return (mainName + hardwareName + std::string(" - Colored") + colorMethod);
}

void Init()
{
	glfwInit();
	glfwWindowHint(GLFW_RESIZABLE, GLFW_FALSE);
	window = glfwCreateWindow(WINDOW_WIDTH, WINDOW_HEIGHT, WindowTitle(mainName).c_str(), NULL, NULL);
	glfwMakeContextCurrent(window);
	glfwSetCharCallback(window, CharacterCallback);
	glewExperimental = GL_TRUE;
	glewInit();
	glViewport(0, 0, WINDOW_WIDTH, WINDOW_HEIGHT);
	glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
	glEnable(GL_DEPTH_TEST);

	ClearFrameBuffer();

	std::string temp = "../resources/" + loadOBJ;
	char* tmp2 = &temp[0];

	std::vector<glm::vec3> vertices;
	std::vector<glm::vec2> texCoords;
	LoadModel(tmp2, vertices, texCoords); //TODO: make .obj a variable
	
	if (!texCoords.empty())
	{
		LoadTexture("../resources/earth.jpg");
		if (texture.empty())
			isTextured = false;
	}
	else
		isTextured = false;
		
	CreateTriangleVector(vertices, texCoords);
	
}



int main()
{	
	std::cout << "please type the name of the file you want to load: ";
	std::cin >> loadOBJ;
	std::cout << std::endl;

	while( (loadOBJ != "bunny.obj") && (loadOBJ != "duck.obj") && (loadOBJ != "sphere.obj")) {
		std::cout << "please enter one of the following: bunny.obj, duck.obj, or sphere.obj: ";
		std::cin >> loadOBJ;
		std::cout << std::endl;
	}
	Init();
	
	while ( glfwWindowShouldClose(window) == 0) 
	{
		for (int i = 0; i < WINDOW_HEIGHT; i++) {
			for (int j = 0; j < WINDOW_WIDTH; j++) {
				depth[i][j] = INT_MAX;
			}
		}
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		Display();
		glfwSwapBuffers(window);
		glfwPollEvents();
		glfwSetWindowTitle(window, WindowTitle(mainName).c_str());
	}

	glfwTerminate();
	return 0;
}