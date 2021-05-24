#include "ShaderTools/VertexArrayObject.h"
#include "ShaderTools/FrameBufferObject.h"

class Pyramid
{
public:
    Pyramid(int width, int height, std::string pullFragmentShaderPath, std::string pushFragmentShaderPath = std::string());
	
	~Pyramid();

    Pyramid* push();
    Pyramid* push(int level);
    Pyramid* pull();
    Pyramid* pull_until(int target_level);
    Pyramid* push_from(int start_level);
	Pyramid* run();

	Pyramid* clear(float r, float g, float b, float a, int level = -1);
	Pyramid* clear(int level = -1);
    Pyramid* texture(std::string name, GLuint textureHandle);
	Pyramid* texture(std::string name, GLuint textureID, GLuint samplerID);
    Pyramid* texture(std::string name, GLuint textureID, GLuint samplerHandle, GLuint target);


	GLuint get(std::string name);
	int getMipmapNumber();

	template <class T>
	Pyramid* update(std::string name, T value) {
		pushShaderProgram->update(name, value);
		if (pullShaderProgram != NULL) {
			pullShaderProgram->update(name, value);
		}
		return this;
	}

private:
	std::map<std::string, GLuint> textureMap;

	// Shader programs to use
	ShaderProgram* pushShaderProgram;
	ShaderProgram* pullShaderProgram;

	// Vertex array object to render 
	VertexArrayObject* vertexArrayObject;

	// mip framebuffer objects to render to
	std::vector<GLuint> fboHandles;
};
