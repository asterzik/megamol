// #include "ShaderTools/VertexArrayObject.h"
// #include "ShaderTools/FrameBufferObject.h"
#include <string>
#include "mmcore/CoreInstance.h"
#include "vislib/graphics/gl/GLSLShader.h"

class Pyramid {
public:
    Pyramid();

    ~Pyramid();

    bool create(std::string name, int width, int height, megamol::core::CoreInstance* ci, std::string pullPath,
        std::string pushPath = "pullpush::pushNormal");
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


    GLuint* get(std::string name);
    int getMipmapNumber();

    vislib::graphics::gl::GLSLShader pushShaderProgram;
    vislib::graphics::gl::GLSLShader pullShaderProgram;
    // template <class T>
    // Pyramid* update(std::string name, T value) {
    // 	pushShaderProgram->update(name, value);
    // 	if (pullShaderProgram != NULL) {
    // 		pullShaderProgram->update(name, value);
    // 	}
    // 	return this;
    // }

private:
    std::map<std::string, GLuint> textureMap;
    std::string textureName;
    // megamol::core::CoreInstance* ci;


    // Shader programs to use
    // std::string vertexPullPath;
    // std::string fragmentPullPath;
    // std::string vertexPushPath;
    // std::string fragmentPushPath;
    // ShaderProgram* pushShaderProgram;
    // ShaderProgram* pullShaderProgram;

    // Vertex array object to render
    // VertexArrayObject* vertexArrayObject;
    GLuint VBO;
    GLuint VAO;

    // mip framebuffer objects to render to
    std::vector<GLuint> fboHandles;
};
