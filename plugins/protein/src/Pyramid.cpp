#include "Pyramid.h"
// #include "ShaderTools/VertexArrayObjects/Quad.h"
using namespace megamol::core::utility::log;


static bool firstFrame = true;

Pyramid::Pyramid()
{
    //intentionally empty
}
bool Pyramid::create(int width, int height, megamol::core::CoreInstance* ci)
{
    vislib::graphics::gl::ShaderSource vertSrc;
    vislib::graphics::gl::ShaderSource fragSrc;

    // Common Full Screen Vertex Shader
    if (!ci->ShaderSourceFactory().MakeShaderSource("protein::contour::vertex", vertSrc)) {
        Log::DefaultLog.WriteMsg(Log::LEVEL_ERROR,
        "%s: Unable to load vertex shader source for pyramid");
        return false;
    }

    // Create PULL Shader
    if (!ci->ShaderSourceFactory().MakeShaderSource("pullpush::pullNormal", fragSrc)) {
        Log::DefaultLog.WriteMsg(Log::LEVEL_ERROR,
            "%s: Unable to load fragment shader source for pull pyramid");
        return false;
    }
    try {
        if (!this->pullShaderProgram.Create(vertSrc.Code(), vertSrc.Count(), fragSrc.Code(), fragSrc.Count())) {
            throw vislib::Exception("Generic creation failure", __FILE__, __LINE__);
        }
    } catch (vislib::Exception e) {
        Log::DefaultLog.WriteMsg(Log::LEVEL_ERROR, "%s: Unable to create pull pyramid shader: %s\n", e.GetMsgA());
        return false;
    }

    // Create PUSH Shader
    if (!ci->ShaderSourceFactory().MakeShaderSource("pullpush::pushNormal", fragSrc)) {
        Log::DefaultLog.WriteMsg(Log::LEVEL_ERROR,
            "%s: Unable to load fragment shader source for push pyramid shader");
        return false;
    }
    try {
        if (!this->pushShaderProgram.Create(vertSrc.Code(), vertSrc.Count(), fragSrc.Code(), fragSrc.Count())) {
            throw vislib::Exception("Generic creation failure", __FILE__, __LINE__);
        }
    } catch (vislib::Exception e) {
        Log::DefaultLog.WriteMsg(Log::LEVEL_ERROR, "%s: Unable to create push pyramid shader: %s\n", e.GetMsgA());
        return false;
    }


    // Create VAO for  screen filling quad
    float quadVertices[] = { // vertex attributes for a quad that fills the entire screen in Normalized Device Coordinates.
        // positions   // texCoords
        -1.0f,  1.0f,  0.0f, 1.0f,
        -1.0f, -1.0f,  0.0f, 0.0f,
        1.0f, -1.0f,  1.0f, 0.0f,

        -1.0f,  1.0f,  0.0f, 1.0f,
        1.0f, -1.0f,  1.0f, 0.0f,
        1.0f,  1.0f,  1.0f, 1.0f
    };
    glGenVertexArrays(1, &VAO);
    glGenBuffers(1, &VBO);
    glBindVertexArray(VAO);
    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(quadVertices), &quadVertices, GL_STATIC_DRAW);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 4 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 4 * sizeof(float), (void*)(2 * sizeof(float)));
    glBindVertexArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, 0);

    // ShaderProgram* temp_SP = NULL;

    // if(pullShaderProgram)
    //     temp_SP = pullShaderProgram;
    // else if(pushShaderProgram)
    //     temp_SP = pushShaderProgram;

    // if (!temp_SP)
    // {
    //     std::cerr << "No shader paths given!" << std::endl;
    //     return;
    // }

    int numTextures = 1;

    std::vector<GLuint> textures(numTextures);
    std::vector<GLuint> drawBuffers(numTextures);

    glGenTextures(numTextures, &textures[0]);

    // Generate texture shader-output.
    // A mipmap is also generated 
    GLuint handle = textures[0];
    glBindTexture(GL_TEXTURE_2D, handle);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F, width, height, 0, GL_RGBA, GL_FLOAT, 0);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_GENERATE_MIPMAP, GL_TRUE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_BORDER);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_BORDER);
    float color[] = { 0.0f, 0.0f, 0.0f, 0.0f };
    glTexParameterfv(GL_TEXTURE_2D, GL_TEXTURE_BORDER_COLOR, color);

    glGenerateMipmap(GL_TEXTURE_2D);

    drawBuffers[0] = GL_COLOR_ATTACHMENT0 + 0;

    glUniform1i(this->pullShaderProgram.ParameterLocation("pyramid_fragNormal"), handle);
    glUniform1i(this->pushShaderProgram.ParameterLocation("pyramid_fragNormal"), handle);
    textureMap["fragNormal"] = handle;

    int mipmapNumber = (int)glm::log2(glm::max<float>(width, height)) + 1;

    // Generate FBOs for each mipmap-level
    fboHandles.resize(mipmapNumber);
    glGenFramebuffers(mipmapNumber, &fboHandles[0]);

    // Bind the mipmap-level i of each texture j to the FBO handle with index i.
    for (int i = 0; i < mipmapNumber; i++) {
        glBindFramebuffer(GL_FRAMEBUFFER, fboHandles[i]);
        for (int j = 0; j < numTextures; j++) {
            glFramebufferTexture2D(GL_FRAMEBUFFER, drawBuffers[j], GL_TEXTURE_2D, textures[j], i);
        }
        glNamedFramebufferDrawBuffers(fboHandles[i], numTextures, &drawBuffers[0]);
    }

    GLenum status;
    status = glCheckFramebufferStatusEXT(GL_FRAMEBUFFER);
    switch(status)
    {
    case GL_FRAMEBUFFER_COMPLETE:
        std::cout << "Pyramid: FBO complete" <<std::endl;
        break;
    default:
        std::cout << "Pyramid: FBO incomplete" <<std::endl;
    }

    glBindFramebuffer(GL_FRAMEBUFFER, 0);
    return true;
}

Pyramid::~Pyramid() {
}

Pyramid* Pyramid::pull() {
    std::cout << "Pull starting" << std::endl;
    pullShaderProgram.Enable();
    for (int level = 0; level < getMipmapNumber(); level++) {
        glBindFramebuffer(GL_FRAMEBUFFER, fboHandles[level]);

        glUniform1i(this->pullShaderProgram.ParameterLocation("level"), level);
        glUniform1d(this->pullShaderProgram.ParameterLocation("lf"), 1.0 / glm::pow(2,getMipmapNumber() - level - 1));
        glBindVertexArray(VAO);
        glDrawArrays(GL_TRIANGLES, 0, 6);
    }
    glBindFramebuffer(GL_FRAMEBUFFER, 0);

    return this;
}

Pyramid *Pyramid::pull_until(int target_level)
{
    std::cout << "Pull until starting" << std::endl;
    pullShaderProgram.Enable();
    for (int level = 0; level <= target_level; level++) {
        glBindFramebuffer(GL_FRAMEBUFFER, fboHandles[level]);

        glUniform1i(this->pullShaderProgram.ParameterLocation("level"), level);
        glUniform1d(this->pullShaderProgram.ParameterLocation("lf"), 1.0 / glm::pow(2,getMipmapNumber() - level - 1));
        glBindVertexArray(VAO);
        glDrawArrays(GL_TRIANGLES, 0, 6);
    }
    glBindFramebuffer(GL_FRAMEBUFFER, 0);

    return this;
}

Pyramid *Pyramid::push_from(int start_level)
{
    std::cout << "Push from starting" << std::endl;
    for(int i = start_level - 1; i >= 0; i--)
        this->push(i);
	return this;
}

Pyramid* Pyramid::push() {
    std::cout << "Push starting" << std::endl;
    pushShaderProgram.Enable();
    glUniform1i(this->pushShaderProgram.ParameterLocation("level_max"), getMipmapNumber());
    for (int level = getMipmapNumber() - 2; level >= 0; level--) {

        glBindFramebuffer(GL_FRAMEBUFFER, fboHandles[level]);
        glUniform1i(this->pushShaderProgram.ParameterLocation("level"), level);
        glUniform1d(this->pushShaderProgram.ParameterLocation("lf"), 1.0 / glm::pow(2,getMipmapNumber() - level - 1));
        glBindVertexArray(VAO);
        glDrawArrays(GL_TRIANGLES, 0, 6);
    }
    glBindFramebuffer(GL_FRAMEBUFFER, 0);
    return this;
}

Pyramid *Pyramid::push(int level)
{
    std::cout << "Push int level" << std::endl;
    if(level > getMipmapNumber() - 1 || level < 0) {
        return this;
    }

    pushShaderProgram.Enable();
    glUniform1i(this->pushShaderProgram.ParameterLocation("level_max"), getMipmapNumber());
    glBindFramebuffer(GL_FRAMEBUFFER, fboHandles[level]);

    glUniform1i(this->pushShaderProgram.ParameterLocation("level"), level);
    glUniform1d(this->pushShaderProgram.ParameterLocation("lf"), 1.0 / glm::pow(2,getMipmapNumber() - level - 1));
    glBindVertexArray(VAO);
    glDrawArrays(GL_TRIANGLES, 0, 6);

    return this;
}

Pyramid* Pyramid::run() {
    std::cout << "run" << std::endl;
    pull();
    push();

    return this;
}

Pyramid* Pyramid::texture(std::string name, GLuint textureHandle) {
    std::cout << "texture" << std::endl;
    if (pullShaderProgram != NULL) {
        glUniform1i(this->pullShaderProgram.ParameterLocation("name"), textureHandle);
    }
    if (pushShaderProgram != NULL) {
        glUniform1i(this->pushShaderProgram.ParameterLocation("name"), textureHandle);
    }
    return this;
}

// Pyramid *Pyramid::texture(std::string name, GLuint textureID, GLuint samplerHandle, GLuint target)
// {
//     if (pullShaderProgram != NULL) {
//         glUniform1i(this->pullShaderProgram.ParameterLocation("name"), textureHandle);
//         pullShaderProgram->texture(name, textureID, samplerHandle, target);
//     }
//     if (pushShaderProgram != NULL) {
//         pushShaderProgram->texture(name, textureID, samplerHandle, target);
//     }
//     return this;
// }

Pyramid* Pyramid::clear(float r, float g, float b, float a, int level) {
    std::cout << "clear" << std::endl;
    glClearColor(r, g, b, a);

    if (level == -1) {
        for (auto f : fboHandles) {
            glBindFramebuffer(GL_FRAMEBUFFER, f);
            glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        }
    } else {
        glBindFramebuffer(GL_FRAMEBUFFER, fboHandles[level]);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    }

    glClearColor(0, 0, 0, 0);
    return this;
}

Pyramid* Pyramid::clear(int level) {
    std::cout << "clear int level" << std::endl;
    clear(0, 0, 0, 0, level);
    return this;
}

GLuint Pyramid::get(std::string name) {
    std::cout << "get" << std::endl;
    return textureMap[name];
}

int Pyramid::getMipmapNumber() {
    std::cout << "MipmapNumber" << std::endl;
    return fboHandles.size();
}