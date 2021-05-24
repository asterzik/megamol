#include "Pyramid.h"
#include "ShaderTools/VertexArrayObjects/Quad.h"


static bool firstFrame = true;

Pyramid::Pyramid(int width, int height, std::string pullFragmentShaderPath, std::string pushFragmentShaderPath)
{
    if(pullFragmentShaderPath == "")
        pullShaderProgram = NULL;
    else
        pullShaderProgram = new ShaderProgram("/ScreenSpaceParameterization/Pyramid/Hybrid/fullscreen.vert", pullFragmentShaderPath);

    if (pushFragmentShaderPath == "") {
        pushShaderProgram = NULL;
    } else {
        pushShaderProgram = new ShaderProgram("/ScreenSpaceParameterization/Pyramid/Hybrid/fullscreen.vert", pushFragmentShaderPath);
    }
    vertexArrayObject = new Quad();

    ShaderProgram* temp_SP = NULL;

    if(pullShaderProgram)
        temp_SP = pullShaderProgram;
    else if(pushShaderProgram)
        temp_SP = pushShaderProgram;

    if (!temp_SP)
    {
        std::cerr << "No shader paths given!" << std::endl;
        return;
    }

    int numTextures = temp_SP->outputMap.size();

    std::vector<GLuint> textures(numTextures);
    std::vector<GLuint> drawBuffers(numTextures);

    glGenTextures(numTextures, &textures[0]);

    // Generate textures for each shader-output.
    // A mipmap is also generated for each texture.
    for (auto e : temp_SP->outputMap) {
        GLuint handle = textures[e.second.location];
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

        drawBuffers[e.second.location] = GL_COLOR_ATTACHMENT0 + e.second.location;

        if (pullShaderProgram != NULL) {
            pullShaderProgram->texture("pyramid_" + e.first, handle);
        }
        if (pushShaderProgram != NULL) {
            pushShaderProgram->texture("pyramid_" + e.first, handle);
        }

        textureMap[e.first] = handle;
    }

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
}

Pyramid::~Pyramid() {
}

Pyramid* Pyramid::pull() {
    if (pullShaderProgram == NULL) {
        return this;
    }
    pullShaderProgram->use();
    for (int level = 0; level < getMipmapNumber(); level++) {
        glBindFramebuffer(GL_FRAMEBUFFER, fboHandles[level]);

        pullShaderProgram->update("level", level);
        pullShaderProgram->update("lf", 1.0 / glm::pow(2,getMipmapNumber() - level - 1));
        vertexArrayObject->draw();
    }
    glBindFramebuffer(GL_FRAMEBUFFER, 0);

    return this;
}

Pyramid *Pyramid::pull_until(int target_level)
{
    if (pullShaderProgram == NULL) {
        return this;
    }
    pullShaderProgram->use();
    for (int level = 0; level <= target_level; level++) {
        glBindFramebuffer(GL_FRAMEBUFFER, fboHandles[level]);

        pullShaderProgram->update("level", level);
        pullShaderProgram->update("lf", 1.0 / glm::pow(2,getMipmapNumber() - level - 1));
        vertexArrayObject->draw();
    }
    glBindFramebuffer(GL_FRAMEBUFFER, 0);

    return this;
}

Pyramid *Pyramid::push_from(int start_level)
{
    for(int i = start_level - 1; i >= 0; i--)
        this->push(i);
	return this;
}

Pyramid* Pyramid::push() {
    if (pushShaderProgram == NULL) {
        return this;
    }

    pushShaderProgram->use();
    pushShaderProgram->update("level_max", getMipmapNumber());
    for (int level = getMipmapNumber() - 2; level >= 0; level--) {
        glBindFramebuffer(GL_FRAMEBUFFER, fboHandles[level]);

        pushShaderProgram->update("level", level);
        pushShaderProgram->update("lf", 1.0 / glm::pow(2,getMipmapNumber() - level - 1));
        vertexArrayObject->draw();
    }
    glBindFramebuffer(GL_FRAMEBUFFER, 0);

    return this;
}

Pyramid *Pyramid::push(int level)
{
    if (pushShaderProgram == NULL
            || level > getMipmapNumber() - 1
            || level < 0) {
        return this;
    }

    pushShaderProgram->use();
    pushShaderProgram->update("level_max", getMipmapNumber());
    glBindFramebuffer(GL_FRAMEBUFFER, fboHandles[level]);

    float lf = 1.0 / glm::pow(2,getMipmapNumber() - level - 1);
    pushShaderProgram->update("level", level);
    pushShaderProgram->update("lf", lf);
    vertexArrayObject->draw();

    return this;
}

Pyramid* Pyramid::run() {
    pull();
    push();

    return this;
}

Pyramid* Pyramid::texture(std::string name, GLuint textureHandle) {
    if (pullShaderProgram != NULL) {
        pullShaderProgram->texture(name, textureHandle);
    }
    if (pushShaderProgram != NULL) {
        pushShaderProgram->texture(name, textureHandle);
    }
    return this;
}

Pyramid *Pyramid::texture(std::string name, GLuint textureID, GLuint samplerHandle, GLuint target)
{
    if (pullShaderProgram != NULL) {
        pullShaderProgram->texture(name, textureID, samplerHandle, target);
    }
    if (pushShaderProgram != NULL) {
        pushShaderProgram->texture(name, textureID, samplerHandle, target);
    }
    return this;
}

Pyramid* Pyramid::clear(float r, float g, float b, float a, int level) {
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
    clear(0, 0, 0, 0, level);
    return this;
}

GLuint Pyramid::get(std::string name) {
    return textureMap[name];
}

int Pyramid::getMipmapNumber() {
    return fboHandles.size();
}
