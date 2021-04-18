/*
 * SimplestSphereRenderer.cpp
 *
 * Copyright (C) 2018 by Karsten Schatz
 * Copyright (C) 2018 by VISUS (Universitaet Stuttgart)
 * Alle Rechte vorbehalten.
 */
#include "stdafx.h"
#include "SimpleContourRenderer.h"
#include "geometry_calls/CallTriMeshData.h"
#include "mmcore/CoreInstance.h"
#include "mmcore/view/CallRender3DGL.h"
#include "mmcore/view/Camera_2.h"
#include "mmcore/param/BoolParam.h"
#include "mmcore/param/FloatParam.h"
#include "mmcore/param/IntParam.h"
#include "vislib/math/Matrix.h"
#include "vislib/math/ShallowMatrix.h"

using namespace megamol;
using namespace megamol::contours;

/*
 * SimplestSphereRenderer::SimplestSphereRenderer
 */
SimpleContourRenderer::SimpleContourRenderer(void) : core::view::Renderer3DModuleGL()
    , getDataSlot("getData", "The input data slot for triangulated mesh data.")
    , sizeScalingSlot("scaling factor", "Scaling factor for the size of the contour") {
    // TUTORIAL: A name and a description for each slot (CallerSlot, CalleeSlot, ParamSlot) has to be given in the
    // constructor initializer list

    // TUTORIAL: For each CallerSlot all compatible calls have to be set
    this->getDataSlot.SetCompatibleCall<megamol::geocalls::CallTriMeshDataDescription>();
    this->MakeSlotAvailable(&this->getDataSlot);

    this->sizeScalingSlot.SetParameter(new core::param::FloatParam(1.0f, 0.01f, 1000.0f));
    this->MakeSlotAvailable(&this->sizeScalingSlot);

    // TUTORIAL: Each slot that shall be visible in the GUI has to be made available by this->MakeSlotAvailable(...)

    VBO = 0;
    VAO = 0;
}

/*
 * SimplestSphereRenderer::~SimplestSphereRenderer
 */
SimpleContourRenderer::~SimpleContourRenderer(void) {
    this->Release();
    // TUTORIAL: this->Release() should be called in each modules' destructor.
}

/*
 * SimplestSphereRenderer::create
 */
bool SimpleContourRenderer::create(void) {

    // TUTORIAL Shader creation should always happen in the create method of a renderer.

    using namespace megamol::core::utility::log;
    using namespace vislib::graphics::gl;

    ShaderSource vertSrc;
    ShaderSource fragSrc;
    if (!this->GetCoreInstance()->ShaderSourceFactory().MakeShaderSource("simpleMesh::vertex", vertSrc)) {
        Log::DefaultLog.WriteMsg(Log::LEVEL_ERROR, "Unable to load vertex shader source for simple mesh shader");
        return false;
    }
    if (!this->GetCoreInstance()->ShaderSourceFactory().MakeShaderSource("simpleMesh::fragment", fragSrc)) {
        Log::DefaultLog.WriteMsg(Log::LEVEL_ERROR, "Unable to load fragment shader source for simple mesh shader");
        return false;
    }
    try {
        if (!this->simpleShader.Create(vertSrc.Code(), vertSrc.Count(), fragSrc.Code(), fragSrc.Count())) {
            throw vislib::Exception("Generic creation failure", __FILE__, __LINE__);
        }
    } catch (vislib::Exception e) {
        Log::DefaultLog.WriteMsg(Log::LEVEL_ERROR, "Unable to create mesh shader: %s\n", e.GetMsgA());
        return false;
    }
    return true;
}

/*
 * SimplestSphereRenderer::GetExtents
 */
bool SimpleContourRenderer::GetExtents(core::view::CallRender3DGL& call) {

    megamol::geocalls::CallTriMeshData *ctmd = this->getDataSlot.CallAs<megamol::geocalls::CallTriMeshData>();
    if (ctmd == NULL) return false;
    ctmd->SetFrameID(static_cast<int>(call.Time()));
    if (!(*ctmd)(1)) return false;

    call.SetTimeFramesCount(ctmd->FrameCount());
    call.AccessBoundingBoxes().Clear();
    call.AccessBoundingBoxes() = ctmd->AccessBoundingBoxes();

    return true;
}

/*
 * SimplestSphereRenderer::release
 */
void SimpleContourRenderer::release(void) {
    if (VAO != 0) {
        glDeleteVertexArrays(1, &VAO);
    }
    if (VBO != 0) {
        glDeleteBuffers(1, &VBO);
    }
}

/*
 * SimplestSphereRenderer::Render
 */
bool SimpleContourRenderer::Render(core::view::CallRender3DGL& call) {
    core::view::CallRender3DGL* cr3d = dynamic_cast<core::view::CallRender3DGL*>(&call);
    if (cr3d == nullptr) return false;
    
    megamol::geocalls::CallTriMeshData *ctmd = this->getDataSlot.CallAs<megamol::geocalls::CallTriMeshData>();
    if (ctmd == NULL) return false;

    ctmd->SetFrameID(static_cast<int>(call.Time()));
    if (!(*ctmd)(1)) return false;

    ctmd->SetFrameID(static_cast<int>(call.Time())); // necessary?
    if (!(*ctmd)(0)) return false;
    for (unsigned int i = 0; i < ctmd->Count(); i++) {
        const megamol::geocalls::CallTriMeshData::Mesh& obj = ctmd->Objects()[i];
        auto vertCount = obj.GetVertexCount();


        // if (obj.HasNormalPointer() != NULL) {
        //     switch (obj.GetNormalDataType()) {
        //         case megamol::geocalls::CallTriMeshData::Mesh::DT_FLOAT:
        //             ::glNormalPointer(GL_FLOAT, 0, obj.GetNormalPointerFloat());
        //             break;
        //         case megamol::geocalls::CallTriMeshData::Mesh::DT_DOUBLE:
        //             ::glNormalPointer(GL_DOUBLE, 0, obj.GetNormalPointerDouble());
        //             break;
        //         default: continue;
        //     }
        // }

    
    glGenVertexArrays(1, &VAO);
    glGenBuffers(1, &VBO);
    glBindVertexArray(VAO);
    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glEnableVertexAttribArray(0);

    switch (obj.GetVertexDataType()) {
        case megamol::geocalls::CallTriMeshData::Mesh::DT_FLOAT:
            glBufferData(GL_ARRAY_BUFFER, sizeof(float) * 3 * vertCount, obj.GetVertexPointerFloat(),
                GL_STATIC_DRAW); 
            glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(float) * 3, 0);
            break;
        case megamol::geocalls::CallTriMeshData::Mesh::DT_DOUBLE:
            glBufferData(GL_ARRAY_BUFFER, sizeof(double) * 3 * vertCount, obj.GetVertexPointerDouble(),
                GL_STATIC_DRAW); 
            glVertexAttribPointer(0, 3, GL_DOUBLE, GL_FALSE, sizeof(double) * 3, 0);
            break;
        default: continue;
    }
    
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindVertexArray(0);

    core::view::Camera_2 localCam;
    cr3d->GetCamera(localCam);

    cam_type::snapshot_type camsnap;
    cam_type::matrix_type viewCam, projCam;
    localCam.calc_matrices(camsnap, viewCam, projCam);

    glm::mat4 view = viewCam;
    glm::mat4 proj = projCam;
    glm::mat4 mvp = projCam * viewCam;

    // start the rendering

    this->simpleShader.Enable();

    glBindVertexArray(VAO);
    glEnableVertexAttribArray(0);

    glUniformMatrix4fv(this->simpleShader.ParameterLocation("mvp"), 1, GL_FALSE, glm::value_ptr(mvp));

    glEnable(GL_DEPTH_TEST);

    if (obj.HasTriIndexPointer() != NULL) {
        switch (obj.GetTriDataType()) {
            case megamol::geocalls::CallTriMeshData::Mesh::DT_BYTE:
                ::glDrawElements(GL_TRIANGLES, obj.GetTriCount() * 3, GL_UNSIGNED_BYTE, obj.GetTriIndexPointerByte());
                break;
            case megamol::geocalls::CallTriMeshData::Mesh::DT_UINT16:
                ::glDrawElements(GL_TRIANGLES, obj.GetTriCount() * 3, GL_UNSIGNED_SHORT, obj.GetTriIndexPointerUInt16());
                break;
            case megamol::geocalls::CallTriMeshData::Mesh::DT_UINT32:
                ::glDrawElements(GL_TRIANGLES, obj.GetTriCount() * 3, GL_UNSIGNED_INT, obj.GetTriIndexPointerUInt32());
                break;
            default: continue;
        }
    } else {
        glDrawArrays(GL_TRIANGLES, 0, obj.GetVertexCount());
    }
    glBindVertexArray(0);

    glDisable(GL_DEPTH_TEST);

    this->simpleShader.Disable();
    }
    return true;
}
