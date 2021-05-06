/*
 * View3DGL.cpp
 *
 * Copyright (C) 2018, 2020 by VISUS (Universitaet Stuttgart).
 * Alle Rechte vorbehalten.
 */

#include "stdafx.h"
#include "mmcore/view/View3DGL.h"

#include "GlobalValueStore.h"
#include "mmcore/param/BoolParam.h"
#include "mmcore/CoreInstance.h"
#include "mmcore/param/FloatParam.h"
#include "mmcore/view/CallRender3DGL.h"
#include "mmcore/view/CallRenderViewGL.h"

using namespace megamol::core;
using namespace megamol::core::view;

/*
 * View3DGL::View3DGL
 */
View3DGL::View3DGL() : view::BaseView<CallRenderViewGL, Camera3DController>() {
    this->_rhsRenderSlot.SetCompatibleCall<CallRender3DGLDescription>();
    this->MakeSlotAvailable(&this->_rhsRenderSlot);
    // Override renderSlot behavior
    this->_lhsRenderSlot.SetCallback(
        view::CallRenderViewGL::ClassName(), InputCall::FunctionName(InputCall::FnOnKey), &AbstractView::OnKeyCallback);
    this->_lhsRenderSlot.SetCallback(view::CallRenderViewGL::ClassName(), InputCall::FunctionName(InputCall::FnOnChar),
        &AbstractView::OnCharCallback);
    this->_lhsRenderSlot.SetCallback(view::CallRenderViewGL::ClassName(),
        InputCall::FunctionName(InputCall::FnOnMouseButton), &AbstractView::OnMouseButtonCallback);
    this->_lhsRenderSlot.SetCallback(view::CallRenderViewGL::ClassName(),
        InputCall::FunctionName(InputCall::FnOnMouseMove), &AbstractView::OnMouseMoveCallback);
    this->_lhsRenderSlot.SetCallback(view::CallRenderViewGL::ClassName(),
        InputCall::FunctionName(InputCall::FnOnMouseScroll), &AbstractView::OnMouseScrollCallback);
    // AbstractCallRender
    this->_lhsRenderSlot.SetCallback(view::CallRenderViewGL::ClassName(),
        AbstractCallRender::FunctionName(AbstractCallRender::FnRender), &AbstractView::OnRenderView);
    this->_lhsRenderSlot.SetCallback(view::CallRenderViewGL::ClassName(),
        AbstractCallRender::FunctionName(AbstractCallRender::FnGetExtents), &AbstractView::GetExtents);
    // CallRenderViewGL
    this->_lhsRenderSlot.SetCallback(view::CallRenderViewGL::ClassName(),
        view::CallRenderViewGL::FunctionName(view::CallRenderViewGL::CALL_RESETVIEW), &AbstractView::OnResetView);
    this->MakeSlotAvailable(&this->_lhsRenderSlot);

    this->_rhsRenderSlot.SetNecessity(megamol::core::AbstractCallSlotPresentation::SLOT_REQUIRED);
}

/*
 * View3DGL::~View3DGL
 */
View3DGL::~View3DGL() {
    this->Release();
}

ImageWrapper megamol::core::view::View3DGL::Render(double time, double instanceTime, bool present_fbo) {
    CallRender3DGL* cr3d = this->_rhsRenderSlot.CallAs<CallRender3DGL>();

    if (cr3d != NULL) {

        BaseView::beforeRender(time, instanceTime);

        // clear fbo before sending it down the rendering call
        // the view is the owner of this fbo and therefore responsible
        // for clearing it at the beginning of a render frame
        _fbo->bind();
        auto bgcol = this->BkgndColour();
        glClearColor(bgcol.r, bgcol.g, bgcol.b, bgcol.a);
        glClearDepth(1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        glBindFramebuffer(GL_FRAMEBUFFER, 0);

        // set camera and fbo in rendering call
        cr3d->SetFramebuffer(_fbo);
        cr3d->SetCamera(this->_camera);

        // call the rendering call
        (*cr3d)(view::CallRender3DGL::FnRender);

        BaseView::afterRender();
    }

    if (present_fbo) {
        // Blit the final image to the default framebuffer of the window.
        // Technically, the view's fbo should always match the size of the window so a blit is fine.
        // Eventually, presenting the fbo will become the frontends job.
        // Bind and blit framebuffer.
        glBindFramebuffer(GL_FRAMEBUFFER, 0);

        _fbo->bindToRead(0);
        glBlitFramebuffer(0, 0, _fbo->getWidth(), _fbo->getHeight(), 0, 0, _fbo->getWidth(), _fbo->getHeight(),
            GL_COLOR_BUFFER_BIT, GL_NEAREST);

        glBindFramebuffer(GL_READ_FRAMEBUFFER, 0);
    }

    return GetRenderingResult();
}

ImageWrapper megamol::core::view::View3DGL::GetRenderingResult() const {
    ImageWrapper::DataChannels channels =
        ImageWrapper::DataChannels::RGBA8; // vislib::graphics::gl::FramebufferObject seems to use RGBA8
    unsigned int fbo_color_buffer_gl_handle =
        _fbo->getColorAttachment(0)->getTextureHandle(); // IS THIS SAFE?? IS THIS THE COLOR BUFFER??
    size_t fbo_width = _fbo->getWidth();
    size_t fbo_height = _fbo->getHeight();

    return frontend_resources::wrap_image({fbo_width, fbo_height}, fbo_color_buffer_gl_handle, channels);
}

void megamol::core::view::View3DGL::Resize(unsigned int width, unsigned int height) {

    BaseView::Resize(width, height);

    bool create_fbo = false;
    if (_fbo == nullptr) {
        create_fbo = true;
    } else if ((_fbo->getWidth() != width) || (_fbo->getHeight() != height)) {
        create_fbo = true;
    }

    if (create_fbo) {
        glBindFramebuffer(GL_FRAMEBUFFER, 0); // better safe then sorry, "unbind" fbo before delting one
        try {
            _fbo = std::make_shared<glowl::FramebufferObject>(width, height);
            _fbo->createColorAttachment(GL_RGBA8, GL_RGBA, GL_UNSIGNED_BYTE);
            // TODO: check completness and throw if not?
        } catch (glowl::FramebufferObjectException const& exc) {
            megamol::core::utility::log::Log::DefaultLog.WriteError(
                "[View3DGL] Unable to create framebuffer object: %s\n", exc.what());
        }
    }
}

/*
 * View3DGL::create
 */
bool View3DGL::create() {
    // intialize fbo with dummy size until the actual size is set during first call to Resize
    this->_fbo = std::make_shared<glowl::FramebufferObject>(1,1);

    const auto arcball_key = "arcball";

    if (!this->GetCoreInstance()->IsmmconsoleFrontendCompatible()) {
        // new frontend has global key-value resource
        auto maybe = this->frontend_resources[0].getResource<megamol::frontend_resources::GlobalValueStore>().maybe_get(
            arcball_key);
        if (maybe.has_value()) {
            this->_camera_controller.setArcballDefault(vislib::CharTraitsA::ParseBool(maybe.value().c_str()));
        }

    } else {
        mmcValueType wpType;
        this->_camera_controller.setArcballDefault(false);
        auto value = this->GetCoreInstance()->Configuration().GetValue(MMC_CFGID_VARIABLE, _T(arcball_key), &wpType);
        if (value != nullptr) {
            try {
                switch (wpType) {
                case MMC_TYPE_BOOL:
                    this->_camera_controller.setArcballDefault(*static_cast<const bool*>(value));
                    break;

                case MMC_TYPE_CSTR:
                    this->_camera_controller.setArcballDefault(
                        vislib::CharTraitsA::ParseBool(static_cast<const char*>(value)));
                    break;

                case MMC_TYPE_WSTR:
                    this->_camera_controller.setArcballDefault(
                        vislib::CharTraitsW::ParseBool(static_cast<const wchar_t*>(value)));
                    break;
                }
            } catch (...) {}
        }
    }

    this->_firstImg = true;

    return true;
}
