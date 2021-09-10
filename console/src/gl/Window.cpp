/*
 * Window.cpp
 *
 * Copyright (C) 2008, 2016 MegaMol Team
 * Alle Rechte vorbehalten.
 */

#include "stdafx.h"
#include "gl/Window.h"
#include "utility/HotFixes.h"
#include "mmcore/utility/log/Log.h"
#include "WindowManager.h"
#include <cassert>
#include <algorithm>
#include <sstream>
#include <string>
#include <iostream>
#include <fstream>
#ifdef _WIN32
#define GLFW_EXPOSE_NATIVE_WIN32 
#include "GLFW/glfw3native.h" 
#endif

#include "mmcore/view/Input.h"

//#include "HotKeyButtonParam.h"
//#include "vislib/RawStorage.h"
//#include "vislib/types.h"
//#include "WindowManager.h"
//#include <chrono>
//#include <vector>
//#include <algorithm>

using namespace megamol;
using namespace megamol::console;

gl::Window::Window(const char* title, const utility::WindowPlacement & placement, GLFWwindow* share)
        : glfw(), 
        hView(), hWnd(nullptr), width(-1), height(-1), uiLayers(), mouseCapture(),
        name(title), fpsCntr(), fps(1000.0f), fpsList(), showFpsInTitle(true), fpsSyncTime(), topMost(false),
        fragmentQuery(0), showFragmentsInTitle(false), showPrimsInTitle(false), firstUpdate(true) {

    if (::memcmp(name.c_str(), WindowManager::TitlePrefix, WindowManager::TitlePrefixLength) == 0) {
        name = name.substr(WindowManager::TitlePrefixLength);
    }
    for (float& f : fpsList) f = 0.0f;

    glfw = glfwInst::Instance(); // we use glfw
    if (glfw->OK()) {
        if (utility::HotFixes::Instance().IsHotFixed("usealphabuffer")) {
            ::glfwWindowHint(GLFW_ALPHA_BITS, 8);
        }

//#ifndef NOWINDOWPOSFIX
//        if (wndX != predictedX || wndY != predictedY ||
//            wndW != predictedWidth || wndH != predictedHeight) {
//            Log::DefaultLog.WriteMsg(Log::LEVEL_WARN, "The actual "
//                "view window location reported by the core (%d, %d), "
//                "size (%d, %d) is "
//                "different from the one predicted. GPU affinity "
//                "may have been set incorrectly.", wndX, wndY, wndW,
//                wndH);
//        }
//#endif

        this->topMost = placement.topMost;
        
        if (!placement.fullScreen) {
            // window mode
            ::glfwWindowHint(GLFW_DECORATED, placement.noDec ? GLFW_FALSE : GLFW_TRUE);
            ::glfwWindowHint(GLFW_VISIBLE, GLFW_FALSE);

            // Hack for window size larger than screen:
            // The width and height for glfwCreateWindow is just a hint for the underlying window manager. Most window
            // managers (Windows + Linux/X11) will just resize and maximise the window on the screen, when a size
            // larger than the screen is requested. When we do not allow resizing, the window managers seems to use the
            // wanted size.
            // But we want a resizable window, therefore we need to set the window resizable again. But the second
            // problem with window managers is, that window creation is an async task (at least with X11). When we
            // immediately after window creation set the window to be resizable again, the automatic resize will happen
            // again. We need to delay this to a later point in time (i.e. rendering of first frame).
            // The async task problem holds (of course) also for querying the window size with glfwGetWindowSize.
            // Calling this right after glfwCreateWindow will just return our values, not the actual window size.
            // Therefore here we cannot test for correct window size.
            ::glfwWindowHint(GLFW_RESIZABLE, GLFW_FALSE);

            int w = placement.w;
            int h = placement.h;
            if (!placement.size || (w <= 0) || (h <= 0)) {
                GLFWmonitor* primary = ::glfwGetPrimaryMonitor();
                const GLFWvidmode* mode = ::glfwGetVideoMode(primary);
                w = mode->width * 3 / 4;
                h = mode->height * 3 / 4;
            }

            // According to glfw docs width and height means the content area of the window without window decorations.
            hWnd = ::glfwCreateWindow(w, h, title, nullptr, share);
            megamol::core::utility::log::Log::DefaultLog.WriteInfo("Console::Window: Create window with size w: %d, h: %d\n", w, h);
            if (hWnd != nullptr) {
                if (placement.pos) ::glfwSetWindowPos(hWnd, placement.x, placement.y);
            }

        } else {
            // fullscreen mode
            int monCnt = 0;
            GLFWmonitor **mons = ::glfwGetMonitors(&monCnt);
            GLFWmonitor *mon = mons[std::min<int>(monCnt - 1, placement.mon)];
            const GLFWvidmode* mode = glfwGetVideoMode(mon);

            if (placement.pos) megamol::core::utility::log::Log::DefaultLog.WriteWarn("Ignoring window placement position when requesting fullscreen.");
            if (placement.size) {
                if ((placement.w != mode->width) || (placement.h != mode->height)) {
                    megamol::core::utility::log::Log::DefaultLog.WriteWarn("Changing screen resolution is currently not supported.");
                }
            }
            if (placement.noDec) megamol::core::utility::log::Log::DefaultLog.WriteWarn("Ignoring no-decorations setting when requesting fullscreen.");

            ::glfwWindowHint(GLFW_DECORATED, GLFW_FALSE);
            ::glfwWindowHint(GLFW_VISIBLE, GLFW_FALSE);
            ::glfwWindowHint(GLFW_RED_BITS, mode->redBits);
            ::glfwWindowHint(GLFW_GREEN_BITS, mode->greenBits);
            ::glfwWindowHint(GLFW_BLUE_BITS, mode->blueBits);
            ::glfwWindowHint(GLFW_REFRESH_RATE, mode->refreshRate);
            // this only works since we are NOT setting a monitor
            ::glfwWindowHint(GLFW_FLOATING, GLFW_TRUE);

            /* note we do not use a real fullscreen mode, since then we would have focus-iconify problems */
            hWnd = ::glfwCreateWindow(mode->width, mode->height, title, nullptr, share);
            megamol::core::utility::log::Log::DefaultLog.WriteInfo("Console::Window: Create window with size w: %d, h: %d\n", mode->width, mode->height);
            int x, y;
            ::glfwGetMonitorPos(mon, &x, &y);
            ::glfwSetWindowPos(hWnd, x, y);
        }

        if (hWnd != nullptr) {
            ::glfwSetWindowUserPointer(hWnd, this); // this is ok, as long as no one derives from Window at this point
            ::glfwShowWindow(hWnd);
            ::glfwMakeContextCurrent(hWnd);
            if ((placement.fullScreen || placement.noDec) && (!utility::HotFixes::Instance().IsHotFixed("DontHideCursor"))) {
                ::glfwSetInputMode(hWnd, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
            }
            vislib::graphics::gl::LoadAllGL();
            ::glfwSetKeyCallback(hWnd, &Window::glfw_onKey_func);
            ::glfwSetMouseButtonCallback(hWnd, &Window::glfw_onMouseButton_func);
            ::glfwSetCursorPosCallback(hWnd, &Window::glfw_onMouseMove_func);
            ::glfwSetScrollCallback(hWnd, &Window::glfw_onMouseWheel_func);
            ::glfwSetCharCallback(hWnd, &Window::glfw_onChar_func);

            GLint vp[4];
            glGetIntegerv(GL_VIEWPORT, vp);
            megamol::core::utility::log::Log::DefaultLog.WriteInfo("Console::Window: viewport size w: %d, h: %d\n", vp[2], vp[3]);

        } else {
            megamol::core::utility::log::Log::DefaultLog.WriteError(
                "Could not create GLFW Window. You probably do not have OpenGL support. Your graphics hardware might "
                "be very old, your drivers could be outdated or you are running in a remote desktop session.");
            // we should do a proper shutdown now, but that is too expensive given the expected lifetime of this front end.
            exit(-1);
        }

        glGenQueries(1, &fragmentQuery);
        glGenQueries(1, &primsQuery);
        fpsSyncTime = std::chrono::system_clock::now();
    }
}

gl::Window::~Window() {
    assert(hWnd == nullptr);
    glDeleteQueries(1, &fragmentQuery);
    glDeleteQueries(1, &primsQuery);
}

void gl::Window::EnableVSync() {
    if (hWnd != nullptr) {
        ::glfwMakeContextCurrent(hWnd);
        ::glfwSwapInterval(0);
    }
}

void gl::Window::AddUILayer(std::shared_ptr<AbstractUILayer> uiLayer) {
    uiLayers.AddUILayer(uiLayer);
}

void gl::Window::RemoveUILayer(std::shared_ptr<AbstractUILayer> uiLayer) {
    uiLayers.RemoveUILayer(uiLayer);
}

void gl::Window::SetShowFPSinTitle(bool show) {
    showFpsInTitle = show;
    if (!showFpsInTitle && !showFragmentsInTitle) {
        ::glfwSetWindowTitle(hWnd, (std::string(WindowManager::TitlePrefix) + name).c_str());
    }
}

void gl::Window::SetShowSamplesinTitle(bool show) {
    showFragmentsInTitle = show;
    if (!showFpsInTitle && !showFragmentsInTitle) {
        ::glfwSetWindowTitle(hWnd, (std::string(WindowManager::TitlePrefix) + name).c_str());
    }
}

void gl::Window::SetShowPrimsinTitle(bool show) {
    showPrimsInTitle = show;
    if (!showFpsInTitle && !showFragmentsInTitle && !showPrimsInTitle) {
        ::glfwSetWindowTitle(hWnd, (std::string(WindowManager::TitlePrefix) + name).c_str());
    }
}

void gl::Window::RequestClose() {
    if (hWnd != nullptr) {
        ::glfwSetWindowShouldClose(hWnd, true);
    }
}

void gl::Window::Update(uint32_t frameID) {
    if (hWnd == nullptr) return;

    // See comment above window creation on window size.
    if (firstUpdate) {
        firstUpdate = false;

        int actualw = 0, actualh = 0;
        glfwGetWindowSize(hWnd, &actualw, &actualh);
        megamol::core::utility::log::Log::DefaultLog.WriteInfo("Console::Window: Actual window size: w: %d, h: %d\n", actualw, actualh);

        glfwSetWindowAttrib(hWnd, GLFW_RESIZABLE, GLFW_TRUE);
    }

    // this also issues the callbacks, which might close this window
    ::glfwPollEvents();

    if (hWnd == nullptr) return;
    if (::glfwWindowShouldClose(hWnd)) {
        uiLayers.ClearUILayers();

        hView.DestroyHandle();

        ::glfwDestroyWindow(hWnd);
        hWnd = nullptr;
        return;
    }

    ::glfwMakeContextCurrent(hWnd);
    int frame_width, frame_height;
    ::glfwGetFramebufferSize(hWnd, &frame_width, &frame_height);
    if ((frame_width != width) || (frame_height != height)) {
        on_resize(frame_width, frame_height);
        width = frame_width;
        height = frame_height;
    }

    fpsCntr.FrameBegin();
    if ((width > 0) && (height > 0)) {
        if (showFragmentsInTitle) glBeginQuery(GL_SAMPLES_PASSED, fragmentQuery);
        if (showPrimsInTitle) glBeginQuery(GL_PRIMITIVES_GENERATED, primsQuery);
        ::mmcRenderView(hView, frameID);
        if (showFragmentsInTitle) glEndQuery(GL_SAMPLES_PASSED);
        if (showPrimsInTitle) glEndQuery(GL_PRIMITIVES_GENERATED);
    }

	this->uiLayers.OnDraw();

    // done rendering. swap and next turn
    ::glfwSwapBuffers(hWnd);
    fpsCntr.FrameEnd();

    std::chrono::system_clock::time_point now = std::chrono::system_clock::now();
    if (now - fpsSyncTime > std::chrono::seconds(1)) {
        on_fps_value(fpsCntr.FPS());
        fpsSyncTime = now;
#ifdef _WIN32
        // TODO fix this for EGL + Win
        if (this->topMost) {
            megamol::core::utility::log::Log::DefaultLog.WriteInfo("Periodic reordering of windows.");
            SetWindowPos(glfwGetWin32Window(this->hWnd), HWND_TOPMOST, 0, 0, 0, 0, SWP_NOSIZE | SWP_NOMOVE);
        }
#endif
    }

}

void gl::Window::glfw_onKey_func(GLFWwindow* wnd, int k, int s, int a, int m) {
    ::glfwMakeContextCurrent(wnd);
    Window* that = static_cast<Window*>(::glfwGetWindowUserPointer(wnd));

    core::view::Key key = static_cast<core::view::Key>(k);

    core::view::KeyAction action(core::view::KeyAction::RELEASE);
    switch (a) {
    case GLFW_PRESS: action = core::view::KeyAction::PRESS; break;
    case GLFW_REPEAT: action = core::view::KeyAction::REPEAT; break;
    case GLFW_RELEASE: action = core::view::KeyAction::RELEASE; break;
    }

    core::view::Modifiers mods;
    // Parameter m is not platform independent, see https://github.com/glfw/glfw/issues/1630.
    // We want here consistent modifier state, in the sense that the modifier is evaluated on the keyboard state after
    // the current key event. A simple solution for this would be to just check glfwGetKey() key here and set the
    // modifiers based on this. The problem with this is, that glfwGetKey only returns a cached state and therefore
    // does not know when the modifier key was pressed outside of the window, when it is not in focus.
    // The new solution is, that we keep using the GLFW modifier state, but in addition we will evaluate the current
    // event to overwrite modifiers when have a event the modifier keys itself.
    if ((m & GLFW_MOD_SHIFT) == GLFW_MOD_SHIFT) mods |= core::view::Modifier::SHIFT;
    if ((m & GLFW_MOD_CONTROL) == GLFW_MOD_CONTROL) mods |= core::view::Modifier::CTRL;
    if ((m & GLFW_MOD_ALT) == GLFW_MOD_ALT) mods |= core::view::Modifier::ALT;

    if (k == GLFW_KEY_LEFT_SHIFT || k == GLFW_KEY_RIGHT_SHIFT) {
        if (a == GLFW_RELEASE) {
            mods.reset(core::view::Modifier::SHIFT);
        } else {
            mods.set(core::view::Modifier::SHIFT);
        }
    }
    if (k == GLFW_KEY_LEFT_CONTROL || k == GLFW_KEY_RIGHT_CONTROL) {
        if (a == GLFW_RELEASE) {
            mods.reset(core::view::Modifier::CTRL);
        } else {
            mods.set(core::view::Modifier::CTRL);
        }
    }
    if (k == GLFW_KEY_LEFT_ALT || k == GLFW_KEY_RIGHT_ALT) {
        if (a == GLFW_RELEASE) {
            mods.reset(core::view::Modifier::ALT);
        } else {
            mods.set(core::view::Modifier::ALT);
        }
    }

    that->uiLayers.OnKey(key, action, mods);
}

void gl::Window::glfw_onChar_func(GLFWwindow* wnd, unsigned int charcode) {
    ::glfwMakeContextCurrent(wnd);
    Window* that = static_cast<Window*>(::glfwGetWindowUserPointer(wnd));
    that->uiLayers.OnChar(charcode);
}

void gl::Window::glfw_onMouseMove_func(GLFWwindow* wnd, double x, double y) {
    ::glfwMakeContextCurrent(wnd);
    Window* that = static_cast<Window*>(::glfwGetWindowUserPointer(wnd));
    if (that->mouseCapture) {
        that->mouseCapture->OnMouseMove(x, y);
    } else {
        that->uiLayers.OnMouseMove(x, y);
    }
}

void gl::Window::glfw_onMouseButton_func(GLFWwindow* wnd, int b, int a, int m) {
    ::glfwMakeContextCurrent(wnd);
    Window* that = static_cast<Window*>(::glfwGetWindowUserPointer(wnd));

    core::view::MouseButton btn = static_cast<core::view::MouseButton>(b);

    core::view::MouseButtonAction action =
        (a == GLFW_PRESS) ? core::view::MouseButtonAction::PRESS
            : core::view::MouseButtonAction::RELEASE;

    core::view::Modifiers mods;
    // Parameter m is not platform independent, see https://github.com/glfw/glfw/issues/1630.
    // But this should only play a role on the key events of the modifier keys itself.
    // Therefore we can ignore this here. See comment in keyboard callback for more details.
    if ((m & GLFW_MOD_SHIFT) == GLFW_MOD_SHIFT) mods |= core::view::Modifier::SHIFT;
    if ((m & GLFW_MOD_CONTROL) == GLFW_MOD_CONTROL) mods |= core::view::Modifier::CTRL;
    if ((m & GLFW_MOD_ALT) == GLFW_MOD_ALT) mods |= core::view::Modifier::ALT;

    if (that->mouseCapture) {
        that->mouseCapture->OnMouseButton(btn, action, mods);
    } else {
		if(that->uiLayers.OnMouseButton(btn, action, mods))
			if (action == core::view::MouseButtonAction::PRESS)
				that->mouseCapture = that->uiLayers.lastEventCaptureUILayer();
    }

    if (that->mouseCapture) {
        bool anyPressed = false;
        for (int mbi = GLFW_MOUSE_BUTTON_1; mbi <= GLFW_MOUSE_BUTTON_LAST; ++mbi) {
            if (::glfwGetMouseButton(wnd, mbi) == GLFW_PRESS) {
                anyPressed = true;
                break;
            }
        }
        if (!anyPressed) {
            that->mouseCapture.reset();
            double x, y;
            ::glfwGetCursorPos(wnd, &x, &y);
            glfw_onMouseMove_func(wnd, x, y); // to inform all of the new location
        }
    }
}

void gl::Window::glfw_onMouseWheel_func(GLFWwindow* wnd, double x, double y) {
    ::glfwMakeContextCurrent(wnd);
    Window* that = static_cast<Window*>(::glfwGetWindowUserPointer(wnd));
    if (that->mouseCapture) {
        that->mouseCapture->OnMouseScroll(x, y);
    } else {
        that->uiLayers.OnMouseScroll(x, y);
    }
}

void gl::Window::on_resize(int w, int h) {
    ::glfwMakeContextCurrent(hWnd);
    if ((w > 0) && (h > 0)) {
        ::glViewport(0, 0, w, h);
        ::mmcResizeView(hView, w, h);
        megamol::core::utility::log::Log::DefaultLog.WriteInfo("Console::Window: Resize window (w: %d, h: %d)\n", w, h);
		uiLayers.OnResize(w, h);
    }
}

void gl::Window::on_fps_value(float fps_val) {
    fps = fps_val;

    auto i1 = fpsList.begin();
    auto i2 = i1 + 1;
    auto e = fpsList.end();
    while (i2 != e) {
        *i1 = *i2;
        ++i1;
        ++i2;
    }
    fpsList[fpsList.size() - 1] = fps;

    //if (showFpsInTitle) {
    //    std::stringstream title;
    //    if (showFragmentsInTitle) {
    //        GLuint64 samp;
    //        glGetQueryObjectui64v(fragmentQuery, GL_QUERY_RESULT, &samp);
    //        title << WindowManager::TitlePrefix << name << " - [" << fps << " fps, " << samp << " samples]";
    //    } else {
    //        title << WindowManager::TitlePrefix << name << " - [" << fps << " fps]";
    //    }
    //    ::glfwSetWindowTitle(hWnd, title.str().c_str());
    //}
    std::stringstream title;
    title.imbue(std::locale(""));
    title << WindowManager::TitlePrefix << name;
    if (showFpsInTitle || showFragmentsInTitle || showPrimsInTitle) title << " - [ ";
    if (showFpsInTitle) {
        title << fps << " fps ";
    }
    if (showFragmentsInTitle) {
        GLuint64 samp;
        glGetQueryObjectui64v(fragmentQuery, GL_QUERY_RESULT, &samp);
        title << samp << " samples ";
    }
    if (showPrimsInTitle) {
        GLuint64 prims;
        glGetQueryObjectui64v(primsQuery, GL_QUERY_RESULT, &prims);
        title << prims << " primitives ";
    }
    if (showFpsInTitle || showFragmentsInTitle || showPrimsInTitle) title << "]";
    ::glfwSetWindowTitle(hWnd, title.str().c_str());
}

