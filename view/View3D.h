/*
 * View3D.h
 *
 * Copyright (C) 2008 - 2010 by VISUS (Universitaet Stuttgart). 
 * Alle Rechte vorbehalten.
 */

#ifndef MEGAMOLCORE_VIEW3D_H_INCLUDED
#define MEGAMOLCORE_VIEW3D_H_INCLUDED
#if (defined(_MSC_VER) && (_MSC_VER > 1000))
#pragma once
#endif /* (defined(_MSC_VER) && (_MSC_VER > 1000)) */

#define ENABLE_KEYBOARD_VIEW_CONTROL 1

#include "BoundingBoxes.h"
#include "view/AbstractCallRender.h"
#include "view/AbstractView3D.h"
#include "CalleeSlot.h"
#include "CallerSlot.h"
#include "param/ParamSlot.h"
#include "vislib/CameraLookAtDist.h"
#include "vislib/CameraOpenGL.h"
#include "vislib/CameraParamsStore.h"
#include "vislib/CameraRotate2D.h"
#include "vislib/CameraRotate2DLookAt.h"
#include "vislib/CameraZoom2DMove.h"
#include "vislib/CameraZoom2DAngle.h"
#include "vislib/Cursor2D.h"
#include "vislib/FpsCounter.h"
#include "vislib/graphicstypes.h"
#include "vislib/InputModifiers.h"
#include "vislib/Cuboid.h"
#include "vislib/PerformanceCounter.h"


namespace megamol {
namespace core {
namespace view {


    /**
     * Base class of rendering graph calls
     */
    class View3D: public AbstractView3D {
    public:

        /**
         * Answer the name of this module.
         *
         * @return The name of this module.
         */
        static const char *ClassName(void) {
            return "View3D";
        }

        /**
         * Answer a human readable description of this module.
         *
         * @return A human readable description of this module.
         */
        static const char *Description(void) {
            return "3D View Module";
        }

        /**
         * Answers whether this module is available on the current system.
         *
         * @return 'true' if the module is available, 'false' otherwise.
         */
        static bool IsAvailable(void) {
            return true;
        }

        /** Ctor. */
        View3D(void);

        /** Dtor. */
        virtual ~View3D(void);

        /**
         * Answer the camera synchronization number.
         *
         * @return The camera synchronization number
         */
        virtual unsigned int GetCameraSyncNumber(void) const;

        /**
         * Serialises the camera of the view
         *
         * @param serialiser Serialises the camera of the view
         */
        virtual void SerialiseCamera(vislib::Serialiser& serialiser) const;

        /**
         * Deserialises the camera of the view
         *
         * @param serialiser Deserialises the camera of the view
         */
        virtual void DeserialiseCamera(vislib::Serialiser& serialiser);

        /**
         * Renders this AbstractView3D in the currently active OpenGL context.
         */
        virtual void Render(void);

        /**
         * Resets the view. This normally sets the camera parameters to
         * default values.
         */
        virtual void ResetView(void);

        /**
         * Resizes the AbstractView3D.
         *
         * @param width The new width.
         * @param height The new height.
         */
        virtual void Resize(unsigned int width, unsigned int height);

        /**
         * Sets the button state of a button of the 2d cursor. See
         * 'vislib::graphics::Cursor2D' for additional information.
         *
         * @param button The button.
         * @param down Flag whether the button is pressed, or not.
         */
        virtual void SetCursor2DButtonState(unsigned int btn, bool down) {
            this->cursor2d.SetButtonState(btn, down);
        }

        /**
         * Sets the position of the 2d cursor. See 'vislib::graphics::Cursor2D'
         * for additional information.
         *
         * @param x The x coordinate
         * @param y The y coordinate
         */
        virtual void SetCursor2DPosition(float x, float y) {
            this->cursor2d.SetPosition(x, y, true);
        }

        /**
         * Sets the state of an input modifier.
         *
         * @param mod The input modifier to be set.
         * @param down The new state of the input modifier.
         */
        virtual void SetInputModifier(mmcInputModifier mod, bool down);

        /**
         * Callback requesting a rendering of this view
         *
         * @param call The calling call
         *
         * @return The return value
         */
        virtual bool OnRenderView(Call& call);

        /**
         * Freezes, updates, or unfreezes the view onto the scene (not the
         * rendering, but camera settings, timing, etc).
         *
         * @param freeze true means freeze or update freezed settings,
         *               false means unfreeze
         */
        virtual void UpdateFreeze(bool freeze);

    protected:

        /**
         * Unpacks the mouse coordinates, which are relative to the virtual
         * viewport size.
         *
         * @param x The x coordinate of the mouse position
         * @param y The y coordinate of the mouse position
         */
        virtual void unpackMouseCoordinates(float &x, float &y);

    private:

        /**
         * internal utility class storing frozen values
         */
        class FrozenValues {
        public:

            /**
             * Ctor
             */
            FrozenValues(void) {
                this->camParams = new vislib::graphics::CameraParamsStore();
                this->time = 0.0f;
                this->freezeCounter = 1;
            }

            /** The camera parameters frozen (does not work at all!) */
            vislib::SmartPtr<vislib::graphics::CameraParameters> camParams;

            /** The frame time frozen */
            float time;

            /** The freezeCounter */
            unsigned int freezeCounter;

        };

        /**
         * Implementation of 'Create'.
         *
         * @return 'true' on success, 'false' otherwise.
         */
        virtual bool create(void);

        /**
         * Implementation of 'Release'.
         */
        virtual void release(void);

        /**
         * Renders the vertices of the bounding box
         */
        inline void renderBBox(void);

        /**
         * Renders the back side of the bounding box
         */
        inline void renderBBoxBackside(void);

        /**
         * Renders the front side of the bounding box
         */
        inline void renderBBoxFrontside(void);

        /**
         * Renders the cross for the look-at point
         */
        void renderLookAt(void);

        /**
         * Renders the soft cursor
         */
        void renderSoftCursor(void);

        /**
         * Stores the current camera settings
         *
         * @param p Must be storeCameraSettingsSlot
         *
         * @return true
         */
        bool onStoreCamera(param::ParamSlot& p);

        /**
         * Restores the camera settings
         *
         * @param p Must be restoreCameraSettingsSlot
         *
         * @return true
         */
        bool onRestoreCamera(param::ParamSlot& p);

        /**
         * Resets the view
         *
         * @param p Must be resetViewSlot
         *
         * @return true
         */
        bool onResetView(param::ParamSlot& p);

#ifdef ENABLE_KEYBOARD_VIEW_CONTROL

        /**
         * Event handler for view keys
         *
         * @param p The parameter slot of the view key hit
         *
         * @return true
         */
        bool viewKeyPressed(param::ParamSlot& p);

#endif /* ENABLE_KEYBOARD_VIEW_CONTROL */

        /**
         * Event animation started/stopped
         *
         * @param p Must be animPlaySlot
         *
         * @return true
         */
        bool onAnimPlayChanged(param::ParamSlot& p);

        /**
         * Event animation speed changed
         *
         * @param p Must be animSpeedSlot
         *
         * @return true
         */
        bool onAnimSpeedChanged(param::ParamSlot& p);

        /** The scene camera */
        vislib::graphics::gl::CameraOpenGL cam;

        /** The normal camera parameters */
        vislib::SmartPtr<vislib::graphics::CameraParameters> camParams;

        /** The camera parameter overrides */
        vislib::SmartPtr<vislib::graphics::CameraParameters> camOverrides;

        /** the 2d cursor of this view */
        vislib::graphics::Cursor2D cursor2d;

        /** the input modifiers corresponding to this cursor. */
        vislib::graphics::InputModifiers modkeys;

        /** camera look at rotator */
        vislib::graphics::CameraRotate2DLookAt rotator1;

        /** camera rotator */
        vislib::graphics::CameraRotate2D rotator2;

        /** camera move zoom */
        vislib::graphics::CameraZoom2DMove zoomer1;

        /** camera angle zoom */
        vislib::graphics::CameraZoom2DAngle zoomer2;

        /** camera look-at distance changer */
        vislib::graphics::CameraLookAtDist lookAtDist;

        /** Slot to call the renderer to render */
        CallerSlot rendererSlot;

        /** The light direction vector (NOT LIGHT POSITION) */
        vislib::graphics::SceneSpaceVector3D lightDir;

        /** flag whether or not the light is a camera relative light */
        bool isCamLight;

        /** The complete scene bounding box */
        BoundingBoxes bboxs;

        /** Bool parameter to play/stop the animation */
        param::ParamSlot animPlaySlot;

        /** Float parameter of animation speed in time frames per second */
        param::ParamSlot animSpeedSlot;

        /** The slot holding the current time */
        param::ParamSlot animTimeSlot;

        /** Slot used to synchronize the animation offset */
        param::ParamSlot animOffsetSlot;

        /** The current time to display */
        float timeFrame;

        /** Bool flag showing the bounding box */
        param::ParamSlot showBBox;

        /** Flag showing the look at point */
        param::ParamSlot showLookAt;

        /** The stored camera settings */
        param::ParamSlot cameraSettingsSlot;

        /** Triggers the storage of the camera settings */
        param::ParamSlot storeCameraSettingsSlot;

        /** Triggers the restore of the camera settings */
        param::ParamSlot restoreCameraSettingsSlot;

        /** Triggers the reset of the view */
        param::ParamSlot resetViewSlot;

        /** The frames per second counter */
        vislib::graphics::FpsCounter fpsCounter;

        /** A timer managing the fps output */
        unsigned int fpsOutputTimer;

        /**
         * Flag if this is the first time an image gets created. Used for 
         * initial camera reset
         */
        bool firstImg;

        /** The frozen values */
        FrozenValues *frozenValues;

        /**
         * Flag whether the light is relative to the camera or to the world 
         * coordinate system
         */
        param::ParamSlot isCamLightSlot;

        /** Direction vector of the light */
        param::ParamSlot lightDirSlot;

        /** Diffuse light colour */
        param::ParamSlot lightColDifSlot;

        /** Ambient light colour */
        param::ParamSlot lightColAmbSlot;

        /** focus distance for stereo projection */
        param::ParamSlot stereoFocusDistSlot;

        /** eye distance for stereo projection */
        param::ParamSlot stereoEyeDistSlot;

        /** The diffuse light colour */
        float lightColDif[4];

        /** The ambient light colour */
        float lightColAmb[4];

        /** The incoming call */
        AbstractCallRender *overrideCall;

#ifdef ENABLE_KEYBOARD_VIEW_CONTROL

        /** The move step size in world coordinates */
        param::ParamSlot viewKeyMoveStepSlot;

        /** The angle rotate step in degrees */
        param::ParamSlot viewKeyAngleStepSlot;

        /** The point around which the view will be roateted */
        param::ParamSlot viewKeyRotPointSlot;

        /** Rotates the view to the left (around the up-axis) */
        param::ParamSlot viewKeyRotLeftSlot;

        /** Rotates the view to the right (around the up-axis) */
        param::ParamSlot viewKeyRotRightSlot;

        /** Rotates the view to the top (around the right-axis) */
        param::ParamSlot viewKeyRotUpSlot;

        /** Rotates the view to the bottom (around the right-axis) */
        param::ParamSlot viewKeyRotDownSlot;

        /** Rotates the view counter-clockwise (around the view-axis) */
        param::ParamSlot viewKeyRollLeftSlot;

        /** Rotates the view clockwise (around the view-axis) */
        param::ParamSlot viewKeyRollRightSlot;

        /** Zooms in */
        param::ParamSlot viewKeyZoomInSlot;

        /** Zooms out */
        param::ParamSlot viewKeyZoomOutSlot;

        /** Moves to the left */
        param::ParamSlot viewKeyMoveLeftSlot;

        /** Moves to the right */
        param::ParamSlot viewKeyMoveRightSlot;

        /** Moves to the top */
        param::ParamSlot viewKeyMoveUpSlot;

        /** Moves to the bottom */
        param::ParamSlot viewKeyMoveDownSlot;

#endif /* ENABLE_KEYBOARD_VIEW_CONTROL */

    };


} /* end namespace view */
} /* end namespace core */
} /* end namespace megamol */

#endif /* MEGAMOLCORE_VIEW3D_H_INCLUDED */
