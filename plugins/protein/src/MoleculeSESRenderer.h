/*
 * MoleculeSESRenderer.h
 *
 * Copyright (C) 2009-2021 by Universitaet Stuttgart (VIS). Alle Rechte vorbehalten.
 */

#ifndef MMPROTEINPLUGIN_MOLSESRENDERER_H_INCLUDED
#define MMPROTEINPLUGIN_MOLSESRENDERER_H_INCLUDED
#if (_MSC_VER > 1000)
#pragma once
#endif /* (_MSC_VER > 1000) */

#include <algorithm>
#include <list>
#include <set>
#include <vector>
#include "Color.h"
#include "Pyramid.h"
#include "ReducedSurface.h"
#include "mmcore/CallerSlot.h"
#include "mmcore/param/ParamSlot.h"
#include "mmcore/view/CallRender3DGL.h"
#include "mmcore/view/Camera_2.h"
#include "mmcore/view/Renderer3DModuleGL.h"
#include "protein_calls/BindingSiteCall.h"
#include "protein_calls/MolecularDataCall.h"
#include "stb_image.h"
#include "vislib/Array.h"
#include "vislib/String.h"
#include "vislib/graphics/gl/GLSLComputeShader.h"
#include "vislib/graphics/gl/GLSLGeometryShader.h"
#include "vislib/graphics/gl/GLSLShader.h"
#include "vislib/graphics/gl/SimpleFont.h"
#include "vislib/math/Quaternion.h"

namespace megamol {
namespace protein {

    /**
     * Molecular Surface Renderer class.
     * Computes and renders the solvent excluded (Connolly) surface.
     */
    class MoleculeSESRenderer : public megamol::core::view::Renderer3DModuleGL {
    public:
        /**
         * Answer the name of this module.
         *
         * @return The name of this module.
         */
        static const char* ClassName(void) {
            return "MoleculeSESRenderer";
        }

        /**
         * Answer a human readable description of this module.
         *
         * @return A human readable description of this module.
         */
        static const char* Description(void) {
            return "Offers protein surface renderings.";
        }

        /**
         * Answers whether this module is available on the current system.
         *
         * @return 'true' if the module is available, 'false' otherwise.
         */
        static bool IsAvailable(void) {
            // return true;
            return vislib::graphics::gl::GLSLShader::AreExtensionsAvailable();
        }

        /** ctor */
        MoleculeSESRenderer(void);

        /** dtor */
        virtual ~MoleculeSESRenderer(void);

        /**********************************************************************
         * 'get'-functions
         **********************************************************************/

        /** Get probe radius */
        const float GetProbeRadius() const {
            return probeRadius;
        };

        /**********************************************************************
         * 'set'-functions
         **********************************************************************/

        /** Set probe radius */
        void SetProbeRadius(const float rad) {
            probeRadius = rad;
        };

    protected:
        /**
         * Implementation of 'Create'.
         *
         * @return 'true' on success, 'false' otherwise.
         */
        virtual bool create(void);

        /**
         * Implementation of 'release'.
         */
        virtual void release(void);

        /**
         * Renders the probe atom at position 'm'.
         *
         * @param m The probe position.
         */
        // void RenderProbe(const vislib::math::Vector<float, 3> m);
        void RenderProbeGPU(const vislib::math::Vector<float, 3> m);

        /**
         * Compute all vertex, attribute and color arrays used for ray casting
         * all molecular surfaces (spheres, spherical triangles, tori).
         */
        void ComputeRaycastingArrays();

        /**
         * Compute all vertex, attribute and color arrays used for ray casting
         * the molecular surface 'ixdRS' (spheres, spherical triangles, tori).
         * @param idxRS The index of the reduced surface.
         */
        void ComputeRaycastingArrays(unsigned int idxRS);

        /**
         * Code a RGB-color into one float.
         * For each color channel, its representation in range 0..255 is computed
         * and stores as follows:
         * rrrgggbbb.0
         * Note that the minimum value for the coded color is 0 and the maximum
         * value is 255255255.0 .
         *
         * @param col Vector containing the color as float [0.0]..[1.0] .
         * @return The coded color value.
         */
        // float CodeColor( const vislib::math::Vector<float, 3> &col) const;
        float CodeColor(const float* col) const;

        /**
         * Decode a coded color to the original RGB-color.
         *
         * @param codedColor Integer value containing the coded color (rrrgggbbb).
         * @return The RGB-color value vector.
         */
        vislib::math::Vector<float, 3> DecodeColor(int codedColor) const;

        /**
         * Creates the frame buffer object and textures needed for offscreen rendering.
         */
        void CreateFBO();

        /**
         * Creates the frame buffer object and textures needed for offscreen rendering.
         */
        void CreateQuadBuffers();

        /**
         * Render the molecular surface using GPU raycasting.
         *
         * @param protein Pointer to the protein data interface.
         */
        void RenderSESGpuRaycasting(const megamol::protein_calls::MolecularDataCall* mol);

        /**
         * Postprocessing: Calculate Suggestive Contours using minima of diffuse shading
         */
        void SuggestiveContours(vislib::graphics::gl::GLSLShader& Shader);
        void Contours(vislib::graphics::gl::GLSLShader& Shader);

        void displayPositions();
        void displayNormalizedPositions();
        void displayNormals();

        /**
         * Postprocessing: Calculate Suggestive Contours from curvature information
         */
        void calculateCurvature(vislib::graphics::gl::GLSLShader& Shader);
        void renderCurvature(vislib::graphics::gl::GLSLShader& Shader);
        void calculateTextureBBX();
        // void calculateRMSF();

        /**
         * Smooth normal texture using pull-push algorithm
         */
        void SmoothNormals(vislib::graphics::gl::GLSLShader& Shader);
        void SmoothCurvature(vislib::graphics::gl::GLSLShader& Shader);
        void SmoothPositions(vislib::graphics::gl::GLSLShader& Shader);

        /**
         * returns the color of the atom 'idx' for the current coloring mode
         *
         * @param idx The index of the atom.
         * @return The color of the atom with the index 'idx'.
         */
        vislib::math::Vector<float, 3> GetProteinAtomColor(unsigned int idx);

        /**
         * Create the singularity textureS which stores for every RS-edge (of all
         * molecular surfaces) the positions of the probes that cut it.
         */
        void CreateSingularityTextures();

        /**
         * Create the singularity texture for the reduced surface 'idxRS' which
         * stores for every RS-edge the positions of the probes that cut it.
         */
        void CreateSingularityTexture(unsigned int idxRS);

    private:
        /**
         * Update all parameter slots.
         *
         * @param mol   Pointer to the data call.
         */
        void UpdateParameters(
            const megamol::protein_calls::MolecularDataCall* mol, const protein_calls::BindingSiteCall* bs = 0);

        bool loadShader(vislib::graphics::gl::GLSLShader& Shader, std::string vertex, std::string fragment);

        /**
         * The get extents callback. The module should set the members of
         * 'call' to tell the caller the extents of its data (bounding boxes
         * and times).
         *
         * @param call The calling call.
         *
         * @return The return value of the function.
         */
        virtual bool GetExtents(megamol::core::view::CallRender3DGL& call);

        /**
         * Open GL Render call.
         *
         * @param call The calling call.
         * @return The return value of the function.
         */
        virtual bool Render(megamol::core::view::CallRender3DGL& call);

        /**
         * Deinitialises this renderer. This is only called if there was a
         * successful call to "initialise" before.
         */
        virtual void deinitialise(void);

        /**********************************************************************
         * variables
         **********************************************************************/

        /** MolecularDataCall caller slot */
        megamol::core::CallerSlot molDataCallerSlot;
        /** BindingSiteCall caller slot */
        megamol::core::CallerSlot bsDataCallerSlot;

        /** camera information */
        // vislib::SmartPtr<vislib::graphics::CameraParameters> cameraInfo;
        core::view::Camera_2 cameraInfo;

        // camera information
        // vislib::SmartPtr<vislib::graphics::CameraParameters> MoleculeSESRenderercameraInfo;
        core::view::Camera_2 MoleculeSESRenderercameraInfo;
        /** parameter slot for coloring mode */
        megamol::core::param::ParamSlot coloringModeParam0;
        /** parameter slot for coloring mode */
        megamol::core::param::ParamSlot coloringModeParam1;
        megamol::core::param::ParamSlot curvatureModeParam;
        megamol::core::param::ParamSlot contourModeParam;
        megamol::core::param::ParamSlot displayedPropertyParam;
        /** parameter slot for coloring mode weighting*/
        megamol::core::param::ParamSlot cmWeightParam;
        /** parameter slot for min color of gradient color mode */
        megamol::core::param::ParamSlot minGradColorParam;
        /** parameter slot for mid color of gradient color mode */
        megamol::core::param::ParamSlot midGradColorParam;
        /** parameter slot for max color of gradient color mode */
        megamol::core::param::ParamSlot maxGradColorParam;
        megamol::core::param::ParamSlot drawSESParam;
        megamol::core::param::ParamSlot drawSASParam;
        megamol::core::param::ParamSlot molIdxListParam;
        /** parameter slot for color table filename */
        megamol::core::param::ParamSlot colorTableFileParam;
        /** Parameter to toggle offscreen rendering */
        megamol::core::param::ParamSlot offscreenRenderingParam;
        megamol::core::param::ParamSlot probeRadiusSlot;
        /** Pyramid Parameters */
        megamol::core::param::ParamSlot pyramidWeightsParam;
        GLfloat pyramidWeight;
        megamol::core::param::ParamSlot smoothNormalsParam;
        GLboolean smoothNormals;
        megamol::core::param::ParamSlot smoothCurvatureParam;
        GLboolean smoothCurvature;
        megamol::core::param::ParamSlot smoothPositionsParam;
        GLboolean smoothPositions;
        megamol::core::param::ParamSlot pyramidLayersParam;
        GLint pyramidLayers;
        megamol::core::param::ParamSlot pyramidGammaParam;
        GLfloat pyramidGamma;
        /** Suggestive Contours Parametes */
        megamol::core::param::ParamSlot SCRadiusParam;
        GLint SCRadius;
        megamol::core::param::ParamSlot SCNeighbourThresholdParam;
        GLfloat SCNeighbourThreshold;
        megamol::core::param::ParamSlot SCDiffThresholdParam;
        GLfloat SCDiffThreshold;
        megamol::core::param::ParamSlot SCMedianFilterParam;
        GLboolean SCMedianFilter;
        megamol::core::param::ParamSlot SCCircularNeighborhoodParam;
        GLboolean SCCircularNeighborhood;
        megamol::core::param::ParamSlot OrthogonalViewParam;
        GLboolean orthogonalView;
        megamol::core::param::ParamSlot OrthoProjParam;
        GLboolean orthoproj;
        megamol::core::param::ParamSlot TestCaseParam;
        GLboolean testcase;
        megamol::core::param::ParamSlot cutOffParam;
        GLfloat cutOff;
        megamol::core::param::ParamSlot nearPlaneParam;
        GLfloat nearplane;
        megamol::core::param::ParamSlot blurParam;
        megamol::core::param::ParamSlot numBlurParam;
        GLint numBlur;
        megamol::core::param::ParamSlot numPosBlurParam;
        GLint numPosBlur;
        megamol::core::param::ParamSlot numCurvBlurParam;
        GLint numCurvBlur;
        megamol::core::param::ParamSlot depthDiffParam;
        GLfloat depthDiff;

        bool drawSES;
        bool drawSAS;
        bool offscreenRendering;

        // Testing
        unsigned int texturePy;

        /** the reduced surface(s) */
        std::vector<std::vector<ReducedSurface*>> reducedSurfaceAllFrames;
        /** the reduced surface(s) */
        std::vector<ReducedSurface*> reducedSurface;

        // The pull-push pyramids
        Pyramid normalPyramid;    // Normal smoothing
        Pyramid positionPyramid;  // Position smoothing
        Pyramid SCpyramid;        // Suggestive Contours using pull-push algorithm
        Pyramid curvaturePyramid; // curvature smoothing

        // BBX pyramids
        Pyramid depthPyramid;  // Calculation of maximum view-space depth
        Pyramid heightPyramid; // Calculation of maximum view-space height
        Pyramid widthPyramid;  // Calculation of maximum view-space width

        int bbx_levelMax; // How many layers do the bbx pyramids have?

        glm::vec4 clear_color;

        // shader for the spheres (raycasting view)
        vislib::graphics::gl::GLSLShader sphereShader;
        vislib::graphics::gl::GLSLShader sphereShaderOR;
        // shader for the spherical triangles (raycasting view)
        vislib::graphics::gl::GLSLShader sphericalTriangleShader;
        vislib::graphics::gl::GLSLShader sphericalTriangleShaderOR;
        // shader for torus (raycasting view)
        vislib::graphics::gl::GLSLShader torusShader;
        vislib::graphics::gl::GLSLShader torusShaderOR;
        // shader for per pixel lighting (polygonal view)
        vislib::graphics::gl::GLSLShader lightShader;
        // shader for contour generation
        vislib::graphics::gl::GLSLShader SC_Shader;
        vislib::graphics::gl::GLSLShader SC_Curvature_Shader;
        vislib::graphics::gl::GLSLShader C_Curvature_Shader;
        vislib::graphics::gl::GLSLShader C_Shader;
        // shader for curvature calculation
        vislib::graphics::gl::GLSLShader curvatureShader;
        vislib::graphics::gl::GLSLShader normalCurvatureShader;
        vislib::graphics::gl::GLSLShader nathanReedCurvatureShader;
        vislib::graphics::gl::GLSLShader meanCurvatureShader;
        vislib::graphics::gl::GLSLShader prantlMeanShader;
        vislib::graphics::gl::GLSLShader prantl2MeanShader;
        vislib::graphics::gl::GLSLShader prantlGaussianShader;
        vislib::graphics::gl::GLSLShader prantl2GaussianShader;
        vislib::graphics::gl::GLSLShader prantlRadialShader;
        vislib::graphics::gl::GLSLShader prantl2RadialShader;
        // pass through Shader sampling from a texture at mipmap level 0
        vislib::graphics::gl::GLSLShader passThroughShader;
        vislib::graphics::gl::GLSLShader normalizePositionsShader;
        vislib::graphics::gl::GLSLShader gaussianBlurShader;
        vislib::graphics::gl::GLSLShader peronaMalikBlurShader;
        vislib::graphics::gl::GLSLShader depthBlurShader;
        vislib::graphics::gl::GLSLShader depthGaussBlurShader;
        ////////////

        // the bounding box of the protein
        vislib::math::Cuboid<float> bBox;

        // epsilon value for float-comparison
        float epsilon;

        // radius of the probe atom
        float probeRadius;

        vislib::Array<float> atomColorTable;
        unsigned int currentArray;

        /** 'true' if the data for the current render mode is computed, 'false' otherwise */
        bool preComputationDone;

        /** The current coloring mode */
        Color::ColoringMode currentColoringMode0;
        Color::ColoringMode currentColoringMode1;

        enum curvatureMode {
            EvansCurvature,
            NormalCurvature,
            NathanReedCurvature,
            MeanCurvature,
            PrantlMean,
            Prantl2Mean,
            PrantlGaussian,
            Prantl2Gaussian,
            PrantlRadial,
            Prantl2Radial
        };
        curvatureMode currentCurvatureMode;
        std::map<curvatureMode, vislib::graphics::gl::GLSLShader*> curvatureShaderMap = {
            {EvansCurvature, &this->curvatureShader}, {NormalCurvature, &this->normalCurvatureShader},
            {NathanReedCurvature, &this->nathanReedCurvatureShader}, {MeanCurvature, &this->meanCurvatureShader},
            {PrantlMean, &this->prantlMeanShader}, {Prantl2Mean, &this->prantl2MeanShader},
            {PrantlGaussian, &this->prantlGaussianShader}, {Prantl2Gaussian, &this->prantl2GaussianShader},
            {PrantlRadial, &this->prantlRadialShader}, {Prantl2Radial, &this->prantl2RadialShader}};
        enum contourMode { Suggestive, SuggestiveAndCurvature, Shading, ShadingAndCurvature };
        std::map<contourMode, vislib::graphics::gl::GLSLShader*> contourShaderMap = {{Suggestive, &this->SC_Shader},
            {SuggestiveAndCurvature, &this->SC_Curvature_Shader}, {Shading, &this->C_Shader},
            {ShadingAndCurvature, &this->C_Curvature_Shader}};
        contourMode currentContourMode;
        enum blurMode { Gaussian, PeronaMalik, DepthSensitive, DepthGaussian };
        std::map<blurMode, vislib::graphics::gl::GLSLShader*> blurShaderMap = {{Gaussian, &this->gaussianBlurShader},
            {PeronaMalik, &this->peronaMalikBlurShader}, {DepthSensitive, &this->depthBlurShader},
            {DepthGaussian, &this->depthGaussBlurShader}};
        blurMode currentBlurMode;
        enum displayedProperty { Position, NormalizedPosition, Normal, Curvature, Contour };
        displayedProperty currentDisplayedProperty;

        /** vertex and attribute arrays for raycasting the tori */
        std::vector<vislib::Array<float>> torusVertexArray;
        std::vector<vislib::Array<float>> torusInParamArray;
        std::vector<vislib::Array<float>> torusQuatCArray;
        std::vector<vislib::Array<float>> torusInSphereArray;
        std::vector<vislib::Array<float>> torusColors;
        std::vector<vislib::Array<float>> torusInCuttingPlaneArray;
        /** vertex ans attribute arrays for raycasting the spherical triangles */
        std::vector<vislib::Array<float>> sphericTriaVertexArray;
        std::vector<vislib::Array<float>> sphericTriaVec1;
        std::vector<vislib::Array<float>> sphericTriaVec2;
        std::vector<vislib::Array<float>> sphericTriaVec3;
        std::vector<vislib::Array<float>> sphericTriaProbe1;
        std::vector<vislib::Array<float>> sphericTriaProbe2;
        std::vector<vislib::Array<float>> sphericTriaProbe3;
        std::vector<vislib::Array<float>> sphericTriaTexCoord1;
        std::vector<vislib::Array<float>> sphericTriaTexCoord2;
        std::vector<vislib::Array<float>> sphericTriaTexCoord3;
        std::vector<vislib::Array<float>> sphericTriaColors;
        /** vertex and color array for raycasting the spheres */
        std::vector<vislib::Array<float>> sphereVertexArray;
        std::vector<vislib::Array<float>> sphereColors;

        // FBOs and textures for postprocessing
        GLuint contourFBO;
        GLuint curvatureFBO;
        GLuint positionFBO[2];
        GLuint normalFBO[2];
        GLuint smoothCurvFBO[2];
        GLuint normalTexture;
        GLuint smoothNormalTexture[2];
        // GLuint* cur_normalTexture;
        // GLuint* cur_positionTexture;
        GLuint curvatureTexture;
        GLuint smoothCurvatureTexture[2];
        GLuint positionTexture;
        GLuint smoothPositionTexture[2];
        GLuint objPositionTexture;
        GLuint contourDepthRBO;

        // boolean for pingpong blurring
        bool horizontal;

        // VAO and VBO for screen filling quad
        unsigned int quadVAO, quadVBO;

        // width and height of view
        unsigned int width;
        unsigned int height;

        /** The color lookup table (for chains, amino acids,...) */
        vislib::Array<vislib::math::Vector<float, 3>> colorLookupTable;
        /** The color lookup table which stores the rainbow colors */
        vislib::Array<vislib::math::Vector<float, 3>> rainbowColors;

        // texture for singularity handling (concave triangles)
        std::vector<GLuint> singularityTexture;
        // sizes of singularity textures
        std::vector<unsigned int> singTexWidth, singTexHeight;
        // data of the singularity texture
        float* singTexData;

        // the list of molecular indices
        vislib::Array<vislib::StringA> molIdxList;
        // flag for SES computation (false = one SES per molecule)
        bool computeSesPerMolecule;
        void RenderTestCase();
        vislib::graphics::gl::GLSLGeometryShader testCaseShader;
    };

} /* end namespace protein */
} /* end namespace megamol */

#endif /* MMPROTEINPLUGIN_MOLSESRENDERER_H_INCLUDED */
