/*
 * Group.h
 *
 * Copyright (C) 2020 by Universitaet Stuttgart (VISUS).
 * Alle Rechte vorbehalten.
 */

#ifndef MEGAMOL_GUI_GRAPH_GROUP_H_INCLUDED
#define MEGAMOL_GUI_GRAPH_GROUP_H_INCLUDED


#include "Module.h"
#include "CallSlot.h"
#include "GUIUtils.h"


namespace megamol {
namespace gui {
namespace configurator {

/**
 * Defines module data structure for graph.
 */
class Group {
public:
  
    Group(ImGuiID uid);
    ~Group();

    const ImGuiID uid;

    // Init when adding group to graph
    std::string name;

    bool AddModule(const ModulePtrType& module_ptr);
    bool DeleteModule(ImGuiID module_uid);

    const ModuleGraphVectorType& GetGroupModules(void) { return this->modules; }

    // GUI Presentation -------------------------------------------------------

    void GUI_Present(StateType& state) {
        this->present.Present(*this, state);
    }

private:

    // VARIABLES --------------------------------------------------------------

    ModuleGraphVectorType modules;

    /**
     * Defines GUI group presentation.
     */
    class Presentation {
    public:
        Presentation(void);

        ~Presentation(void);

        void Present(Group& inout_group, StateType& state);

        void UpdateSize(Group& inout_group, const CanvasType& in_canvas);

    private:
        // Relative position without considering canvas offset and zooming
        ImVec2 position;
        // Relative size without considering zooming
        ImVec2 size;
        GUIUtils utils;
        std::string name_label;
        bool minimized_view;
        bool selected;
        bool update_once;        

    } present;

    // FUNCTIONS --------------------------------------------------------------


};


} // namespace configurator
} // namespace gui
} // namespace megamol

#endif // MEGAMOL_GUI_GRAPH_GROUP_H_INCLUDED