// SPDX-License-Identifier: LGPL-3.0-or-later
//
// This file is part of MOO - Modelica / Model Optimizer
// Copyright (C) 2025 University of Applied Sciences and Arts
// Bielefeld, Faculty of Engineering and Mathematics
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//

#ifndef MOO_TIMING_H
#define MOO_TIMING_H

#include <base/util.h>
#include <base/log.h>

#include <string>
#include <vector>
#include <memory>
#include <chrono>
#include <unordered_map>

class TimingNode;

using visitor_fn = std::function<void(const TimingNode* node, const TimingNode* parent, const TimingNode* next)>;

class MOO_EXPORT TimingNode {
    friend class TimingTree;

public:
    std::string name;
    f64 duration_nano = 0.0;
    std::vector<std::unique_ptr<TimingNode>> children;
    TimingNode* parent = nullptr;

    TimingNode(std::string n, TimingNode* p = nullptr);
    virtual ~TimingNode() = default;

    TimingNode* add_child(std::unique_ptr<TimingNode> child);

    // called on destruction
    // use it to set #steps, #failures, etc for the node
    // use it to insert VirtualTimer's (e.g. for 3rdParty libs)
    virtual void finalize();

private:
    void print_table_tree(const FixedTableFormat<4>& tf, const std::string& prefix = "", bool is_last  = true) const;
};

class MOO_EXPORT TimingTree {
public:
    static TimingTree& instance();
    TimingTree(const TimingTree&) = delete;
    TimingTree& operator=(const TimingTree&) = delete;

    TimingNode* current_node();
    void set_current(TimingNode* n);
    TimingNode* add_child(std::unique_ptr<TimingNode> child);

    void print_tree_table() const;
    void traverse(visitor_fn visitor) const;

private:
    TimingTree();

    std::unique_ptr<TimingNode> root;
    TimingNode* current;

    void traverse_node(const TimingNode* node,
                       const TimingNode* parent,
                       const TimingNode* next,
                       visitor_fn& visitor) const;
};

template <typename NodeType = TimingNode, typename... Args>
class MOO_EXPORT ScopedTimer {
public:
    ScopedTimer(const std::string& name, Args&&... args)
        : start(std::chrono::high_resolution_clock::now())
    {
        auto& tree = TimingTree::instance();
        auto parent = tree.current_node();
        node = parent->add_child(std::make_unique<NodeType>(name, parent, std::forward<Args>(args)...));
        tree.set_current(node);
    }

    ~ScopedTimer() {
        f64 nano = std::chrono::duration<f64, std::nano>(std::chrono::high_resolution_clock::now() - start).count();
        node->duration_nano = nano;
        node->finalize();

        TimingTree::instance().set_current(node->parent);
    }

private:
    TimingNode* node;
    std::chrono::high_resolution_clock::time_point start;
};

template <typename NodeType = TimingNode, typename... Args>
class MOO_EXPORT VirtualTimer {
public:
    VirtualTimer(const std::string& name, f64 duration_nano, Args&&... args)
    {
        auto& tree = TimingTree::instance();
        auto parent = tree.current_node();
        node = parent->add_child(std::make_unique<NodeType>(name, parent, std::forward<Args>(args)...));
        node->duration_nano = duration_nano;
        tree.set_current(node);
    }

    ~VirtualTimer() {
        node->finalize();
        TimingTree::instance().set_current(node->parent);
    }

private:
    TimingNode* node;
};

class MOO_EXPORT CountedTimingNode : public TimingNode {
public:
    int cnt = 0;

    CountedTimingNode(std::string n, TimingNode* p, int cnt);
};

namespace Timing {

constexpr f64 s_to_nano(f64 seconds) { return seconds * 1e9; }
constexpr f64 nano_to_ms(f64 seconds) { return seconds * 1e-6; }

// TODO: remove strings and make these type-safe - via enum or so?
std::unordered_map<std::string, f64> accumulate_prefix(std::string prefix);
std::vector<f64> accumulate_blocks(const std::string& prefix_start,
                                   const std::string& prefix_break);

} // namespace Timing

#endif // MOO_TIMING_H
