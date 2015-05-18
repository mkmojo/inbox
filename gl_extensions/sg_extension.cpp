#include <graphlab/sdk/toolkit_function_macros.hpp>
#include <graphlab/sdk/gl_sgraph.hpp>
using namespace graphlab;

gl_sgraph construct_graph(){
    gl_sgraph g;
    gl_sframe vertices { {"vid", {1,2,3}} };
    gl_sframe edges {{"src", {1, 3}}, {"dst", {2, 2}}};
    return g = g.add_vertices(vertices, "vid").add_edges(edges, "src", "dst");
}

BEGIN_FUNCTION_REGISTRATION
REGISTER_FUNCTION(construct_graph);
END_FUNCTION_REGISTRATION
