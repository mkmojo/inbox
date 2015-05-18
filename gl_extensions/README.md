1. FLEX_UNDEFINED is the macro that will later be used as None in Python.
2. flexible_type is the magical type that can take many python structures.
3. When you want to publish functions, can use macro to do so.

sgraph example
====

When trying to write extensions to GraphLab, need to wrap program into 
functions. 

For instance cannot have this directly shows up in the source,
```C++
#innclude <graphlab/sdk/toolkit_function_macros.hpp>
#include <graphlab/sdk/gl_sgraph.hpp>

gl_sgraph g;  // empty graph
gl_sframe vertices { {"vid", {1,2,3} } };
gl_sframe edges { {"src", {1, 3}}, {"dst", {2, 2}} };
g = g.add_vertices(vertices, "vid").add_edges(edges, "src", "dst");
```
Instead, we should give it a *NAME*.

```C++
gl_sgraph& construct_graph(){
    gl_sgraph g;
    gl_sframe vertices { {"vid", {1,2,3}} };
    gl_sframe edges {{"src", {1, 3}}, {"dst", {2, 2}}};
    return g = g.add_vertices(vertices, "vid").add_edges(edges, "src", "dst");
}
```
