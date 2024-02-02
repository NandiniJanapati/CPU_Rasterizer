This is the third assignment I completed for my Computer Graphics class I took at Texas A&M.
The goal was to create a triangle rasterizer using the CPU.

--------
-There are three different .obj files available to draw: sphere.obj, bunny.obj, and duck.obj. The user must type in the file they wish to load at the beginning.
-There are 3 different color modes that can be switched between using '0', '1', and '2'.
    0: Color every triangle in the model with a randomly generated solid color
    1: Color every triangle in the model with 3 randomly generated colors - one color for each vertex - and interpolate between the vertices
    2: Color every pixel with its z-depth buffer value.
-Zoom in / out with 'w' and 's'
-rotate left / right with 'a' and 'd'
-Toggle between using the CPU and GPU with the space bar (the window title will tell the user which they are using) (but also the render time haha)
-If there are texture coordinates available on the obj, 't' will toggle between mapping a texture onto the object or using the selected colormode.
      (since the texture image is a world map, it is recommended to load the sphere.obj to view the texture mapping)
-If the texture mapping is toggled on, 'n', 'm', and 'l' changes how the texture is mapped
    'n': Map texture using nearest neighbor
    'l': Map texture using bilinear interpolation
    'm': Map texture using mip-mapping
-'q' to quit
