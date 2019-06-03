/** @file mainpage.h
* @brief Quadrics visualization and texturing
*
* $Header: /mainpage.html,v 1.1.1.1 2017/12/12 15:03:16 chuckp Exp $
*/
/** @mainpage Quadrics rendering
*
* @author Paulo Roma Cavalcanti
*
* @section Goals
* This program provides code templates for use by OpenGL developers. 
* Although simple, the code supplies everything a beginner in Computer Graphics
* should know:
* 
* - Basic <A HREF="../../OpenGL/surfaces.pdf">math</A>
* -# Surface <A HREF="http://www.hao-li.com/cs599-ss2015/slides/Lecture02.1.pdf">Implicit representation</A>
* -# Surface <A HREF="http://tutorial.math.lamar.edu/Classes/CalcIII/ParametricSurfaces.aspx">Parametric representation</A>
* -# <A HREF="https://www.khanacademy.org/math/linear-algebra/vectors-and-spaces/dot-cross-products/v/normal-vector-from-plane-equation">Normal vector</A> calculation
* -# <A HREF="https://en.wikipedia.org/wiki/Quadric">Quadrics</A>
* -# <A HREF="https://blogs.scientificamerican.com/roots-of-unity/a-few-of-my-favorite-spaces-the-torus/">Torus</A> (ring, horn or spindle)
* - Several render styles:
* -# <A HREF="https://knowledge.autodesk.com/support/autocad/getting-started/caas/CloudHelp/cloudhelp/2016/ENU/AutoCAD-Core/files/GUID-84E193D7-A18D-4EE2-B978-19E4AFBCAEEC-htm.html">Wireframe</A>
* -# <A HREF="http://graphics.wikia.com/wiki/Gouraud_shading">Smooth shading</A>
* -# <A HREF="https://www.textures.com">Texture</A> (procedural or image based)
* - <A HREF="http://mathworld.wolfram.com/NormalVector.html">Surface normal</A> rendering
* -# drawn proportionally to the bounding box diagonal
* -# pointing outside the object
* - <A HREF="https://en.wikipedia.org/wiki/Minimum_bounding_box">Bounding Box</A> computation
* - <A HREF="https://msdn.microsoft.com/en-us/library/hh920765(v=vs.85).aspx">Timer</A> for controlling basic animation
*
*
* Also note the existence of the following directories:
* - Arcball
* -# Rotation <A HREF="../../../WebGL/extras/doc/Arcball.pdf">paradigm</A> for CG.
* - doc
* -# <A HREF="torus_8cpp.html">Torus.cpp</A>
*
*
* As you prepare to develop code for the course, please be sure you are aware 
* of our current
* <A HREF="https://www.gnu.org/prep/standards/standards.html"> Coding Standards </A>
*
*
* If using the code in this package as an example - please modify the comments
* as appropriate for your own specific code.
*
* <hr>
* @section notes release.notes
* This program runs on linux 64 or 32 bits (Fedora, Ubuntu), and MacOS.
* To generate executables, just run:
* - for linux:
* -# make
* - for MacOS:
* -# make OS=osx
* - for cross-compiling for windows:
* -# make OS=cross
*
* <hr>
* @section requirements 
*
* - C++ compiler
* - <A HREF="https://www.talisman.org/opengl-1.1/Reference.html">OpenGL 1.1</A>
* - <A HREF="http://openil.sourceforge.net/">DevIL</A>
* - <A HREF="http://eigen.tuxfamily.org/index.php?title=Main_Page">Eigen3</A>
*
* @verbinclude requirements
* <hr> 
* @todo [graphical interface, Java Implementation]
* @todo Qt4 or <A HREF="https://www.qt.io">Qt5</A> based
* @todo <A HREF="http://www.cs.cornell.edu/courses/cs4620/2011fa/lectures/practicum01.pdf">Introduction to OpenGL</A>, <A HREF="http://jogamp.org/jogl/www/">JOGL</A>
*
*
* <hr>
* @section pictures Screen Shots
*
* @imageSize{torus.png,height:272px;width:267px;}
* @imageSize{torus.same-ratio.png,height:272px;width:267px;}
* @imageSize{torus.tube.gt.circle.png,height:272px;width:267px;}
* \htmlonly <style>div.image img[src="paraboloid_hiperbolic.png"]{width:545px;}</style> \endhtmlonly 
* @image html torus.png "Ring Torus" 
* @image html torus.same-ratio.png "Horn Torus" 
* @image html torus.tube.gt.circle.png "Spindle Torus" 
* @image html paraboloid_hiperbolic.png "Hyperbolic Paraboloid"
* \htmlonly
* <p style="float: left; font-size: 12pt; text-align: center; width: 30%; margin-right: 1%; margin-bottom: 0.5em;"><img src="torus.png" style="width: 100%">Ring Torus</p>
* <p style="float: left; font-size: 12pt; text-align: center; width: 30%; margin-right: 1%; margin-bottom: 0.5em;"><img src="torus.same-ratio.png" style="width: 100%">Horn Torus</p>
* <p style="float: left; font-size: 12pt; text-align: center; width: 30%; margin-right: 1%; margin-bottom: 0.5em;"><img src="torus.tube.gt.circle.png" style="width: 100%">Spindle Torus</p>
* <p style="clear: both;">
* \endhtmlonly
*
*/
