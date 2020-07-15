============
Introduction
============

---------------------------------------------------
Tools for image forming optical design and analysis
--------------------------------------------------- 

The **ray-optics** project has several goals, both in the optical domain and in the software development domain  

* Rethink how image forming optical calculations are done absent historical constraints on computer speed and memory
* Investigate interactive graphics for optical design
* Serve as a reference implementation of basic ray tracing and wavefront analysis algorithms
* Leverage extensive Python libraries to avoid reinventing the wheel
* Look for architectures that can be deployed on the desktop (e.g. using Qt), using (Jupyter) notebooks, and in mobile environments

Image forming optical design and analysis was one of the first domain areas to be transferred from the realm of human computers to the first electronic computing machines. The calculations were so voluminous and tedious that automation was a tremendous benefit. Because the capabilities of electronic computers were extremely limited, the computational workflow was optimized to be maximally efficent for computational effort, paralleling the practice with human computers. 

Computers are vastly more powerful now than they were when the venerable CODE V program was initially developed. Zemax and Oslo date from the early IBM PC days, when both speed and memory were limiting factors. The programs were developed to be comprehensive tools. The software tools available then were limited as well. In order to gain acceptable performance, compiled languages such as C and FORTRAN were required. Graphical user interfaces were also expensive to develop and were often considered secondary in importance to developing a comprehensive feature set.

The **ray-optics** project doesn't aspire to be comprehensive. Instead, it tries to do several dominant use cases very well:

* Rotationally symmetric optical systems
* Bilaterally symmetric optical systems

These two categories represent a vast space of practical imaging and/or coherent optical systems. The reason to restrict the program in this way is simplify the development of an interactive graphic user interface focused on optical modeling. These use cases can be represented in 2D with no loss in generality. This offers the potential of significently simplified and/or more focused tools for manipulation of optical properties.

The rotationally symmetric case has a rich set of paraxial or first order design methods, in addition to a well developed aberration theory. The interactive optical layout view uses the ideas behind paraxial ray solves to implement graphical editing operations directly with paraxial rays. Interactive |ybar| and |nubar| diagrams provide alternative forms of paraxial design.

In an interactive design environment, the graphical renderings mostly correspond with optical elements, not individual surfaces (excepting mirrors). The virtue of having an element-based description of the optical model, in parallel with a sequential interface description, is highly leveraged when graphical tools become a focus of the software. It is important that either the sequence or the element description be able to generate the other description. 

The core analysis operation is geometrical ray tracing through a sequence of optical interfaces.

Optical calculation technology can be considered a mature field. There is a long history in the literature of investigations tying optical theory to calculation approaches that maintain accuracy, handle corner cases, etc. Much of this knowledge is embedded in production optical design software; having it available in an open source implementation brings this knowledge out of the literature and makes it accessible to students and researchers.

The advent of scripting environments like Python make a fresh look at optical design calculations worthwhile. Python has many advantages for scientific and engineering applications, including libraries for numerical calculations, a variety of graphics options including a very full featured plotting library, and support for data science applications that look promising for managing data generated during extensive analysis of optical systems. There is also good support for unit testing of Python software, as well as debugging and performance profiling.

Finally, computational notebooks, such as those of the Jupyter project, provide the ability to document **reproducibly** the implementation of an algorithm and results produced by the algorithm.
