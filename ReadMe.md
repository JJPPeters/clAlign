# clAlign

[GitHub page](http://jjppeters.github.io/clAlign)

A C++ plugin for Gatan Digital Micrograph to align and sum stack images, using OpenCL where possible.
Initial alignments will use a method of solving an overdetermined equation as detailed by [X. Li et al.](http://doi.org/10.1038/NMETH.2472). Small parts of their code have been reused here (made available under the terms of the GNU Public License version 3 (GPLv3))

A lot of the linear algebra part has been taken or reworked from Zhang Ming's work on [SP++](https://code.google.com/p/tspl/) (GNU Public License version 2 (GPLv2))
