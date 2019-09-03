The one-dimensional Godunov scheme of the first order for the equations of fluid dynamics
has been implemented via using OpenMP technology for multi-threading.

The two versions of the scheme are enabled: with nonlinear (1959) and linear (2018) solution of the Riemann problem on the boundaries.
The user can switch them in the code by itself (see definitions.h)

Linux:
For better performance, please compile the code with Intel Compiler:
icc -ipo -O2 -qopenmp -o a.out main.cpp functions.cpp iteration.cpp gnuplot.cpp arrays.cpp

Run the appilcation as
./a.out

Windows:
Just create the project in Visual Studio IDE with the source folder, compile and run it.
Do not forget to set option /Qopenmp in the settings of the project to run OpenMP version.
For better performance, turn on keys /Qipo /O2 /QxHost in the settings too.

Before running, set the number (n) of OpenMP threads with command:
set OMP_NUM_THREADS=n (for Windows)
export OMP_NUM_THREADS=n (for Linux)

This application generates output files and print the results with GNUPLOT application (if respective #define was turn on)
ATTENTION! Output directory must have a special structure of folders. Please, see gnuplot.cpp for the further information.

Please, ask your questions at dmitriy_klyuchinskiy@mail.ru
