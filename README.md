# LoopClosure-Detector
C++ implementation of a LoopClosure Detector for a 2D laser scans sensor.


Libraries needed for compiling:
- Eigen3
- OpenCV

loop-detector-2d.txt contains the output of a laser scanner.

trajectoryLaser.g2o contains some key nodes out from the total poses where scans have been made, together with their pose (extracted from odometry).

To compile:

mkdir build  
cd build  
cmake ..  
make   

To execute:

./display_test This executable allows to display all the scans: use "A" and "D" keys to browse between them.


./alignment_test This executable allows to run optimization between a pair of scans. The 2 scans can be given as input arguments (e.g. ./alignment_test 32975 31308).
You can run 1 optimization round using "spacebar", the full optimization until stopping criterion using "enter", display the current best transformation using "X".


./closure_finder This will automatically produce a list of candidate closures, written into a .txt file. 

