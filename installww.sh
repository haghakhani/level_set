make
make install
#export LD_LIBRARY_PATH="/home/hossein/petsc/arch-linux2-c/lib"
cd bin/cc/
cp ../titan .
./titan
#mpirun -np 2 ./titan
