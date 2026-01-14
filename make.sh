#run it from src folder; Need to copy the .cpp to KSPACE folder
cp pair_vspbks.* KSPACE/
make
make yes-KSPACE
make mpi
