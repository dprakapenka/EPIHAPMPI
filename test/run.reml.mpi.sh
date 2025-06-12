mpirun -n 2 ../bin/reml \
    -p 7,../example/example.phen \
    -i ./saved/grm.mpi \
    --va 3 \
    --vd 1 \
    --ve 1 \
    -o out/mpi \
    -I 10
#mpirun -n 2 ../bin/reml -p 7,../example/example.phen -i ./saved/grm.mpi --map ../example/map.txt -a 8.237762e+00 -d 2 -e 1.462395e+01 -o out -I 3
