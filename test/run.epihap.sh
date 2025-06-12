#!/bin/bash
module load intel-oneapi-mkl

/project/dairyxbreed/bin/epihap.solve \
    --gfile ../example/example.dat \
    --gmap ../example/example.map \
    --hfile ../example/example.hap \
    --create-grm \
    --egrm-mtd 1 \
    --nthreads 16 \
    --out saved/grm.single &&

echo "finished grms"

/project/dairyxbreed/prakapenka/projects/merge_epihap1_epihap2-atlas/bin/epihap \
    --grmfile saved/grm.single \
    --pfile ../example/example.phen \
    --missing-phen-val -9999111 \
    --trait-col 7 \
    --reml-ce \
    --va 3 \
    --vd 1 \
    --ve 1 \
    --ai-reml-iter-start 3 \
    --nthreads 16 \
    --epochs 10 \
    --out out/epihap
