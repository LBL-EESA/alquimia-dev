#!/bin/sh

git clone https://bitbucket.org/crunchflow/crunchtope-dev.git --branch petsc-upgrade $CRUNCHTOPE_DIR
cd $CRUNCHTOPE_DIR/source
git apply MakefileForAlquimia.patch
make libcrunchchem.a
