#!/bin/sh

git clone https://bitbucket.org/pflotran/pflotran $PFLOTRAN_DIR
cd $PFLOTRAN_DIR/src/pflotran
git checkout v3.0
#git apply $ALQUIMIA_DIR/.travis/pflotran_version_patch.txt
make libpflotranchem.a
