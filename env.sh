#!/bin/bash


export PATH=$PATH:/home/caixx/git/ncplugin-PiXiu/testcode/utils/cache/ncrystal/install/bin/
. testcode/utils/bootstrap.sh


export NCRYSTAL_PLUGIN_LIST="/home/caixx/git/ncplugin-PiXiu/testcode/utils/cache/bld/libNCPlugin_PiXiu.so:$NCRYSTAL_PLUGIN_LIST"
ncrystal-config --setup

#export NCrystal_DIR=/home/caixx/git/ncplugin-PiXiu/venv/lib/python3.8/site-packages/NCrystal/ncrystal_pyinst_data/lib/cmake/NCrystal
#export LD_LIBRARY_PATH=/home/caixx/git/ncplugin-PiXiu/venv/lib/python3.8/site-packages/NCrystal/ncrystal_pyinst_data/lib:/home/caixx/git/ncplugin-PiXiu/testcode/utils/cache/install/lib:/home/caixx/git/ncplugin-PiXiu/testcode/utils/cache/bld
