#!/bin/bash
######################################################
##
## MSOT TdF v1.0 Universidad Nacional de
##               Tierra del Fuego,
##               Ushuaia, ARGENTINA
##
## Abril 2024
##
## Daniel Badagnani,               Monica Manceñido
## dani.en.villa.rica@gmail.com    yavamoni@gmail.com
## 
#######################################################
export MSOT_HOME=`pwd`
if [[ ":$PATH:" != *":$MSOT_HOME:"* ]]
then
    export PATH="$PATH:$MSOT_HOME"
fi
