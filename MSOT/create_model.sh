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
if [ $(find . -maxdepth 0 -empty) ]; then
    echo
    echo "######################################################"
    echo "##"
    echo "## MSOT TdF v1.0 Universidad Nacional de"
    echo "##               Tierra del Fuego,"
    echo "##               Ushuaia, ARGENTINA"
    echo "##"
    echo "## Abril 2024"
    echo "##"
    echo "## Daniel Badagnani,               Monica Manceñido"
    echo "## dani.en.villa.rica@gmail.com    yavamoni@gmail.com"
    echo "## "
    echo "#######################################################"
    echo
    echo ">> COPIANDO ARCHIVOS TEMPLATE"
    cp -fr $MSOT_HOME/TEMPLATE/* .
    echo
    echo ">> Listo. ¡Que te diviertas!"
else
    echo
    echo "ERROR: Esta carpeta no está vacía."
    echo "No se copió ningún archivo"
    echo
fi
