#/!bin/bash
# Este script inicia la corrida del contenedor y monta la carpeta con MSOT
#
# Daniel Badagnani, Ushuaia, octubre 2023
# dani.en.villa.rica@gmail.com
#
docker run -it -v /home/daniel/2_AREAS/UNTDF/SIMU_OCEANO/CROCOmodels/TdF_MSOT_docker/MSOT:/opt/MSOT domarcroco/images-for-croco:base_croco_msot-1.0.0 /bin/bash
