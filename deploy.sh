#!/bin/bash
# Este script crea la imagen con todos los prerequisitos 
# para compilar y ejecutar MSOT
#
# Daniel Badagnani, Ushuaia, octubre 2023
# dani.en.villa.rica@gmail.com
#

#docker build --no-cache -t base_croco_msot .
docker pull domarcroco/images-for-croco:base_croco_msot-1.0.0

echo '#/!bin/bash' > run.sh
echo "# Este script inicia la corrida del contenedor y monta la carpeta con MSOT" >> run.sh
echo "#"                                                                          >> run.sh
echo "# Daniel Badagnani, Ushuaia, octubre 2023"                                  >> run.sh
echo "# dani.en.villa.rica@gmail.com"                                             >> run.sh
echo "#"                                                                          >> run.sh
echo "docker run -it -v $PWD/MSOT:/opt/MSOT domarcroco/images-for-croco:base_croco_msot-1.0.0 /bin/bash"      >> run.sh

chmod +x run.sh
