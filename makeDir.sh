
#!/bin/bash

printf "\nMake structure directory: data, barcodes refgenome, results, logs\n"
for i in data barcodes refgenome results logs
do
mkdir $i
done
