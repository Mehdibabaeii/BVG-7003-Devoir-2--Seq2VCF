
#!/bin/bash

printf "\nMake structure directory: data, barcodes refgenome, results, log\n"
for i in data reject barcodes refgenome results stat
do
mkdir $i
done
