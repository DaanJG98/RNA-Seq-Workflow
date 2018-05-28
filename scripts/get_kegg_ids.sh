#!/bin/sh
KEGG_IDs=$(cut -f2 $1)
echo $KEGG_IDs
for ID in $KEGG_IDs
do
  curl -s http://rest.kegg.jp/get/${ID}
done
