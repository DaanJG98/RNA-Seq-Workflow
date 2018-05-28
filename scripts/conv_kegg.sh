#!/bin/sh
while read NCBI_ID
do
  curl -s http://rest.kegg.jp/conv/genes/ncbi-geneid:${NCBI_ID} >> $2
done
