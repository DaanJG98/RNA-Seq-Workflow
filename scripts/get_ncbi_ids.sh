#!/bin/bash
while read id
do
  data=$(curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gene&term=${id}")
  title=$(grep -oPm1 "(?<=<Id>)[^<]+" <<< "$data")
  echo $title >> $2
done
