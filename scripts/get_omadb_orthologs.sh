#!/bin/bash
while read id
do
  data=$(curl -s "https://omabrowser.org/api/protein/${id}/orthologs/")
  id=">$id"
  echo $id
  
  for row in $(echo "${sample}" | jq -r '.[] | @base64'); do
    _jq() {
     echo ${row} | base64 --decode | jq -r ${1}
    }

   if [ -n "$(_jq '.canonicalid')" ]; then
	echo $(_jq '.canonicalid')
   fi
done