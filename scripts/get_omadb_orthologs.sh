#!/bin/sh
while read id
do
  data=$(curl -s "https://omabrowser.org/api/protein/${id}/orthologs/" | jq -r '.[0:10]')
  print_id=">$id"
  echo $print_id >> $2

  for row in $(echo "${data}" | jq -r '.[] | @base64'); do
    _jq() {
      echo ${row} | base64 --decode | jq -r ${1}
    }
    if [ -n "$(_jq '.canonicalid')" ]; then
      echo $(_jq '.canonicalid') >> $2
    fi
  done
done
