### Tip 1: R does not supoort sample names with "-", replace it with "_"
### Tip 2: Make sure that all files have the sample ID column with the same name eg.: OTUID 


## .biom to .tsv
biom convert \
  -i "${OUTPUT_DIR}/8-export/*.biom" \
  -o "${OUTPUT_DIR}/8-export/*.tsv" \
  --to-tsv 


## .tsv to .biom
biom convert \
  -i "${OUTPUT_DIR}/9-export/*.tsv" \
  -o "${OUTPUT_DIR}/9-export/*.biom" \
  --table-type="OTU table" \
  --to-hdf5 