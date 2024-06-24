## .csv to .tsv and back
cat "${METADATA}/metadata.csv" | sed 's/;/\t/g' > "${METADATA}/metadata.tsv" 
cat "${METADATA}/metadata.tsv" | sed 's/;/\t/g' > "${METADATA}/metadata.csv" 


