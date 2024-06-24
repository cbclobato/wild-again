## .csv to .fasta
cat "${BC}/bcpr-fw.csv" | sed 's/^/>/' | tr ";" "\n" > "${BC}/bcpr-fw.fasta"
cat "${BC}/bcpr-rv.csv" | sed 's/^/>/' | tr ";" "\n" > "${BC}/bcpr-rv.fasta"

