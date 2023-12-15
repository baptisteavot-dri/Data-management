IFS=, read -ra fields <<<"${1}"
# Loop over fields from the back, down to the 2nd.
for (( i = 0; i <= 100; i++ )); do
    if [[ "${fields[i]}" == *'synapse_mirror/'* ]]; then
        concat_file="${fields[i]}"
	count=$i
    fi
done

IFS=, read -ra fields <<<"${1}"
for (( i = 0; i < $count ; i++ )); do
    if [[ "${fields[$count]}" == *'synapse_mirror/'* ]]; then
	zcat "${fields[i]}" | gzip >> $concat_file
  	#echo $concat_file  
    fi
done

