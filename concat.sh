if [[ $3 == *'synapse_mirror/'* ]]; then
	zcat $1 $2 | gzip > ${3}
fi

if [[ $9 == *'synapse_mirror/'* ]]; then
	zcat $1 $2 $3 $4 $5 $6 $7 $8 | gzip > ${9}
fi

if [[ ${6} == *'synapse_mirror/'* ]]; then
        zcat $1 $2 $3 $4 $5 | gzip > ${6}
fi

if [[ $5 == *'synapse_mirror/'* ]]; then
        zcat $1 $2 $3 $4 | gzip > ${5}
fi

if [[ $7 == *'synapse_mirror/'* ]]; then
        zcat $1 $2 $3 $4 $5 $6 | gzip > ${7}
fi

if [[ $8 == *'synapse_mirror/'* ]]; then
        zcat $1 $2 $3 $4 $5 $6 $7 | gzip > ${8}
fi

