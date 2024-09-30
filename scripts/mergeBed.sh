#!/bin/bash
INFILE=$1 # wildcard input bed6 (+,- strand) paths or paths joined by space or a txt file containing paths
cov_num=$2 # num of sample freq.: [1, ]
chrSize=$3 # genome or tx length file
OUTFILE=$4 # out bed file path

#bedtools v2.30
echo $INFILE
echo $cov_num
echo $chrSize
echo $OUTFILE
reName="TRUE"
echo $reName

# Check if INFILE is a text file
if [[ $INFILE == *.txt ]]; then
  echo "bed path .txt as input"
  PATHS=$(cat $INFILE)
else
  echo "wildercard *.bed as inputs"
  PATHS=$INFILE
fi

# Split strand
for strand in $(echo $PATHS | tr ' ' '\n' | xargs -n1 cat | cut -f6 | sort | uniq); do
  # strand="+"
  echo "start ${strand}"
  if [ $reName == "FALSE" ]; then
    echo $PATHS \
      | xargs -n1 cat \
      | bedtools sort \
      | bedtools genomecov -i - -strand ${strand} -g $chrSize -bg \
      | awk -v s=${strand} 'BEGIN{{OFS="\t";FS="\t"}}{{print $1,$2,$3,"X",$4,s}}' \
      | awk -v c=${cov_num} '$5 >= c' \
      | LC_COLLATE=C sort -k1,1 -k2,2n \
      | bedtools merge -s -c 2,3,5,6 -o collapse,collapse,collapse,collapse \
      | awk 'BEGIN{{OFS="\t";FS="\t"}}
        {{split($4,a,/,/); split($5,b,/,/); split($6,c,/,/); split($7,d,/,/);
        cov=0.0; for(i=1;i<=length(a);i++){{cov+=c[i]*(b[i]-a[i]);}}
        cov /= $3-$2;
        print $1,$2,$3,$1":"$2"-"$3"_"d[1],cov,d[1]
      }}' > ${OUTFILE}.${strand}
  elif [ $reName == "TRUE" ]; then
    echo $PATHS \
      | xargs -n1 cat \
      | bedtools sort \
      | bedtools genomecov -i - -strand ${strand} -g $chrSize -bg \
      | awk -v s=${strand} 'BEGIN{{OFS="\t";FS="\t"}}{{print $1,$2,$3,"X",$4,s}}' \
      | awk -v c=${cov_num} '$5 >= c' \
      | LC_COLLATE=C sort -k1,1 -k2,2n \
      | bedtools merge -s -c 2,3,5,6 -o collapse,collapse,collapse,collapse \
      | awk 'BEGIN{{OFS="\t";FS="\t"}}
        {{split($4,a,/,/); split($5,b,/,/); split($6,c,/,/); split($7,d,/,/);
        cov=0.0; for(i=1;i<=length(a);i++){{cov+=c[i]*(b[i]-a[i]);}}
        cov /= $3-$2;
        print $1,$2,$3,"peak_"NR,cov,d[1]
      }}' > ${OUTFILE}.${strand}
  fi
done

# Merge + and - strands
echo "merge +,-"
if [ $reName == "FALSE" ]; then
  cat ${OUTFILE}.[+-] | bedtools sort > ${OUTFILE}
elif [ $reName == "TRUE" ]; then
  cat ${OUTFILE}.[+-] | bedtools sort | awk 'BEGIN{{OFS="\t";FS="\t"}}{{print $1,$2,$3,"peak_"NR,$5,$6}}' > ${OUTFILE}
fi
rm -f ${OUTFILE}.[+-]
