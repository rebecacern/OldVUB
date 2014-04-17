#!/bin/bash

maxRetries=5
timeout=200
nThreads=3

files=()

# checking if sample is staged on https://mtop.iihe.ac.be/samples for fall back download option

if [ -f "folders" ];then
    rm -v folders
fi

name="folders_$(date +%s).txt"
curl -k -o $name https://mtop.wn.iihe.ac.be/samples/folders >/dev/null 2>&1

mtopstaged=$(grep -R $1 $name | wc -l)

rm -v $name

# get the final list of files to download
for i in $(cd $1; ls *.root);do

    # option 2 allows to give a list of files so we can download a subset
    if [ $# -gt 1 ];then
	
	IFS=',' read -a subset <<< "$2"

	for (( j=0; j<${#subset[@]}; ++j ))
	  do
	  
	  #echo $i ${subset[$j]}
	  if [ "$i" == "${subset[$j]}" ];then

	      files+=($i)

	  fi
	  
	done
    
    else

	files+=($i)
	
    fi

done

#echo ${#files[@]}

#exit

for (( i=0; i<${#files[@]}; i+=$nThreads ))
  do

  tfiles=()
  pid=()
  start=()
  retries=()

  let max=$i+nThreads

  if [ "$max" -gt "${#files[@]}" ];then
      max=${#files[@]}
  fi

  #echo $i $max

  for (( j=$i; j<$max; ++j ))
    do

    let n=$j+1
    echo "$(date) Downloading file: ${files[$j]} ($n/${#files[@]})"

    dccp dcap://maite.iihe.ac.be"$1"${files[$j]} ./${files[$j]} &

    pid+=($!)

    start+=($(date +%s))

    retries+=(0)

    tfiles+=(${files[$j]})
  
  done

  stop=0

  while (( $stop == 0 ))
    do

    anyRunning=0
    for (( k=0; k<${#pid[@]}; ++k ))
      do

      kill -0 ${pid[$k]} >/dev/null 2>&1

      pidrunning=$?

      if [ $pidrunning == 0 ];then
	  anyRunning=1
      else
	  continue
      fi

      now=$(date +%s)
      let running=$now-${start[$k]}

      #echo ""
      #echo "PID IS : ${pid[$k]}"
      #echo "RUNNING SINCE : $running"
      #echo "RETRIES IS : ${retries[$k]}"
      #echo ""

      if [ "$running" -gt "$timeout" ]; then

	  echo "$(date)  !-->! Download of ${tfiles[$k]} got stuck, killing process"

	  kill -9 ${pid[$k]} >/dev/null 2>&1
	  
	  if [ "${retries[$k]}" -lt "$maxRetries" ];then

	      let retries[$k]=${retries[$k]}+1
	      let start[$k]=$(date +%s)

	      echo "$(date)    !-->! Retrying download, attempt ${retries[$k]}/$maxRetries"

	      if [ "$mtopstaged" == "1" ];then
		  echo "$(date)      !-->! Sample is staged on HTTP, trying to download from mtop"
		  IFS='/' read -a array <<< "$1"
		  let pos=${#array[@]}-2
		  folder=${array[$pos]}
		  curl -k -o ${tfiles[$k]} https://mtop.wn.iihe.ac.be/samples/$folder/${tfiles[$k]}
	      else
		  #sleep 200 &
		  dccp dcap://maite.iihe.ac.be"$1"${tfiles[$k]} ./${tfiles[$k]} &
		  let pid[$k]=$!
	      fi
	      
	  fi

      fi

    done

    if [ "$anyRunning" -eq 0 ];then
	stop=1
    else
	sleep 5
    fi

  done

done
