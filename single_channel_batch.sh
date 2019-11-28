file=$1
rm -f "${file%.in}.out" "${file%.in}.chi"
< "$1" xargs -t -n 8 -P 4 sh -c './single_channel "$1" "$2" "$3" "$4" "$5" "$6" "$7" "$8" ' argv0
