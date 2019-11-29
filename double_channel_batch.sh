file=$1
rm -f "${file%.in}.out" "${file%.in}.chi"
< "$1" xargs -t -n 11 -P $(($(nproc)-1)) sh -c './double_channel "$1" "$2" "$3" "$4" "$5" "$6" "$7" "$8" "$9" "$10" "$11" ' argv0
