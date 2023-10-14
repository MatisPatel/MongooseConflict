block=0
n_per_job=2
while (( block * n_per_job < $(wc -l < "params.tab") )); do
  while read param1 param2; do
    echo "julia hello_world.jl ${param1} ${param2}"
  done<<EOF > "runlist${block}.txt"
    $(awk "NR > $(( block * n_per_job )) && NR <= $(( (block + 1) * n_per_job ))" "params.tab")
EOF
  parallel < runlist${block}.txt
  block=$(( block + 1 ))
done