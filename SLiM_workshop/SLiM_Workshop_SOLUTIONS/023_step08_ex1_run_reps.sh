for N in 100 200 500 1000 ; do
	for SEED in 1 2 3 ; do
		echo "Running with N ==" ${N}:
		slim -s ${SEED} -d MU=1e-7 -d R=1e-8 -d N=${N} -d G="N*10" rep.slim > out_N${N}_SEED${SEED}.txt
		echo
	done
done
