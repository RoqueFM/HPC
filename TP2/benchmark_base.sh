#
# This script is used by the run_* scripts and executed all benchmarks
#

# STRING to gather output in table format (for postprocessing)
CSVTABLE=""
CSVTABLE_HEADER="benchmark\\N"

# This script is intended to be called by other scripts after preparing the parameters!

FIRST_VARIANT=true
for VARIANT in $VARIANT_; do

	CSVTABLE+="\n"
	
	FIRST_N=true
	for N in $N_; do
		# Prepare execution
		CACHE_BLOCKING_SIZE=$(echo "scale=4; sqrt($N)" | bc)
		EXEC="$PROGRAM $VARIANT $N $CACHE_BLOCKING_SIZE"

		if [ ! -z $DEBUG_EXEC ]; then
			if $DEBUG_EXEC; then
				echo -- "$EXEC"
			fi
		fi

		# Execute and gather output
		OUTPUT=`$EXEC` || exit 1

		# Gather header information for CSV output
		if $FIRST_VARIANT; then
			CSVTABLE_HEADER+="\t$N"
		fi

		if $FIRST_N; then
			KERNEL_NAME=$(echo -- "$OUTPUT" | grep " kernel: " | sed "s/.* kernel: //")
			CSVTABLE+="$KERNEL_NAME"
			FIRST_N=false

			echo "**********************************************"
			echo "$KERNEL_NAME"
			echo "**********************************************"

			echo -e "N\tGFLOP/s\tITERATIONS"
		fi

		# Get data from output
		N=$(echo -n "$OUTPUT" | grep " N: " | sed "s/.* N: //" | xargs echo -n)
		echo -en "$N"

		GFLOPS=$(echo -n "$OUTPUT" | grep " g_num_flops/s: " | sed "s/.* g_num_flops\/s: //")
		echo -en "\t$GFLOPS"

		NUM_ITERATIONS=$(echo -- "$OUTPUT" | grep " num_iterations: " | sed "s/.* num_iterations: //")
		echo -en "\t$NUM_ITERATIONS"

		echo ""

		CSVTABLE+="\t$GFLOPS"
	done
	FIRST_VARIANT=false
done


CSVTABLE="$CSVTABLE_HEADER$CSVTABLE"


echo "***"
echo -en "$CSVTABLE"
echo ""
echo ""

OUTPUTFILE="${0/.sh/}.csv"
OUTPUTFILE="output_$(basename $OUTPUTFILE)"

echo "Writing data to file $OUTPUTFILE"
echo -en "$CSVTABLE" > $OUTPUTFILE
