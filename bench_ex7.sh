#!/usr/bin/env bash

set -euo pipefail

CONFIG=${CONFIG:-config.txt}
OUT=${OUT:-ex7_results.csv}
MODE=${MODE:-strong}
PROCS=${PROCS:-"1 2 4 8"}
EXERCISES=${EXERCISES:-"1 2 3 4 5 6"}
REPEATS=${REPEATS:-1}
MPIRUN=${MPIRUN:-mpirun}
MPI_ARGS=${MPI_ARGS:-"--use-hwthread-cpus"}

if [[ "${MODE}" != "strong" && "${MODE}" != "weak" ]]; then
	echo "MODE must be either strong or weak" >&2
	exit 1
fi

make

echo "mode,exercise,processes,repeat,scale,total_time_seconds" > "${OUT}"

for exercise in ${EXERCISES}; do
	for processes in ${PROCS}; do
		scale=1
		if [[ "${MODE}" == "weak" ]]; then
			scale=${processes}
		fi

		for repeat in $(seq 1 "${REPEATS}"); do
			echo "Running mode=${MODE} exercise=${exercise} np=${processes} repeat=${repeat} scale=${scale}" >&2
			output=$(${MPIRUN} ${MPI_ARGS} -np "${processes}" ./lbm -n -c "${CONFIG}" -e "${exercise}" -s "${scale}")
			total_time=$(printf '%s\n' "${output}" | awk '/Total time:/ {print $3}')

			if [[ -z "${total_time}" ]]; then
				echo "Could not parse total time for exercise=${exercise}, np=${processes}" >&2
				exit 1
			fi

			echo "${MODE},${exercise},${processes},${repeat},${scale},${total_time}" | tee -a "${OUT}"
		done
	done
done

echo "Wrote ${OUT}" >&2
