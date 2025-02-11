#!/bin/zsh
# shell script to launch matlab codes in background

wd=$(pwd)
preamble="cd ${wd};"
# echo $preamble

os=$(uname)

if [[ $os = Darwin ]]; then
	alias matlab='~/Applications/MATLAB_R2023b.app/bin/matlab'
	caffeinate -iw $$ &
fi

foobatch='foodobatch.m'
rm -f $foobatch

doSamStartSPF=false

doSPFquarterlyOnly=false # note: doSPFquarterlyOnly=true makes choices for doY1Q4 obsolete

DATALABELS="{'RGDP','UNRATE','PGDP','CPI'}"

samStart=1 # 52=1981:Q3

MCMCdraws=3000
quicky=false # to be safe, set to false

echo DATALABELS=$DATALABELS
echo MCMCdraws=$MCMCdraws
echo doSamStartSPF=$doSamStartSPF
echo doSPFquarterlyOnly=$doSPFquarterlyOnly


for samStart in 1; do
	for doY1Q4 in true false; do
		if [[ "$doY1Q4" = false ]] && [[ "$doSPFquarterlyOnly" = true ]]; then
			continue
		fi
		for mfile in $@; do
			thismfile=$(basename $mfile .m)
			# echo doY1Q4=$doY1Q4
			for NGAP in "BOP"; do                              # "Ny"
				foofile=foo${thismfile}${NGAP}samStart${samStart} 
				if [[ "$doY1Q4" = false ]]; then
					foofile="${foofile}EXy1q4"
				fi
				if [[ "$doSPFquarterlyOnly" = true ]]; then
					foofile="${foofile}SPFq"
				fi
				length=${#foofile}
				# echo "Length of foofile: $length"
				if [[ $length -gt 63 ]]; then
					echo "**Warning**: Length of foofile exceeds 63 characters."
				fi

				echo $foofile

				# overwrite parameters
				sed '/SED-PARAMETERS-HERE/a\ 
        DATALABELS='$DATALABELS'; \
		MCMCdraws='$MCMCdraws'; \
		quicky='$quicky'; \
        NGAP='\'$NGAP\''; \
		samStart='$samStart'; \
		doY1Q4='$doY1Q4'; \
		doSPFquarterlyOnly='$doSPFquarterlyOnly'; \
        doSamStartSPF='$doSamStartSPF'; \
        ' \
					$mfile >$foofile.m

				# add foofile to batch list
				echo $foofile >>$foobatch

			done # NGAP
		done  # doSPFquarterlyOnly
	done   # samStart
done    # mfiles
