#!/bin/bash


echo "argument(s) passed: $1 $2"
#could read seed from file here instead, and then pass to R script:
#seed=$(head -n 1 seed.flag)

DATE=$(date +'%d/%b')
NOW=$(date +'%H:%M:%S')
# this was Wodehoues: cxf
#OUT="${1}${2}.out" # brackets (string concatenation) to prevent '.out' being interpreted as part of var name
OUT="${1}.out" # brackets (string concatenation) to prevent '.out' being interpreted as part of var name
#TEST=5 
if [ "$2" = "test" ]; then 
	TEST="test-${1}.out"
	echo ">>> running as test; saving to $TEST"
	nohup Rscript /home/wodehouse/Projects/fate_glm/isim.R >> "$TEST" 2>&1 "$1" "s333" "m5" "r3" "j1" "deb" "test" &
	echo ">>> args: $1 s333 m5 r3 j1 deb test"
	echo -e "\t\t>>> PID: $! \n"
	echo -e "\n\n\n[][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][]\n" >> "$TEST"  
	echo -e "\tDATE: $DATE \t\tTIME: $NOW \t\tPID: $! \n" >> "$OUT"
	#echo -e "\t\tDATE: $NOW \t\tPID: $! \n" >> "$TEST"
	echo -e "[][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][]\n\n" >> "$TEST"  

else
	IMP=70
	REP=100
	SAVE=50
	nohup Rscript /home/wodehouse/Projects/fate_glm/isim.R >> "$OUT" 2>&1 "m$IMP" "$1" "r$REP" "j$SAVE" & 
	echo "running with $1 as response variable: $IMP imputations, $REP reps; output to $OUT"
	echo -e "\t\t>>> PID: $! \n"
	# -e allows \n to be interpreted as newline
	echo -e "\n\n\n[][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][]\n" >> "$OUT"  
	#echo -e "\t\tPID: $! \n" >> "$OUT"
	echo -e "\tDATE: $DATE \t\tTIME: $NOW \t\tPID: $! \n" >> "$OUT"
	echo -e "[][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][]\n\n" >> "$OUT"  
fi

# use flag file if you want something to run only once per day (like printing DATE)
#echo "$DATE" > isim_flag.out
#if [ "$1"="HF_mis" ]; then
#	echo "running with HF_mis as response variable: 20 imputations, 500 reps "
#	nohup Rscript /home/wodehouse/Projects/fate_glm/isim.R >> hf_mis.out 2>&1 m20 "HF_mis" & 
#	echo "PID $NOW: $!" >> hf_mis.out

#elif [ "$1"="is_u" ]; then
#	echo "running with HF_mis as response variable: 20 imputations, 500 reps "
#	nohup Rscript /home/wodehouse/Projects/fate_glm/isim.R >> is_u.out 2>&1 m20 "is_u" & 
#	echo "PID $NOW: $!" >> is_u.out
#else
#	echo "error: response variable argument should be one of: 'HF_mis' 'is_u'"
#fi

#echo $OUT

# If I could figure out how to create an 'outfile' object, could be much shorter:

# outfile =?
