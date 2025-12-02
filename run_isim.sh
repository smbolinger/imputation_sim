#!/bin/bash


echo "argument(s) passed: $1 $2"
#could read seed from file here instead, and then pass to R script:
#seed=$(head -n 1 seed.flag)

DATE=$(date +'%d/%b')
NOW=$(date +'%H:%M:%S')
# this was Wodehoues: cxf
#OUT="${1}${2}.out" # brackets (string concatenation) to prevent '.out' being interpreted as part of var name
#TEST=5 
if [ "$2" = "test" ]; then 
	OUT="test-${1}.out"
        nImp=5
        nRun=3
        nSave=1
	#echo ">>> running as test; appending output to $TEST"
	echo "\n>>> running as test; seed=$1, $nImp imputations, $nRun reps; appending output to $OUT"
        deb=1
        testOn="test"
else
	#echo "running with $1 as name: $IMP imputations, $REP reps; output to $OUT"
        OUT="${1}.out" # brackets (string concatenation) to prevent '.out' being interpreted as part of var name
        nImp=25
        nRun=100
        nSave=50
        testOn=""

        if [ "$2" = "debug" ]; then
            deb=1
        else
            deb=0
        fi   
	echo "\n>>> seed=$1, $nImp imputations, $nRep reps; output appended to $OUT"
fi

#nohup Rscript /home/wodehouse/Projects/fate_glm/isim.R >> "$OUT" 2>&1 "s$1" "m5" "r3" "j1" "deb" "xdeb" "test" &
# most of these args get handled in the isim.R script - not sure if that's better or worse than here:
#if [ "$deb"==1 ]; then
if [ "$deb" -eq 1 ]; then
    nohup Rscript /home/wodehouse/Projects/fate_glm/isim.R >> "$OUT" 2>&1 "s$1" "m$nImp" "r$nRun" "j$nSave" "deb" "$testOn"&
elif [ "$deb" -eq 2 ]; then
    echo -e "extra debug!\n"
    nohup Rscript /home/wodehouse/Projects/fate_glm/isim.R >> "$OUT" 2>&1 "s$1" "m$nImp" "r$nRun" "j$nSave" "deb" "xdeb" "$testOn"&
else
    nohup Rscript /home/wodehouse/Projects/fate_glm/isim.R >> "$OUT" 2>&1 "s$1" "m$nImp" "r$nRun" "j$nSave" "$testOn"&
fi

echo -e "\n\n\n[][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][]\n" >> "$OUT"  
echo -e "\tDATE: $DATE \t\tTIME: $NOW \t\tPID: $! \n" >> "$OUT"
echo -e "[][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][]\n\n" >> "$OUT"  

echo -e "\n\t>>> PID: $! " # -e allows \n to be interpreted as newline
echo -e "\t>>>>>> *DEBUG LEVEL*: $deb\n"
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
