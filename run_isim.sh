#!/bin/bash


if [ $# -eq 0 ]; then
    echo ">> No arguments provided!"
    echo ">> Usage: $0 - [seed] [test | no] [other arguments passed to isim.R]"
    exit 1
else
    #echo -e "\nargument(s) passed: $1 $2 $3"
    echo -e "\nargument(s) passed to shell script: $@"
fi
#could read seed from file here instead, and then pass to R script:
#seed=$(head -n 1 seed.flag)

DATE=$(date +'%d/%b')
NOW=$(date +'%H:%M:%S')
FILE="/home/wodehouse/Projects/fate_glm/isim.R"
# this was Wodehoues: cxf
#OUT="${1}${2}.out" # brackets (string concatenation) to prevent '.out' being interpreted as part of var name
#TEST=5 
if [ "$2" = "test" ]; then 
    OUT="test-${1}.out"
    testOn="test"
    ampWt=3
    suff="aw$ampWt"
    nImp=5
    nRun=3
    nSave=1
    vb=2
    for arg in "$@"; do
        if [[ "$arg" =~ [v] ]]; then
            #vb="$arg"
            vb=${arg#*$"v"} # try to extract only the number (everything following v)
            #echo -e "\n\t>>> PID: $! \t\t>>>>>> *VERBOSITY LEVEL*: $arg\n" # -e allows \n to be interpreted as newline
        #else
            #echo -e "\n\t>>> PID: $! \t\t>>>>>> *VERBOSITY LEVEL*: default (0)\n" # -e allows \n to be interpreted as newline
        fi
    done
    #ampWt=1
    #echo ">>> running as test; appending output to $TEST"
    #echo -e "\n>>> running as test; seed=$1, $nImp imputations, $nRun reps; appending output to $OUT"
    #deb=2
    #if [ ! -z "$3" ]; then
    if [ ! -z "$4" ]; then
        echo -e "\n>>> running **AS TEST** with seed $1; other args passed to script: s$1, $testOn, $@; appending output to $OUT"
        nohup Rscript /home/wodehouse/Projects/fate_glm/isim.R >> "$OUT" 2>&1 "s$1" "$testOn" "$@" & # pass all remaining arguments to script
    else # use defaults
        echo -e "\n>>> running **AS TEST** with $nRun runs, $nImp imputations, & seed $1; appending output to $OUT"
        nohup Rscript /home/wodehouse/Projects/fate_glm/isim.R >> "$OUT" 2>&1 "s$1" "suff_$suff" "aw$ampWt" "m$nImp" "r$nRun" "j$nSave" "v$vb" "$testOn"&
    fi

    #if [ "$3" = "v2" ]; then
    #    deb=2
    #elif [ "$3" = "v3" ]; then
    #    deb=3
    #else
    #    deb=1
    #fi
    #if [ ! -z "$4" ]; then
    #    nRun="$4"
    #fi
else
	#echo "running with $1 as name: $IMP imputations, $REP reps; output to $OUT"
        OUT="${1}.out" # brackets (string concatenation) to prevent '.out' being interpreted as part of var name
        ampWt=3
        suff="aw$ampWt"
        nImp=25
        nRun=100
        nSave=100
        testOn=""
        #deb=0
        #if [ "$2" = "v1" ]; then
        #    deb=1
        #else
        #    deb=0
        #fi   
        vb="default (0)"
        for arg in "$@"; do
            if [[ "$arg" =~ [v] ]]; then
                #vb="$arg"
                vb=${arg#*$"v"} # try to extract only the number (everything following v)
                #echo -e "\n\t>>> PID: $! \t\t>>>>>> *VERBOSITY LEVEL*: $arg\n" # -e allows \n to be interpreted as newline
            #else
                #echo -e "\n\t>>> PID: $! \t\t>>>>>> *VERBOSITY LEVEL*: default (0)\n" # -e allows \n to be interpreted as newline
            fi
        done
        if [ ! -z "$2" ]; then
            echo -e "\n>>> seed=$1; other args passed to script: s$1, $@; output appended to $OUT"
            nohup Rscript /home/wodehouse/Projects/fate_glm/isim.R >> "$OUT" 2>&1 "s$1" "$@" & # pass all remaining arguments to script
        else # use defaults
            echo -e "\n>>> seed=$1, $nImp imputations, $nRun reps; output appended to $OUT"
            nohup Rscript /home/wodehouse/Projects/fate_glm/isim.R >> "$OUT" 2>&1 "s$1" "suff_$suff" "aw$ampWt" "m$nImp" "r$nRun" "j$nSave" "v$deb" "$testOn"&
        fi

fi

#nohup Rscript /home/wodehouse/Projects/fate_glm/isim.R >> "$OUT" 2>&1 "s$1" "m5" "r3" "j1" "deb" "xdeb" "test" &
# most of these args get handled in the isim.R script - not sure if that's better or worse than here:
#if [ "$deb"==1 ]; then
#if [ "$deb" -eq 1 ]; then
#    nohup Rscript /home/wodehouse/Projects/fate_glm/isim.R >> "$OUT" 2>&1 "s$1" "suff_$suff" "aw$ampWt" "m$nImp" "r$nRun" "j$nSave" "$testOn"&
#elif [ "$deb" -eq 2 ]; then
#    nohup Rscript /home/wodehouse/Projects/fate_glm/isim.R >> "$OUT" 2>&1 "s$1" "suff_$suff" "aw$ampWt" "m$nImp" "r$nRun" "j$nSave" "v$deb" "$testOn"&
#elif [ "$deb" -eq 3 ]; then
#    #echo -e "extra debug!\n"
#    nohup Rscript /home/wodehouse/Projects/fate_glm/isim.R >> "$OUT" 2>&1 "s$1" "suff_$suff" "aw$ampWt" "m$nImp" "r$nRun" "j$nSave" "v$deb" "xdeb" "$testOn"&
#else

#fi

echo -e "\n\n\n[] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] []\n" >> "$OUT"  
#echo -e "\tDATE: $DATE \t\tTIME: $NOW \t\tPID: $! \n" >> "$OUT"
#echo -e "\tFILE: $FILE\n" >> "$OUT"
echo -e "\tDATE: $DATE \t\tTIME: $NOW \t\tPID: $!\tFILE: $FILE\n" >> "$OUT"
echo -e "[] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] []\n" >> "$OUT"  
#echo -e "[][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][]\n\n" >> "$OUT"  

#echo -e "\n\t>>> PID: $! \t\t>>>>>> *VERBOSITY LEVEL*: $deb\n" # -e allows \n to be interpreted as newline
#vb=0

#echo -e "\n\t>>> PID: $! \t\t>>>>>> *VERBOSITY LEVEL*: $deb\n" # -e allows \n to be interpreted as newline
echo -e "\n\t>>> PID: $! \t\t>>>>>> *VERBOSITY LEVEL*: $vb\n" # -e allows \n to be interpreted as newline
#echo -e "\t>>>>>> *VERBOSITY LEVEL*: $deb\n"
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
