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

# PART 1: setup
# ----------------------------------------------------------------------------------
#simul(seed, outf, testOn, nImp, nRun, vb, default) {
#simul(seed, testOn, outf, vb, argList, default) { # don't think this syntax is correct
simul() {
    local seed="$1"
    local testOn="$2"
    local outf="$3"
    local vb="$4"
    local argList="$5"
    local default="$6"
    #nohup Rscript /home/wodehouse/Projects/fate_glm/isim.R >> "$out" 2>&1 "s$1" "suff_$suff" "aw$ampWt" "m$nImp" "r$nRun" "j$nSave" "v$vb" "$testOn"&
    if [ "$testOn" == "test" ]; then
        #msg="\n>>> running *AS TEST* with seed = $seed & $default args: $argList ; *** verbosity = $vb ***"
        echo -e "\n>>> running *AS TEST* with seed = $seed & $default args: $argList ; *** verbosity = $vb ***"
    else
        #msg="\n>>> running with seed = $seed & $default args: $argList ; *** verbosity = $vb ***"
        echo -e "\n>>> running with seed = $seed & $default args: $argList ; *** verbosity = $vb ***"
        #msg="\n>>> running with $default args: $argList"
        #echo -e "\n>>> running **AS TEST** with seed $1; other args passed to script: s$1, $testOn, $@; appending output to $OUT"
    fi
    #nohup Rscript /home/wodehouse/Projects/fate_glm/isim.R >> "$outf" 2>&1 "s$seed" "aw3" "m$nImp" "r$nRun" "j1" "v$vb" "$testOn"&
    nohup Rscript /home/wodehouse/Projects/fate_glm/isim.R >> "$outf" 2>&1 "s$seed" "aw3" "j1" "v$vb" "$testOn" "$argList"&
}

calcb() {
    local modDir="$1"
    local nrun="$2"
    local vb="$3"
    local biasOut="$4"

    nohup Rscript /home/wodehouse/Projects/fate_glm/calc_bias.R >> "$biasOut" 2>&1 "$modDir" "$nrun" "$vb"
}

if [ "$2" = "test" ]; then 
    #OUT="test-${1}.out" # all caps vars should be reserved for system vars?

    ###### Default vals for test: ##########################################################################
    out="test-${1}.out"
    biasOut="tbias-${1}.out"
    testOn="test"
    seed="s$1"
    ampWt="aw3"
    #suff="aw$ampWt"
    nImp="m5"
    nRun="r3"
    nSave="j1"
    vb="v2"

    # check to see if another verbosity value was provided in args:
    for arg in "$@"; do
        if [[ "$arg" =~ [v] ]]; then
            #vb="$arg"
            vb=${arg#*$"v"} # try to extract only the number (everything following v)
            #echo -e "\n\t>>> PID: $! \t\t>>>>>> *VERBOSITY LEVEL*: $arg\n" # -e allows \n to be interpreted as newline
        #else
            #echo -e "\n\t>>> PID: $! \t\t>>>>>> *VERBOSITY LEVEL*: default (0)\n" # -e allows \n to be interpreted as newline
        fi
    done

    if [ ! -z "$4" ]; then
        #simul && calcb
        argList="$@" # do I need to declare it as local? only if repeating the command?
        def="CLI"
    else
        argList=("$nImp", "$nRun", "$suff", "$nSave", "$ampWt") # default vals
        def="default"
    fi
    simul "$seed" "$testOn" "$out" "$vb" "$argList" "$def" && calcb "$out" "$nRun" "$vb" "$biasOut"
    #ampWt=1
    #echo ">>> running as test; appending output to $TEST"
    #echo -e "\n>>> running as test; seed=$1, $nImp imputations, $nRun reps; appending output to $OUT"
    #deb=2
    #if [ ! -z "$3" ]; then
    #if [ ! -z "$4" ]; then
    #    echo -e "\n>>> running **AS TEST** with seed $1; other args passed to script: s$1, $testOn, $@; appending output to $OUT"
    #    nohup Rscript /home/wodehouse/Projects/fate_glm/isim.R >> "$OUT" 2>&1 "s$1" "$testOn" "$@" & # pass all remaining arguments to script
    #else # use defaults
    #    echo -e "\n>>> running **AS TEST** with $nRun runs, $nImp imputations, & seed $1; appending output to $OUT"
    #    nohup Rscript /home/wodehouse/Projects/fate_glm/isim.R >> "$OUT" 2>&1 "s$1" "suff_$suff" "aw$ampWt" "m$nImp" "r$nRun" "j$nSave" "v$vb" "$testOn"&
    #fi

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

        ###### Default vals for non-test: ##########################################################################
    out="${1}.out" # brackets (string concatenation) to prevent '.out' being interpreted as part of var name
    biasOut="bias-${1}.out"
    testOn=""
    seed="s$1"
    ampWt="aw3"
    #suff="aw$ampWt"
    nImp="m25"
    nRun="r100"
    nSave="j100"
    vb="default (0)"

        # check to see if another verbosity value was provided in args:
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
            argList="$@"
            def="CLI"
            #echo -e "\n>>> seed=$1; other args passed to script: s$1, $@; output appended to $OUT"
            #nohup Rscript /home/wodehouse/Projects/fate_glm/isim.R >> "$OUT" 2>&1 "s$1" "$@" & # pass all remaining arguments to script
        else # use defaults
            argList=("$nImp", "$nRun", "$suff", "$nSave", "$ampWt") # default vals
            def="default"
            #echo -e "\n>>> seed=$1, $nImp imputations, $nRun reps; output appended to $OUT"
            #nohup Rscript /home/wodehouse/Projects/fate_glm/isim.R >> "$OUT" 2>&1 "s$1" "suff_$suff" "aw$ampWt" "m$nImp" "r$nRun" "j$nSave" "v$deb" "$testOn"&
        fi
        simul "$seed" "$testOn" "$out" "$vb" "$argList" "$def" && calcb "$out" "$nRun" "$vb" "$biasOut"
fi

echo -e "\n[] [] [] [] [] [] [] [] output >> $out [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] []\n" 
echo -e "\n[] [] [] [] [] [] [] [] bias >> $biasOut [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] []\n" 

echo -e "\n\n\n[] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] []\n" >> "$out"  
echo -e "\tDATE: $DATE \t\tTIME: $NOW \t\tPID: $!\tFILE: $FILE\n" >> "$out"
echo -e "[] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] []\n" >> "$out"  
echo -e "\n\t>>> PID: $! \t\t>>>>>> *VERBOSITY LEVEL*: $vb\n" # -e allows \n to be interpreted as newline

# PART 2: run the commands
# ----------------------------------------------------------------------------------

#$simulation && 
#nohup Rscript /home/wodehouse/Projects/fate_glm/calc_bias.R >> "$biasOut" 2>&1 "$OUT" "$nRun" v1
# need to wait until other script is finished
#echo -e "\n[] [] [] [] [] [] [] [] bias >> $biasOut ; PID = $! [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] []\n" 


#nohup Rscript /home/wodehouse/Projects/fate_glm/calc_bias.R >> "$biasOut" 2>&1 "$OUT" "$nRun" v1
#echo -e "\n\n\n[] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] []\n" >> "$biasOut"  
#echo -e ">> calculating bias"  >> "$biasOut"
#echo -e "\n[] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] []\n" >> "$biasOut"  

