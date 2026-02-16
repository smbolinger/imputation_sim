#!/bin/bash


if [ $# -eq 0 ]; then
    echo ">> No arguments provided!"
    echo ">> Usage: $0 - [seed] [test | no] [verbosity level] [other arguments passed to isim.R]"
    exit 1
else
    #echo -e "\nargument(s) passed: $1 $2 $3"
    #echo -e -n "\nargument(s) passed to shell script: $@" # 'n' tells it not to add newline at end
    echo -e -n "\n[] [] [] " # 'n' tells it not to add newline at end
fi
#could read seed from file here instead, and then pass to R script:
#seed=$(head -n 1 seed.flag)

#date=$(date +'%d/%b')
#date=$(date +'%d-%b')
date=$(date +'%d%b')
now=$(date +'%H:%M:%S')
FILE="/home/wodehouse/Projects/fate_glm/isim.R"
# this was Wodehoues: cxf
#OUT="${1}${2}.out" # brackets (string concatenation) to prevent '.out' being interpreted as part of var name

# ----------------------------------------------------------------------------------
# 1: functions to run the scripts:
# ----------------------------------------------------------------------------------
simul() {
    local seed="$1"
    local testOn="$2"
    local outf="$3"
    local vb="$4"
    local default="$5"
    #local argList="$5" 
    #local argList=("$@:6") # reassemble args from array as local array, starting at arg 6
    local argList=() # reassemble args from array as local array, starting at arg 6, and add only non-empty vals
    argList+=("$seed" "$vb" "$testOn")
    #echo "args passed to simul: ${@:6}"
    for val in "${@:6}"; do
        #echo "$val"
        if [[ -n "$val" ]]; then
            argList+=("$val")
        fi
    done

    #nohup Rscript /home/wodehouse/Projects/fate_glm/isim.R >> "$out" 2>&1 "s$1" "suff_$suff" "aw$ampWt" "m$nImp" "r$nRun" "j$nSave" "v$vb" "$testOn"&
    #echo "args to pass: ${argList[@]}"

    if [ "$testOn" == "test" ]; then
        #msg="\n>>> running *AS TEST* with seed = $seed & $default args: $argList ; *** verbosity = $vb ***"
        #echo -e "\n>>> running *AS TEST* with seed = $seed & $default args: ${argList[@]} ; *** verbosity = $vb ***"
        echo -e "\n>>> running ***AS TEST*** with: \tSEED = $seed \t& $default ARGS: ${@:6} ; \t*** VERBOSITY = $vb ***\n"
    else
        #msg="\n>>> running with seed = $seed & $default args: $argList ; *** verbosity = $vb ***" # '$argList' is equivalent to $argList[0] - need to expand
        echo -e "\n>>> running with: \tSEED = $seed \t& $default ARGS: ${argList[@]} ; \t*** VERBOSITY = $vb ***\n"
        #msg="\n>>> running with $default args: $argList"
        #echo -e "\n>>> running **AS TEST** with seed $1; other args passed to script: s$1, $testOn, $@; appending output to $OUT"
    fi
    #nohup Rscript /home/wodehouse/Projects/fate_glm/isim.R >> "$outf" 2>&1 "s$seed" "aw3" "m$nImp" "r$nRun" "j1" "v$vb" "$testOn"&
    #nohup Rscript /home/wodehouse/Projects/fate_glm/isim.R >> "$outf" 2>&1 "$seed" "aw3" "j1" "v$vb" "$testOn" "$argList" &
    #nohup Rscript /home/wodehouse/Projects/fate_glm/isim.R >> "$outf" 2>&1 "$seed" "aw3" "j1" "v$vb" "$testOn" "${argList[@]}" &
    #
    # this one isn't expanding the array properly:
    #nohup Rscript /home/wodehouse/Projects/fate_glm/isim.R >> "$outf" 2>&1 $seed aw3 j1 v$vb $testOn "${argList[@]}" &

    #nohup Rscript /home/wodehouse/Projects/fate_glm/isim.R >> "$outf" 2>&1 "${argList[@]}" &
    nohup Rscript /home/wodehouse/Projects/fate_glm/isim.R >> "$outf" 2>&1 "${argList[@]}" 
}

calcb() {
    local modDir="$1"
    local nrun="$2"
    local vb="$3"
    local biasOut="$4"

    nohup Rscript /home/wodehouse/Projects/fate_glm/calc_bias.R >> "$biasOut" 2>&1 "$modDir" "$nrun" "$vb"
}

# ----------------------------------------------------------------------------------
# 2: default values:
# ----------------------------------------------------------------------------------
out="test-${1}.out" # why is test the default?
biasOut="tbias-${1}.out"
mdir="t-${date}"
seed="s$1"
ampWt="aw3"
#suff="aw$ampWt"
testOn="$2"

#if [ "$2" = "test" ]; then 
# ----------------------------------------------------------------------------------
# 3: values for test:
# ----------------------------------------------------------------------------------
if [ "$testOn" = "test" ]; then 
    ###### Default vals for test: ##########################################################################
    nImp="m5"
    nRun="r3"
    nSave="j1"
    vb="v2"

    # check to see if another verbosity value was provided in args:
    for arg in "$@"; do
        if [[ "$arg" =~ [v] ]]; then
            vb="$arg"
            #vb=${arg#*$"v"} # try to extract only the number (everything following v)
            #echo -e "\n\t>>> PID: $! \t\t>>>>>> *VERBOSITY LEVEL*: $arg\n" # -e allows \n to be interpreted as newline
        #else
            #echo -e "\n\t>>> PID: $! \t\t>>>>>> *VERBOSITY LEVEL*: default (0)\n" # -e allows \n to be interpreted as newline
        fi
    done

    if [ ! -z "$4" ]; then
        #simul && calcb
        # do I need to declare argList as local? only if repeating the command?
        argList=("${@:2}") #  start @ arg 2
        #argList+=("$seed")
        def="CLI"
    else
        #argList=("$nImp", "$nRun", "$suff", "$nSave", "$ampWt") # default vals
        argList=("$nImp" "$nRun" "$suff" "$nSave" "$ampWt") # default vals - NO COMMAS!

        #argList+=("$seed")
        def="default"
    fi
    # calcb starts too early - not sure why. maybe best to include bias calculation in isim.R. this happens bc simul was running in the background
    #echo "$argList"
    #(simul "$seed" "$testOn" "$out" "$vb" "$def" "${argList[@]}" && calcb "$out" "$nRun" "$vb" "$biasOut") &
    #
    # this doesn't pass the correct nRun value to calc_bias!
    (simul "$seed" "$testOn" "$out" "$vb" "$def" "${argList[@]}" && calcb "$mdir" "$nRun" "$vb" "$biasOut") &
    #simul "$seed" "$testOn" "$out" "$vb" "$argList" "$def" 
else
    ###### Default vals for non-test: ##########################################################################
    #can prob move them out of the else statement
    out="${1}.out" # why is test the default?
    biasOut="bias-${1}.out"
    mdir="${date}"
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
        #argList="$@"
        # do I need to declare argList as local? only if repeating the command?
        #argList="${@:2}" #  start @ arg 4
        argList=("${@:2}") #  start @ arg 2
        #argList+=("$seed")
        def="CLI"
        #echo -e "\n>>> seed=$1; other args passed to script: s$1, $@; output appended to $OUT"
        #nohup Rscript /home/wodehouse/Projects/fate_glm/isim.R >> "$OUT" 2>&1 "s$1" "$@" & # pass all remaining arguments to script
    else # use defaults
        #argList=("$nImp", "$nRun", "$suff", "$nSave", "$ampWt") # default vals
        argList=("$nImp" "$nRun" "$suff" "$nSave" "$ampWt") # default vals - NO COMMAS!
        #argList+=("$seed")
        def="default"
        #echo -e "\n>>> seed=$1, $nImp imputations, $nRun reps; output appended to $OUT"
        #nohup Rscript /home/wodehouse/Projects/fate_glm/isim.R >> "$OUT" 2>&1 "s$1" "suff_$suff" "aw$ampWt" "m$nImp" "r$nRun" "j$nSave" "v$deb" "$testOn"&
    fi
    #simul "$seed" "$testOn" "$out" "$vb" "$argList" "$def" && calcb "$out" "$nRun" "$vb" "$biasOut"
    #echo "$argList"
    #(simul "$seed" "$testOn" "$out" "$vb" "$argList" "$def" && calcb "$out" "$nRun" "$vb" "$biasOut") &
    (simul "$seed" "$testOn" "$out" "$vb" "$def" "${argList[@]}" && calcb "$mdir" "$nRun" "$vb" "$biasOut") &
    #simul "$seed" "$testOn" "$out" "$vb" "$argList" "$def" 
fi

#sleep 1
echo -e "\n\n\n[] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] []\n" >> "$out"  
#echo -e "\tDATE: $date \t\tTIME: $now \t\tPID: $!\tFILE: $FILE\n" >> "$out"
dash="-"
datestr="${date:0:2}${-}${date:2}"
echo -e "\tDATE: $datestr \t\tTIME: $now \t\tPID: $!\tFILE: $FILE\n" >> "$out"
echo -e "[] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] []\n" >> "$out"  
#echo -e "\n\t>>> PID: $! \t\t>>>>>> *VERBOSITY LEVEL*: $vb\n" # -e allows \n to be interpreted as newline
#echo -e "\n\n\t>>> PID: $! \t\t>>>>>> *VERBOSITY LEVEL*: $vb\n" # -e allows \n to be interpreted as newline
#echo -e -n "\t>>> PID: $! " # -e allows \n to be interpreted as newline, \t as tab, etc

#sleep 1
#echo -e "[] [] [] [] [] [] [] [] output >> $out [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] [] []\n" 
echo -e "[] [] [] [] [] [] [] PID: $! [] [] [] [] [] [] bias >> $biasOut [] [] [] [] [] [] []  output >> $out [] [] [] [] [] [] [] [] [] [] [] [] [] \n" 

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

