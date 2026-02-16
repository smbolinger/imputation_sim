
#!/bin/bash

if [ $# -eq 0 ]; then
    echo ">> No arguments provided!"
    echo ">> Usage: $0 - [directory in out/] [nrun] [verbosity]"
    exit 1
else
    #echo -e "\nargument(s) passed: $1 $2 $3"
    #echo -e -n "\nargument(s) passed to shell script: $@" # 'n' tells it not to add newline at end
    echo -e -n "\n[] [] [] " # 'n' tells it not to add newline at end
fi

date=$(date +'%d%b')
now=$(date +'%H:%M:%S')
file="/home/wodehouse/Projects/fate_glm/calc_bias.R"

calcb() {
    local modDir="$1"
    local nrun="$2"
    local vb="$3"
    local biasOut="$4"

    nohup Rscript /home/wodehouse/Projects/fate_glm/calc_bias.R >> "$biasOut" 2>&1 "$modDir" "$nrun" "$vb"
}

