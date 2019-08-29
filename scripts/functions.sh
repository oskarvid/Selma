#!/bin/bash

# Function that prints help message
usage () {
echo "${GREEN}""Usage:${RESET} $0 [ -h ] [ -i input directory ] [ -t tsv file ] [ -o output directory ] [ -r genome version] [ -l interval file ] \


-h : Print this help message \

-r : Human genome reference: b37 or hg38, required \

-i : Path of the input file directory, path must be absolute, required \

-t : Tab separated file with fastq info, path must be absolute, required \

-o : Path of the output folder, path must be absolute, required \

-l : Optional interval file, path must be absolute \
" 1>&2; exit 0; }

# Make your terminal text output beautiful <3
# The following variables are used to define the color in the functions below
BLACK=$(tput setaf 0)
RED=$(tput setaf 1)
GREEN=$(tput setaf 2)
YELLOW=$(tput setaf 3)
BLUE=$(tput setaf 4)
MAGENTA=$(tput setaf 5)
CYAN=$(tput setaf 6)
WHITE=$(tput setaf 7)
RESET=$(tput sgr 0)

showdate () {
	local DATE
	DATE=$(date "+%F %T")
	printf "$DATE"
}

err () {
	local DATE
	DATE=$(showdate)
	printf "[$DATE][${RED}ERROR${RESET}]: $1\n" 1>&2
}

warn () {
	local DATE
	DATE=$(showdate)
	printf "[$DATE][${YELLOW}WARNING${RESET}]: $1\n"
}

inf () {
	local DATE
	DATE=$(showdate)
	printf "[$DATE][INFO]: $1\n"
}

succ () {
	local DATE
	DATE=$(showdate)
	printf "[$DATE][${GREEN}SUCCESS${RESET}]: $1\n"
}

rprog () {
	rsync -ah --progress $@
}

selma-ascii () {
inf "Starting Selma"
cat << "EOF"
                 _..--+~/@-~--.
             _-=~      (  .   "}
          _-~     _.--=.\ \""""
        _~      _-       \ \_\
       =      _=          '--'
      '      =                             .
     :      :       ____                   '=_. ___
___  |      ;                            ____ '~--.~.
     ;      ;                               _____  } |
  ___=       \ ___ __     __..-...__           ___/__/__
     :        =_     _.-~~          ~~--.__
_____ \         ~-+-~                   ___~=_______
     ~@#~~ == ...______ __ ___ _--~~--_
EOF
}
