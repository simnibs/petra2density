#!/bin/bash

if [[ $# -lt 3 ]]; then
    echo "Usage: ./petra2density.sh [-p] subject_id T1_file PETRA_file" >&2
    echo "Please provide the subject_id, T1_file, and PETRA_file" >&2
    echo "use -p to register to the petra. default is to register to the t1" >&2
    exit 2
fi



PETRA2DENSITY_DIR="$(dirname "$(which "$0")")"
REGISTER_TO_PETRA=0

while getopts "h:p" opt; do
  case ${opt} in
    h) echo "./petra2density.sh [-p] subject_id T1_file PETRA_file"; exit 1;;
    p) echo "Registering to petra image"
       REGISTER_TO_PETRA=1
       ;;
    ?) echo "Invalid option: -${OPTARG}."; exit 1 ;;
  esac
done

shift $((OPTIND - 1))

if [ "$REGISTER_TO_PETRA" -eq 1 ]; then
    charm "$1" "$3" "$2" --usesettings $PETRA2DENSITY_DIR/charm.ini --forceqform
    python $PETRA2DENSITY_DIR/petra2density.py "./m2m_$1" --register_to_petra
    exit
fi

charm "$1" "$2" "$3" --usesettings $PETRA2DENSITY_DIR/charm.ini --forceqform
python $PETRA2DENSITY_DIR/petra2density.py "./m2m_$1"

