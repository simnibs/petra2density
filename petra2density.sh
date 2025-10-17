#!/bin/bash

if [[ $# -lt 3 ]]; then
    echo "Usage: ./petra2density.sh [-c /calibration_file.csv] [-p] subject_id T1_file PETRA_file" >&2
    echo "Please provide the subject_id, T1_file, and PETRA_file" >&2
    echo "use -p to register to the petra. default is to register to the t1" >&2
    echo "use -c /path/to/calibration.csv to to use a custom calibration file." >&2
    exit 2
fi

CT2DENSITY_CALIBRATION_FILE='none'
PETRA2DENSITY_DIR="$(dirname "$(which "$0")")"
REGISTER_TO_PETRA=0

while getopts ":hpc:" opt; do
  case ${opt} in
    h) echo "./petra2density.sh [-p] subject_id T1_file PETRA_file"; exit 1;;
    p) echo "Registering to petra image"
       REGISTER_TO_PETRA=1
       ;;
    c) echo "Using calibration file"
       CT2DENSITY_CALIBRATION_FILE="$OPTARG"
       ;;
    ?) echo "Invalid option: -${OPTARG}."; exit 1 ;;
  esac
done

shift $((OPTIND - 1))

if [ "$REGISTER_TO_PETRA" -eq 1 ]; then
    charm "$1" "$3" "$2" --usesettings $PETRA2DENSITY_DIR/charm.ini --forceqform
    python $PETRA2DENSITY_DIR/petra2density.py "./m2m_$1" --register_to_petra --ct2density_calibration_file $CT2DENSITY_CALIBRATION_FILE
    exit
fi

charm "$1" "$2" "$3" --usesettings $PETRA2DENSITY_DIR/charm.ini --forceqform
python $PETRA2DENSITY_DIR/petra2density.py "./m2m_$1" --ct2density_calibration_file $CT2DENSITY_CALIBRATION_FILE

