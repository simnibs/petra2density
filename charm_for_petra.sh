#!/bin/bash

PETRA2DENSITY_DIR="$(dirname "$(which "$0")")"

charm "$1" "$2" "$3" --usesettings $PETRA2DENSITY_DIR/charm.ini --forceqform
python $PETRA2DENSITY_DIR/petra2density.py "./m2m_$1"
