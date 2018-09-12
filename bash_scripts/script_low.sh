#!/bin/bash

set -e
set -o pipefail
python2 ValidateStandardizedData.py low
bash /Users/sergiocamelo/Dropbox/Sergio-Joann/StandardizedData/instances/bash_daily.sh
bash /Users/sergiocamelo/Dropbox/Sergio-Joann/StandardizedData/instances/bash_spatial.sh
bash /Users/sergiocamelo/Dropbox/Sergio-Joann/StandardizedData/instances/bash_temporal.sh