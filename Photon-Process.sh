#!/bin/sh

if [ -z "$BASH_VERSION" ]; then
  # Not bash, so rely on sourcing from correct location
  if [ ! -f geant4.sh ]; then
    echo 'ERROR: Photon-Process.sh could NOT self-locate Photon-Process installation'
    echo 'This is most likely because you are using ksh, zsh or similar'
    echo 'To fix this issue, cd to the directory containing this script'
    echo 'and source it in that directory.'
    return 1
  fi
  envbindir=$(pwd)
else
  sourced_dir=$(dirname ${BASH_ARGV[0]})
  envbindir=$(cd $sourced_dir > /dev/null ; pwd)
fi

export PATH=$envbindir:$PATH
export Photon_Process_Data_Dir=$envbindir/../share/PhotonProcess/DataTables
