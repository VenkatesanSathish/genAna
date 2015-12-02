#!/bin/bash

# first file of sequences
ifnameA="${1}"
# second file of sequences
ifnameB="${2}"
# number of permutation to choose (sampling)
nbperm="${3}"
# prefix of the output file
prefofn="${4}"


PERM_CMDC="./permgen_w 1 2"

tdirn="/tmp/${USER}_$$"

## Create the temporary folder

mkdir -p "${tdirn}"

# sqn: sequence number
# sqc: pair of sub-sequence counter
declare -i sqc sqn tt ii

## Read the first file
sqn=1
sqc=0
while IFS="" read -r line; do
  if [[ "${line:0:1}" == ">" ]]; then
    sqc=$(($sqc+1))
    fname="${tdirn}/seq_${sqn}_${sqc}.fa"
    echo "${line}" > "${fname}"
  else
    echo "${line}" >> "${fname}"
  fi
done < <(cat "${ifnameA}")
tt="${sqc}"

## Read the second file
sqn=2
sqc=0
while IFS="" read -r line; do
  if [[ "${line:0:1}" == ">" ]]; then
    sqc=$(($sqc+1))
    fname="${tdirn}/seq_${sqn}_${sqc}.fa"
    echo "${line}" > "${fname}"
  else
    echo "${line}" >> "${fname}"
  fi
done < <(cat "${ifnameB}")

## First check-up (number of sequences in input files agrees)
if [[ "${tt}" -ne "${sqc}" ]]; then
  rm -rf "${tdirn}"
  echo "Error: the files ${ifnameA} and ${ifnameB} have different number of sequences"
  exit 1
fi

## Generate the permutations files (2 files of ${nbperm} lines)
kk=1
while [[ "${kk}" -le "${nbperm}" ]]; do
  ii=1
  while [[ "${ii}" -le "${sqc}" ]]; do
    ${PERM_CMDC} "${tdirn}/seq_1_${ii}.fa" "${tdirn}/seq_2_${ii}.fa" \
       "${tdirn}/pseq_k${kk}"
    ii=$(($ii+1))
  done
  kk=$(($kk+1))
done

## Second check-up (do we have repetition of permutations)
declare -i rn un
rn="$(cat ${tdirn}/pseq_k*_s1 | grep ">" | cut -f3 -d'|' | sort | wc -l)"
un="$(cat ${tdirn}/pseq_k*_s1 | grep ">" | cut -f3 -d'|' | sort -u | wc -l)"
if [[ "${rn}" -ne "${un}" ]]; then
  echo "Warning (not critical):"
  echo "  There were $((${rn}-${un})) repetitions !"
fi

## Copy the results into the current folder
kk=1
while [[ "${kk}" -le "${nbperm}" ]]; do
  ii=1
  while [[ "${ii}" -le "${sqc}" ]]; do
    cp -f "${tdirn}/pseq_k${kk}_s1" "${prefofn}_k${kk}_s1.fa"
    cp -f "${tdirn}/pseq_k${kk}_s2" "${prefofn}_k${kk}_s2.fa"
    ii=$(($ii+1))
  done
  kk=$(($kk+1))
done

## All done, clean-up
rm -rf "${tdirn}"

