#!/usr/bin/env bash

if [[ -z "${DST_FOLDER}" ]]; then
	echo "DST_FOLDER must be defined"
	exit 1;
fi

if [[ -z "${LEVMAR_FOLDER}" ]]; then
	echo "LEVMAR_FOLDER must be defined"
	exit 1;
fi

gcc \
	-c \
	-fpic \
	-O3 \
	./gauss_2d.c \
	-I "${LEVMAR_FOLDER}" \
	-o "${DST_FOLDER}/gauss_2d.o"

gcc \
	-shared \
	-o "${DST_FOLDER}/liblmfits.so" \
	"${DST_FOLDER}/gauss_2d.o" \
	-L "${LEVMAR_FOLDER}" \
	-llevmar \
	-lm \
	-llapack \
	-lblas

#gcc \
#	-c \
#	./lmfits_test.c \
#	-O3 \
#	-I "${LEVMAR_FOLDER}" \
#	-o "${DST_FOLDER}/lmfits_test.o"
#
#gcc \
#	-c \
#	./gauss_2d.c \
#	-O3 \
#	-I "${LEVMAR_FOLDER}" \
#	-o "${DST_FOLDER}/_gauss_2d.o"
#
#gcc \
#	-o "${DST_FOLDER}/lmfits_test" \
#	"${DST_FOLDER}/_gauss_2d.o" \
#	"${DST_FOLDER}/lmfits_test.o" \
#	"${LEVMAR_FOLDER}/liblevmar.a" \
#	-lm \
#	-llapack \
#	-lblas
