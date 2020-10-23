#!/usr/bin/env bash

pwd
ls -l ./levmar-2.6/*.so

if [[ -z "${DST_FOLDER}" ]]; then
	echo "DST_FOLDER must be defined"
	exit 1;
fi

gcc \
	-c \
	-fpic \
	-O3 \
	./gauss_2d.c \
	-I ./levmar-2.6 \
	-o "${DST_FOLDER}/gauss_2d.o"

gcc \
	-shared \
	-o "${DST_FOLDER}/liblmfits.so" \
	"${DST_FOLDER}/gauss_2d.o" \
	-L ./levmar-2.6 \
	-llevmar \
	-lm \
	-llapack \
	-lblas
