#!/bin/bash

echo "Pulling version from git"

GIT_COMMIT_HASH=$(git log -1 --format=%h-%ci | cut -d' ' -f1,2 --output-delimiter -)

GIT_STATUS=$(git status -s)
if [ "${GIT_STATUS}" != "" ]
then
    GIT_STATUS="D-"
fi

CMAKE_VERSION=$2

mkdir -p $( dirname $3 )

# @todo check if file needs to be overridden
cat $1 | sed "s/@CMAKE_VERSION@/$CMAKE_VERSION/g" \
       | sed "s/@GIT_COMMIT_HASH@/$GIT_COMMIT_HASH/g" \
       | sed "s/@GIT_STATUS@/$GIT_STATUS/g" \
       > $3