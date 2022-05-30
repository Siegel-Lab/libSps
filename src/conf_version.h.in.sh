#!/bin/bash

GIT_COMMIT_HASH=$(git log -1 --format=%h)
GIT_COUNT=$(git rev-list --all --count)

GIT_STATUS=$(git status -s)
if [ "${GIT_STATUS}" != "" ]
then
    GIT_STATUS="-D"
fi

CMAKE_VERSION=$2

mkdir -p $(dirname $3)

cat $1 | sed "s/@CMAKE_VERSION@/$CMAKE_VERSION/g" \
       | sed "s/@GIT_COUNT@/$GIT_COUNT/g" \
       | sed "s/@GIT_COMMIT_HASH@/$GIT_COMMIT_HASH/g" \
       | sed "s/@GIT_STATUS@/$GIT_STATUS/g" \
       > $3