#!/bin/bash

function log() {
    echo -n "$1... "
}

function validate() {
    if [[ "$1" != "" ]]; then
        echo "done";
    else
        echo "failed";
        exit
    fi
}
