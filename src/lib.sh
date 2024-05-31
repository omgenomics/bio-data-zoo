function log() {
    echo -n "$1... "
}

function validate() {
    if [[ "$1" != "" ]]; then
        echo "ok";
    else
        echo "failed";
        exit
    fi
}
