#!/bin/bash

USER_ID=${LOCAL_USER_ID:-9001}
useradd --shell /bin/bash -u $USER_ID -o -c "" -m user
export HOME=/home/user
chmod o+w /dev/stdout
chmod o+w /dev/stderr
exec /usr/local/bin/gosu user "$@"
