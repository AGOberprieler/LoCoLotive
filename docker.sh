#!/bin/bash

docker run -e LOCAL_USER_ID=$(id -u $USER) -v "$PWD":/wd ul90/locolotive "$@"
