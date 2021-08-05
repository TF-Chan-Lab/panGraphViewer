#!/bin/bash

set -e

# remove existing container
container_id=`docker ps -aqf "name=pangraph1"`
if [ ! -z "$container_id" ]
then
  docker stop $container_id > /dev/null
  docker rm $container_id > /dev/null
fi

docker run -d -P --name pangraph1 pangraph > /dev/null

port=$( docker port pangraph1 8000 | cut -d':' -f2 )
echo "Port for web is $port"
