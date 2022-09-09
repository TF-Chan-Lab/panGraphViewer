#!/bin/bash

#mode=PROD

if [[ $mode == "PROD" ]]
then
   echo "Generating static files ..."
   python manage.py collectstatic --noinput

   echo "Starting server ...."
   python manage.py runserver 0.0.0.0:8000
else
   echo "Starting server ...."
   python manage.py runserver --insecure 0.0.0.0:8000
fi
