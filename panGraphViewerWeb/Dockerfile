FROM python:3.7-slim

RUN apt-get update && apt-get install build-essential zlib1g-dev libbz2-dev liblzma-dev -y

# set work directory
WORKDIR /usr/src/app

# install dependencies
RUN pip install --upgrade pip
COPY ./requirements.txt /usr/src/app
RUN pip install -r requirements.txt

# copy project
COPY . /usr/src/app

EXPOSE 8000:8000

CMD ./run.sh
